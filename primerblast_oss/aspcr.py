"""Allele-specific PCR / ARMS and tetra-primer ARMS design helpers."""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

from .genome import revcomp
from .gel import genotype_pattern_discrimination, recommend_gel_percent, recommend_ladder


_BASES = "ACGT"


def _gc(sequence: str) -> float:
    sequence = sequence.upper()
    return (
        100.0 * sum(base in "GC" for base in sequence) / len(sequence)
        if sequence else 0.0
    )


def _tm(sequence: str) -> float:
    try:
        import primer3
        return round(float(primer3.calc_tm(sequence)), 1)
    except Exception:
        sequence = sequence.upper()
        return float(
            2 * sum(base in "AT" for base in sequence)
            + 4 * sum(base in "GC" for base in sequence)
        )


def _allele_window(
    sequence: str,
    snp_index: int,
    allele: str,
    role: str,
    length: int,
) -> Optional[Tuple[str, int, int]]:
    """Return exact-match allele-specific primer before deliberate mismatch."""
    sequence = sequence.upper()
    allele = allele.upper()
    if role == "F":
        start = snp_index - length + 1
        end = snp_index
        if start < 0 or end >= len(sequence):
            return None
        target = list(sequence[start:end + 1])
        target[-1] = allele
        return "".join(target), start, end
    if role == "R":
        start = snp_index
        end = snp_index + length - 1
        if start < 0 or end >= len(sequence):
            return None
        target = list(sequence[start:end + 1])
        target[0] = allele
        return revcomp("".join(target)), start, end
    raise ValueError("role must be F or R")


def _target_for_allele(
    sequence: str,
    snp_index: int,
    allele: str,
    role: str,
    length: int,
) -> Optional[str]:
    result = _allele_window(sequence, snp_index, allele, role, length)
    return None if result is None else result[0]


def _mismatch_profile(primer: str, target: str) -> Dict:
    if len(primer) != len(target):
        raise ValueError("primer and allele target must have equal length")
    mismatches = [
        index for index, (p, t) in enumerate(zip(primer, target)) if p != t
    ]
    n = len(primer)
    return {
        "total": len(mismatches),
        "terminal_match": primer[-1] == target[-1],
        "terminal_mismatch": primer[-1] != target[-1],
        "mismatch_positions_from_3prime": [
            n - index for index in mismatches
        ],
        "last3": sum(index >= n - 3 for index in mismatches),
        "last5": sum(index >= n - 5 for index in mismatches),
    }


def _discrimination_score(
    intended: Dict,
    non_target: Dict,
    tm: float,
    common_tm: Optional[float],
    deliberate_position: Optional[int],
) -> float:
    """Interpretable ranking score, not a probability of allele specificity."""
    score = 50.0
    if non_target["terminal_mismatch"]:
        score += 25.0
    score += min(15.0, 7.5 * max(
        0, non_target["last3"] - intended["last3"]))
    score += min(10.0, 5.0 * max(
        0, non_target["last5"] - intended["last5"]))
    if deliberate_position == 2:
        score += 5.0
    elif deliberate_position == 3:
        score += 3.0
    score -= 12.0 * max(0, intended["last3"] - 1)
    if common_tm is not None:
        score -= min(15.0, 2.5 * abs(tm - common_tm))
    return round(max(0.0, min(100.0, score)), 1)


def generate_as_primers(
    sequence: str,
    snp_index: int,
    ref_base: str,
    alt_base: str,
    *,
    roles: Sequence[str] = ("F", "R"),
    min_length: int = 18,
    max_length: int = 26,
    opt_length: int = 20,
    deliberate_positions: Sequence[Optional[int]] = (2, 3, None),
    opt_tm: float = 60.0,
    common_tm_by_role: Optional[Dict[str, float]] = None,
) -> List[Dict]:
    """Generate ref/alt allele-specific primers terminating at the SNP.

    Deliberate mismatch positions are counted from the 3' end (2=penultimate,
    3=antepenultimate). All alternate bases are enumerated and ranked.
    """
    sequence = sequence.upper()
    ref_base = ref_base.upper()
    alt_base = alt_base.upper()
    if ref_base == alt_base:
        raise ValueError("ref and alt alleles must differ")
    if not (0 <= snp_index < len(sequence)):
        raise ValueError("SNP index outside sequence")
    if sequence[snp_index] not in _BASES:
        raise ValueError("reference sequence at SNP is not A/C/G/T")

    candidates: List[Dict] = []
    for role in roles:
        common_role = "R" if role == "F" else "F"
        common_tm = (common_tm_by_role or {}).get(common_role)
        for allele_name, allele, other in (
            ("ref", ref_base, alt_base),
            ("alt", alt_base, ref_base),
        ):
            for length in range(min_length, max_length + 1):
                built = _allele_window(
                    sequence, snp_index, allele, role, length)
                if built is None:
                    continue
                exact_primer, start, end = built
                other_target = _target_for_allele(
                    sequence, snp_index, other, role, length)
                if other_target is None:
                    continue
                for position in deliberate_positions:
                    if position is None:
                        variants = [(None, exact_primer)]
                    else:
                        if position < 2 or position > length:
                            continue
                        index = length - position
                        original = exact_primer[index]
                        variants = []
                        for base in _BASES:
                            if base == original:
                                continue
                            primer = (
                                exact_primer[:index]
                                + base
                                + exact_primer[index + 1:]
                            )
                            variants.append(({
                                "position_from_3prime": position,
                                "primer_index": index,
                                "from": original,
                                "to": base,
                            }, primer))
                    for deliberate, primer in variants:
                        intended_target = _target_for_allele(
                            sequence, snp_index, allele, role, length)
                        intended = _mismatch_profile(primer, intended_target)
                        non_target = _mismatch_profile(primer, other_target)
                        tm = _tm(primer)
                        score = _discrimination_score(
                            intended,
                            non_target,
                            tm,
                            common_tm,
                            None if deliberate is None
                            else deliberate["position_from_3prime"],
                        )
                        # Mild preference for usual primer length/Tm, without
                        # making those hard filters.
                        score -= min(8.0, 0.8 * abs(length - opt_length))
                        score -= min(8.0, 0.8 * abs(tm - opt_tm))
                        if role == "F":
                            primer_5p = start
                            primer_3p = end
                        else:
                            primer_5p = end
                            primer_3p = start
                        candidates.append({
                            "allele": allele_name,
                            "allele_base": allele,
                            "non_target_base": other,
                            "role": role,
                            "primer": primer,
                            "length": length,
                            "primer_start": start,
                            "primer_end": end,
                            "primer_5p": primer_5p,
                            "primer_3p": primer_3p,
                            "snp_index": snp_index,
                            "deliberate_mismatch": deliberate,
                            "tm": round(tm, 1),
                            "gc": round(_gc(primer), 1),
                            "intended_profile": intended,
                            "non_target_profile": non_target,
                            "discrimination_score": round(
                                max(0.0, min(100.0, score)), 1),
                            "evidence_note": (
                                "3'-terminal mismatch discrimination is "
                                "polymerase/condition dependent; wet-lab "
                                "optimization may be required"
                            ),
                        })
    candidates.sort(key=lambda item: (
        -item["discrimination_score"],
        abs(item["tm"] - opt_tm),
        abs(item["length"] - opt_length),
        item["allele"],
        item["role"],
    ))
    return candidates


def classical_assays_from_pair(
    sequence: str,
    parent_pair,
    snp_index: int,
    ref_base: str,
    alt_base: str,
    *,
    top_per_allele: int = 6,
    min_length: int = 18,
    max_length: int = 26,
    opt_length: int = 20,
    deliberate_positions: Sequence[Optional[int]] = (2, 3, None),
    min_product: int = 50,
    max_product: int = 1500,
) -> List[Dict]:
    """Pair AS primers with the opposite common primer from a parent pair."""
    candidates = generate_as_primers(
        sequence,
        snp_index,
        ref_base,
        alt_base,
        min_length=min_length,
        max_length=max_length,
        opt_length=opt_length,
        deliberate_positions=deliberate_positions,
        common_tm_by_role={"F": parent_pair.tm_f, "R": parent_pair.tm_r},
    )
    assays = []
    for candidate in candidates:
        if candidate["role"] == "F":
            product_size = parent_pair.right_start - candidate["primer_5p"] + 1
            common_role = "R"
            common_primer = parent_pair.reverse
        else:
            product_size = candidate["primer_5p"] - parent_pair.left_start + 1
            common_role = "F"
            common_primer = parent_pair.forward
        if not (min_product <= product_size <= max_product):
            continue
        assays.append({
            **candidate,
            "common_role": common_role,
            "common_primer": common_primer,
            "product_size": product_size,
            "reaction": "%s-specific + common-%s" % (
                candidate["allele"], common_role),
        })

    selected = []
    for allele in ("ref", "alt"):
        selected.extend([
            item for item in assays if item["allele"] == allele
        ][:top_per_allele])
    selected.sort(key=lambda item: (
        item["allele"], -item["discrimination_score"]))
    return selected


def _tetra_gel_score(pattern: Dict) -> float:
    gaps = []
    for comparison in pattern["comparisons"].values():
        gap = comparison.get("best_unique_log10_gap")
        if gap is not None:
            gaps.append(float(gap))
    if not gaps:
        return 0.0
    return round(min(100.0, 1000.0 * min(gaps)), 1)


def tetra_arms_from_pair(
    sequence: str,
    parent_pair,
    snp_index: int,
    ref_base: str,
    alt_base: str,
    *,
    candidates_per_orientation: int = 5,
    min_length: int = 18,
    max_length: int = 26,
    opt_length: int = 20,
    deliberate_positions: Sequence[Optional[int]] = (2, 3),
    min_product: int = 50,
) -> List[Dict]:
    """Build single-tube tetra-ARMS candidates from a parent outer pair."""
    all_candidates = generate_as_primers(
        sequence,
        snp_index,
        ref_base,
        alt_base,
        min_length=min_length,
        max_length=max_length,
        opt_length=opt_length,
        deliberate_positions=deliberate_positions,
        common_tm_by_role={"F": parent_pair.tm_f, "R": parent_pair.tm_r},
    )

    by_key: Dict[Tuple[str, str], List[Dict]] = {}
    for candidate in all_candidates:
        by_key.setdefault(
            (candidate["allele"], candidate["role"]), []).append(candidate)
    for key in by_key:
        by_key[key] = by_key[key][:candidates_per_orientation]

    output = []
    for ref_role, alt_role in (("F", "R"), ("R", "F")):
        for ref_candidate in by_key.get(("ref", ref_role), []):
            for alt_candidate in by_key.get(("alt", alt_role), []):
                if ref_role == "F":
                    ref_size = (
                        parent_pair.right_start
                        - ref_candidate["primer_5p"] + 1)
                    alt_size = (
                        alt_candidate["primer_5p"]
                        - parent_pair.left_start + 1)
                else:
                    ref_size = (
                        ref_candidate["primer_5p"]
                        - parent_pair.left_start + 1)
                    alt_size = (
                        parent_pair.right_start
                        - alt_candidate["primer_5p"] + 1)
                control_size = parent_pair.product_size
                if min(ref_size, alt_size, control_size) < min_product:
                    continue
                # Allele-specific bands should be smaller than outer control.
                if ref_size >= control_size or alt_size >= control_size:
                    continue

                pattern = genotype_pattern_discrimination(
                    [control_size, ref_size],
                    [control_size, alt_size],
                    min_log10_gap=0.04,
                )
                ladder_name, ladder_bands = recommend_ladder(
                    [control_size, ref_size, alt_size])
                gel_percent = recommend_gel_percent(
                    [control_size, ref_size, alt_size])
                gel_score = _tetra_gel_score(pattern)
                combined = round(
                    0.35 * ref_candidate["discrimination_score"]
                    + 0.35 * alt_candidate["discrimination_score"]
                    + 0.30 * gel_score,
                    1,
                )
                output.append({
                    "mode": "tetra-ARMS",
                    "outer_forward": parent_pair.forward,
                    "outer_reverse": parent_pair.reverse,
                    "control_product_size": control_size,
                    "ref_inner": ref_candidate,
                    "alt_inner": alt_candidate,
                    "ref_product_size": ref_size,
                    "alt_product_size": alt_size,
                    "genotype_bands": {
                        "AA_ref": sorted(
                            {control_size, ref_size}, reverse=True),
                        "AB": sorted(
                            {control_size, ref_size, alt_size}, reverse=True),
                        "BB_alt": sorted(
                            {control_size, alt_size}, reverse=True),
                    },
                    "gel_discrimination": pattern,
                    "gel_scorable": pattern["distinguishable"],
                    "gel_score": gel_score,
                    "ladder": ladder_name,
                    "ladder_bands": ladder_bands,
                    "gel_percent": gel_percent,
                    "score": combined,
                    "evidence_note": (
                        "Tetra-ARMS band predictions assume allele-specific "
                        "extension behaves as designed; wet optimization may "
                        "still be required"
                    ),
                })
    output.sort(key=lambda item: (
        not item["gel_scorable"],
        -item["score"],
        abs(item["ref_product_size"] - item["alt_product_size"]),
    ))
    return output


def build_aspcr(
    sequence: str,
    parent_pair,
    snp_index: int,
    alt_base: str,
    *,
    max_classical_per_allele: int = 4,
    max_tetra: int = 6,
    min_length: int = 18,
    max_length: int = 26,
    opt_length: int = 20,
) -> Dict:
    """Build classical AS-PCR and tetra-ARMS options for one parent pair."""
    sequence = sequence.upper()
    if not (0 <= snp_index < len(sequence)):
        return {
            "status": "snp_outside_template",
            "classical": [],
            "tetra_arms": [],
        }
    ref_base = sequence[snp_index]
    if ref_base not in _BASES or alt_base.upper() not in _BASES:
        return {
            "status": "non_snv",
            "classical": [],
            "tetra_arms": [],
        }
    classical = classical_assays_from_pair(
        sequence,
        parent_pair,
        snp_index,
        ref_base,
        alt_base,
        top_per_allele=max_classical_per_allele,
        min_length=min_length,
        max_length=max_length,
        opt_length=opt_length,
    )
    tetra = tetra_arms_from_pair(
        sequence,
        parent_pair,
        snp_index,
        ref_base,
        alt_base,
        min_length=min_length,
        max_length=max_length,
        opt_length=opt_length,
    )[:max_tetra]
    return {
        "status": "candidates_found" if classical or tetra else "no_candidate",
        "ref_base": ref_base,
        "alt_base": alt_base.upper(),
        "snp_index": snp_index,
        "classical": classical,
        "tetra_arms": tetra,
        "best_classical_ref": next(
            (item for item in classical if item["allele"] == "ref"), None),
        "best_classical_alt": next(
            (item for item in classical if item["allele"] == "alt"), None),
        "best_tetra": tetra[0] if tetra else None,
        "warning": (
            "Allele-specific PCR selectivity is polymerase and condition "
            "dependent; 3'-terminal mismatch logic is a design heuristic and "
            "requires empirical validation."
        ),
    }
