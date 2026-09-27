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
    if (len(ref_base) != 1 or ref_base not in _BASES
            or len(alt_base) != 1 or alt_base not in _BASES):
        raise ValueError("ref and alt must be single A/C/G/T bases")
    if sequence[snp_index] != ref_base:
        raise ValueError("ref allele does not match the reference sequence")

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
                    "outer_left_start": parent_pair.left_start,
                    "outer_right_start": parent_pair.right_start,
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


def _amplicon_dict(amplicon) -> Dict:
    return {
        "subject": amplicon.subject,
        "start": amplicon.start,
        "end": amplicon.end,
        "size": amplicon.size,
        "orientation": amplicon.orientation,
        "fwd_mismatch": amplicon.fwd_mismatch,
        "rev_mismatch": amplicon.rev_mismatch,
    }


def _screen_primer_pool(
    primers: Dict[str, str],
    databases: Sequence[str],
    *,
    spec_params=None,
    blastn_bin=None,
    genomes_by_db=None,
    thermo_params=None,
    thermo_gate: bool = True,
    expected_sizes: Sequence[int] = (),
    size_tolerance: int = 10,
    expected_products: Sequence[Dict] = (),
    design_db: Optional[str] = None,
    required_orientation: Optional[str] = None,
    dimer_params=None,
) -> Dict:
    from .pipeline import resolve_genome_for_database, thermo_metadata
    from .specificity import SpecParams, in_silico_pcr

    sp = spec_params or SpecParams()
    from . import dimers
    from dataclasses import asdict
    structures = []
    if dimers.available():
        for name, sequence in primers.items():
            structures.extend([dimers.hairpin(name, sequence, dimer_params),
                               dimers.self_dimer(name, sequence, dimer_params)])
        cross = dimers.analyze_multiplex(list(primers.items()), dimer_params) or {}
        structures.extend(cross.get("concerning", []))
    structures = [asdict(structure) for structure in structures if structure is not None]
    structure_summary = {
        "status": "evaluated" if dimers.available() else "not_evaluated",
        "n_concerning": sum(bool(structure["concerning"]) for structure in structures),
        "structures": structures,
    }
    per_db = []
    for database in databases:
        genome, association = resolve_genome_for_database(
            database, databases, genomes_by_db=genomes_by_db)
        result = in_silico_pcr(
            primers,
            database,
            sp=sp,
            blastn_bin=blastn_bin,
            genome=genome,
            thermo_params=thermo_params,
            thermo_gate=thermo_gate,
        )
        result.update(thermo_metadata(
            genome, thermo_params, thermo_gate, association,
            result.get("thermo_site_stats")))
        products = [_amplicon_dict(product) for product in result["products"]]
        expected_matches = [
            product for product in products
            if any(abs(product["size"] - size) <= size_tolerance
                   for size in expected_sizes)
        ]
        anchored = bool(expected_products) and database == design_db
        intended = [
            product for product in products
            if anchored and any(
                all(product[key] == expected[key]
                    for key in ("subject", "start", "end", "orientation"))
                for expected in expected_products)
        ]
        unexpected = ([product for product in products if product not in intended]
                      if anchored else [])
        required_observed = required_orientation is None or any(
            product["orientation"] == required_orientation for product in intended)
        per_db.append({
            "db": database,
            "n_products": len(products),
            "products": products,
            "expected_size_matches": expected_matches,
            "intended_products": intended,
            "classification": "genomic_anchor" if anchored else "unverified_size_only",
            "required_product_observed": required_observed,
            "unique_intended_products": len(intended) == len({
                (p["subject"], p["start"], p["end"], p["orientation"]) for p in intended}),
            "n_unexpected_products": len(unexpected) if anchored else None,
            "unexpected_products": unexpected,
            "unverified_products": products if not anchored else [],
            "search_completeness": result.get("search_completeness"),
            "search_complete": result.get("search_complete"),
            "completeness_recommendation": result.get(
                "completeness_recommendation"),
            "thermo_status": result.get("thermo_status"),
        })
    return {
        "status": "screened",
        "primer_structures": structure_summary,
        "search_complete_all_db": (
            bool(per_db)
            and all(view.get("search_complete") is True for view in per_db)
        ),
        "max_unexpected_products": max(
            (view["n_unexpected_products"] for view in per_db
             if view["n_unexpected_products"] is not None), default=None),
        "genome_screen_acceptable": bool(per_db) and all(
            view["classification"] == "genomic_anchor"
            and view["search_complete"] is True
            and view["n_unexpected_products"] == 0
            and view["required_product_observed"]
            and view["unique_intended_products"]
            for view in per_db),
        "per_db": per_db,
    }


def _expected_product(template, start, end, left_name, right_name) -> Dict:
    g1, g2 = template.to_genomic(start), template.to_genomic(end)
    if g1 > g2:
        left_name, right_name = right_name, left_name
    return {"subject": template.region.chrom, "start": min(g1, g2),
            "end": max(g1, g2), "orientation": left_name + "/" + right_name}


def screen_aspcr(
    aspcr: Dict,
    databases: Sequence[str],
    *,
    spec_params=None,
    blastn_bin=None,
    genomes_by_db=None,
    thermo_params=None,
    thermo_gate: bool = True,
    size_tolerance: int = 10,
    template=None,
    design_db: Optional[str] = None,
    dimer_params=None,
) -> Dict:
    """BLAST-screen the best classical and tetra-ARMS primer sets."""
    output = dict(aspcr)
    screens = {}
    for allele, key in (
        ("ref", "best_classical_ref"),
        ("alt", "best_classical_alt"),
    ):
        candidate = output.get(key)
        if not candidate:
            continue
        expected = []
        if template is not None:
            if candidate["role"] == "F":
                start = candidate["primer_5p"]
                end = start + candidate["product_size"] - 1
                left_name, right_name = allele + "_AS", "common_R"
            else:
                end = candidate["primer_5p"]
                start = end - candidate["product_size"] + 1
                left_name, right_name = "common_F", allele + "_AS"
            expected = [_expected_product(template, start, end, left_name, right_name)]
        screens["classical_" + allele] = _screen_primer_pool(
            {
                "%s_AS" % allele: candidate["primer"],
                "common_%s" % candidate["common_role"]: candidate["common_primer"],
            },
            databases,
            spec_params=spec_params,
            blastn_bin=blastn_bin,
            genomes_by_db=genomes_by_db,
            thermo_params=thermo_params,
            thermo_gate=thermo_gate,
            expected_sizes=[candidate["product_size"]],
            size_tolerance=size_tolerance,
            expected_products=expected,
            design_db=design_db,
            required_orientation=expected[0]["orientation"] if expected and allele == "ref" else None,
            dimer_params=dimer_params,
        )

    tetra = output.get("best_tetra")
    if tetra:
        expected = []
        if template is not None:
            left, right = tetra["outer_left_start"], tetra["outer_right_start"]
            expected.append(_expected_product(template, left, right, "outer_F", "outer_R"))
            for allele in ("ref", "alt"):
                inner = tetra[allele + "_inner"]
                name = allele + "_inner_" + inner["role"]
                if inner["role"] == "F":
                    expected.append(_expected_product(template, inner["primer_5p"], right, name, "outer_R"))
                else:
                    expected.append(_expected_product(template, left, inner["primer_5p"], "outer_F", name))
        screens["tetra_arms"] = _screen_primer_pool(
            {
                "outer_F": tetra["outer_forward"],
                "outer_R": tetra["outer_reverse"],
                "ref_inner_%s" % tetra["ref_inner"]["role"]:
                    tetra["ref_inner"]["primer"],
                "alt_inner_%s" % tetra["alt_inner"]["role"]:
                    tetra["alt_inner"]["primer"],
            },
            databases,
            spec_params=spec_params,
            blastn_bin=blastn_bin,
            genomes_by_db=genomes_by_db,
            thermo_params=thermo_params,
            thermo_gate=thermo_gate,
            expected_sizes=[
                tetra["control_product_size"],
                tetra["ref_product_size"],
                tetra["alt_product_size"],
            ],
            size_tolerance=size_tolerance,
            expected_products=expected,
            design_db=design_db,
            required_orientation=expected[0]["orientation"] if expected else None,
            dimer_params=dimer_params,
        )

    output["specificity_screen"] = {
        "status": "screened" if screens else "no_candidate",
        "sets": screens,
        "search_complete_all_sets": (
            bool(screens)
            and all(screen.get("search_complete_all_db") is True
                    for screen in screens.values())
        ),
    }
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
