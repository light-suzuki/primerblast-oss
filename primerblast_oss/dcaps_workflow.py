"""End-to-end dCAPS candidate construction and validation."""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .caps import (
    ENZYME_METADATA,
    apply_engineered_changes,
    caps_scan,
    dcaps_candidates,
    materialize_dcaps_primers,
    result_to_dict,
)
from .design import PrimerPair
from .pipeline import (
    _score_pair,
    resolve_genome_for_database,
    thermo_metadata,
)
from .specificity import SpecParams, pair_specificity


def _gc(sequence: str) -> float:
    sequence = sequence.upper()
    if not sequence:
        return 0.0
    return 100.0 * sum(base in "GC" for base in sequence) / len(sequence)


def _tm(sequence: str) -> float:
    """Primer3 Tm when available; deterministic Wallace fallback otherwise."""
    try:
        import primer3
        return round(float(primer3.calc_tm(sequence)), 1)
    except Exception:  # pragma: no cover - optional dependency
        sequence = sequence.upper()
        return float(2 * sum(base in "AT" for base in sequence)
                     + 4 * sum(base in "GC" for base in sequence))


def _candidate_pair(original: PrimerPair, candidate: Dict,
                    amplicon_start: int) -> PrimerPair:
    role = candidate["primer_role"]
    local_start = amplicon_start + int(candidate["primer_start"])
    local_end = amplicon_start + int(candidate["primer_end"])
    primer = candidate["primer_sequence"]
    if role == "F":
        left_start = local_start
        right_start = original.right_start
        forward, reverse = primer, original.reverse
        left_len, right_len = len(primer), original.right_len
        tm_f, tm_r = _tm(primer), original.tm_r
        gc_f, gc_r = _gc(primer), original.gc_r
    else:
        left_start = original.left_start
        right_start = local_end
        forward, reverse = original.forward, primer
        left_len, right_len = original.left_len, len(primer)
        tm_f, tm_r = original.tm_f, _tm(primer)
        gc_f, gc_r = original.gc_f, _gc(primer)
    return PrimerPair(
        index=original.index,
        template_id=original.template_id + "_dCAPS",
        forward=forward,
        reverse=reverse,
        left_start=left_start,
        left_len=left_len,
        right_start=right_start,
        right_len=right_len,
        product_size=right_start - left_start + 1,
        tm_f=tm_f,
        tm_r=tm_r,
        gc_f=gc_f,
        gc_r=gc_r,
        penalty=float(candidate.get("mismatches", 0)),
    )


def _engineered_products(template_sequence: str, pair: PrimerPair,
                         snp_local: int, ref_base: str, alt_base: str,
                         changes_global: Sequence[Dict]) -> Tuple[str, str]:
    low, high = pair.left_start, pair.right_start
    reference = list(template_sequence[low:high + 1].upper())
    alternate = list(reference)
    snp_relative = snp_local - low
    if not (0 <= snp_relative < len(reference)):
        return "", ""
    reference[snp_relative] = ref_base.upper()
    alternate[snp_relative] = alt_base.upper()
    relative_changes = []
    for change in changes_global:
        relative = int(change["position"]) - low
        if 0 <= relative < len(reference):
            relative_changes.append({**change, "position": relative})
    return (
        apply_engineered_changes("".join(reference), relative_changes),
        apply_engineered_changes("".join(alternate), relative_changes),
    )


def evaluate_dcaps_candidates(
    template,
    original_pair: PrimerPair,
    snp_local: int,
    alt_base: str,
    databases: Sequence[str],
    spec_params: Optional[SpecParams] = None,
    blastn_bin: Optional[str] = None,
    genomes_by_db: Optional[Mapping[str, object]] = None,
    thermo_params=None,
    thermo_gate: bool = True,
    dimer_params=None,
    variants: Optional[Sequence] = None,
    gel_min_gap: int = 25,
    max_engineered_mismatches: int = 2,
    max_candidates_to_screen: int = 12,
) -> Dict:
    """Build, digest, BLAST-screen and risk-rank dCAPS primer candidates."""
    from .assay import analyze_pair

    specificity = spec_params or SpecParams()
    amplicon_start = original_pair.left_start
    amplicon_end = original_pair.right_start
    amplicon = template.seq[amplicon_start:amplicon_end + 1].upper()
    snp_in_amplicon = snp_local - amplicon_start
    if not (0 <= snp_in_amplicon < len(amplicon)):
        return {
            "status": "snp_outside_amplicon",
            "n_frames": 0,
            "n_materialized": 0,
            "n_orderable": 0,
            "candidates": [],
        }
    ref_base = amplicon[snp_in_amplicon]
    frames = dcaps_candidates(
        amplicon,
        snp_in_amplicon,
        ref_base,
        alt_base,
        max_primer_mismatch=max_engineered_mismatches,
    )
    materialized = materialize_dcaps_primers(
        amplicon, snp_in_amplicon, frames)

    evaluated: List[Dict] = []
    for candidate in materialized[:max_candidates_to_screen]:
        candidate_pair = _candidate_pair(
            original_pair, candidate, amplicon_start)
        global_changes = [
            {**change, "position": amplicon_start + int(change["position"])}
            for change in candidate["engineered_changes"]
        ]
        reference_product, alternate_product = _engineered_products(
            template.seq,
            candidate_pair,
            snp_local,
            ref_base,
            alt_base,
            global_changes,
        )
        if not reference_product or not alternate_product:
            continue
        enzyme = ENZYME_METADATA[candidate["enzyme"]]
        digest_results = caps_scan(
            reference_product,
            alternate_product,
            enzymes={enzyme.name: enzyme},
            gel_min_gap=gel_min_gap,
        )
        digest = digest_results[0] if digest_results else None
        if digest is None:
            continue

        per_database = []
        for database in databases:
            genome, association = resolve_genome_for_database(
                database, databases, genomes_by_db=genomes_by_db)
            result = pair_specificity(
                candidate_pair.forward,
                candidate_pair.reverse,
                database,
                designed_size=candidate_pair.product_size,
                sp=specificity,
                blastn_bin=blastn_bin,
                genome=genome,
                thermo_params=thermo_params,
                thermo_gate=thermo_gate,
            )
            result.update(thermo_metadata(
                genome, thermo_params, thermo_gate, association))
            per_database.append(result)

        from . import dimers as dimer_module
        dimer = (
            dimer_module.analyze_pair(
                candidate_pair.forward, candidate_pair.reverse, dimer_params)
            if dimer_module.available() else None
        )
        _score_pair(candidate_pair, per_database, dimer)
        assay_summary = analyze_pair(
            candidate_pair,
            per_database,
            databases[0],
            template,
            variants or [],
            None,
            gel_min_gap=specificity.gel_min_gap_bp,
            dimer_params=dimer_params,
        )
        digest_dict = result_to_dict(digest)
        orderable = (
            digest.distinguishable
            and assay_summary.get("specific") is True
            and assay_summary.get("risk") in ("low", "medium")
            and assay_summary.get("search_complete_all_db") is True
        )
        if not digest.distinguishable:
            recommendation_status = "digest_not_resolvable"
        elif assay_summary.get("specificity_status") == "indeterminate":
            recommendation_status = "rerun_specificity_exhaustively"
        elif assay_summary.get("specific") is not True:
            recommendation_status = "not_specific"
        elif assay_summary.get("risk") == "high":
            recommendation_status = "high_risk"
        else:
            recommendation_status = "orderable"

        evaluated.append({
            "marker_type": "dCAPS",
            "enzyme": candidate["enzyme"],
            "recognition": candidate["recognition"],
            "recognition_orientation": candidate["orientation"],
            "present_in": candidate["present_in"],
            "engineered_mismatches": candidate["mismatches"],
            "engineered_changes_amplicon": candidate["engineered_changes"],
            "engineered_changes_template": global_changes,
            "modified_primer_role": candidate["primer_role"],
            "modified_primer": candidate["primer_sequence"],
            "forward": candidate_pair.forward,
            "reverse": candidate_pair.reverse,
            "product_size": candidate_pair.product_size,
            "tm_f": candidate_pair.tm_f,
            "tm_r": candidate_pair.tm_r,
            "gc_f": candidate_pair.gc_f,
            "gc_r": candidate_pair.gc_r,
            "digest": digest_dict,
            "specificity": assay_summary,
            "orderable": orderable,
            "recommendation_status": recommendation_status,
        })

    risk_order = {"low": 0, "medium": 1, "high": 2}
    evaluated.sort(key=lambda candidate: (
        not candidate["orderable"],
        risk_order.get(candidate["specificity"].get("risk"), 3),
        candidate["engineered_mismatches"],
        abs(candidate["tm_f"] - candidate["tm_r"]),
        -candidate["digest"].get("min_gel_gap", 0),
    ))
    return {
        "status": "candidates_found" if evaluated else "no_valid_candidate",
        "n_frames": len(frames),
        "n_materialized": len(materialized),
        "n_screened": min(len(materialized), max_candidates_to_screen),
        "n_evaluated": len(evaluated),
        "n_orderable": sum(candidate["orderable"] for candidate in evaluated),
        "best": evaluated[0] if evaluated else None,
        "candidates": evaluated,
    }
