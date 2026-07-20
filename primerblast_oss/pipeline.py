"""End-to-end pipeline: design -> multi-database specificity -> scoring."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from .design import DesignParams, PrimerPair, design_primers
from .specificity import (
    SEARCH_COMPLETE,
    SpecParams,
    combine_search_completeness,
    pair_specificity,
)


@dataclass
class PipelineResult:
    template_id: str
    template_len: int
    pairs: List[PrimerPair]
    primer3_explain: str
    databases: List[str]
    params: Dict = field(default_factory=dict)


def _score_pair(pair: PrimerPair, per_db: Sequence[Dict],
                dimer: Optional[Dict] = None) -> None:
    """Rank a pair while separating observed products from search completeness."""
    total_off = sum(result["n_off_target"] for result in per_db)
    total_on = sum(result["n_on_target"] for result in per_db)
    total_comigrating = sum(result.get("n_comigrating", 0) for result in per_db)
    total_distinguishable_off = total_off - total_comigrating

    completeness = combine_search_completeness([
        result.get("search_completeness", SEARCH_COMPLETE) for result in per_db
    ])
    incomplete_databases = [
        result["db"] for result in per_db
        if result.get("search_completeness", SEARCH_COMPLETE) != SEARCH_COMPLETE
    ]
    search_complete_all = completeness == SEARCH_COMPLETE and len(per_db) > 0
    observed_specific_all = (
        all(result.get("specific_observed", result.get("specific") is True)
            for result in per_db)
        and len(per_db) > 0
    )
    explicit_non_specific = any(result.get("specific") is False for result in per_db)
    specific_all = (
        all(result.get("specific") is True for result in per_db)
        and len(per_db) > 0
    )
    gel_clean = (
        all(result.get("gel_distinguishable", True) for result in per_db)
        and len(per_db) > 0
    )

    if explicit_non_specific:
        specificity_status = "non_specific"
    elif not search_complete_all and observed_specific_all:
        specificity_status = "indeterminate"
    elif specific_all:
        specificity_status = "specific"
    else:
        specificity_status = "non_specific"

    score = 100.0
    score -= 25.0 * total_comigrating
    score -= 4.0 * total_distinguishable_off
    if total_on == 0:
        score -= 40.0
    if specificity_status == "indeterminate":
        score -= 15.0
    score -= min(10.0, abs(pair.tm_f - pair.tm_r) * 2.0)
    for gc in (pair.gc_f, pair.gc_r):
        if gc < 30.0 or gc > 70.0:
            score -= 3.0
    score -= min(10.0, pair.self_end_th / 5.0)

    dimer_concern = bool(dimer and dimer.get("n_concerning", 0) > 0)
    if dimer_concern:
        score -= 12.0
    score = max(0.0, min(100.0, score))

    # I is an explicit non-orderable/needs-rerun rank. It avoids disguising an
    # incomplete clean search as B while still sorting ahead of genuinely poor
    # products by score when users inspect the full output.
    if specificity_status == "indeterminate":
        rank = "I"
    elif dimer_concern:
        rank = "C" if (total_comigrating == 0 and total_on > 0) else "D"
    elif specific_all and score >= 85.0:
        rank = "A"
    elif total_comigrating == 0 and total_on > 0:
        rank = "B"
    elif total_comigrating <= 1:
        rank = "C"
    else:
        rank = "D"

    pair.specificity = {
        "per_db": list(per_db),
        "total_off_target": total_off,
        "total_on_target": total_on,
        "total_comigrating": total_comigrating,
        "specific_all_db": specific_all,
        "specific_observed_all_db": observed_specific_all,
        "specificity_status": specificity_status,
        "search_completeness": completeness,
        "search_complete_all_db": search_complete_all,
        "incomplete_databases": incomplete_databases,
        "gel_distinguishable": gel_clean,
        "score": round(score, 1),
        "rank": rank,
        "dimers": ({
            "worst_dg": dimer["worst_dg"],
            "cross_dimer_dg": dimer["cross_dimer_dg"],
            "n_concerning": dimer["n_concerning"],
            "ok": dimer["ok"],
            "concerning": [
                {"kind": structure.kind, "a": structure.a, "b": structure.b,
                 "tm": structure.tm, "dg": structure.dg}
                for structure in dimer["structures"] if structure.concerning
            ],
        } if dimer else None),
    }


def run_pipeline(
    template_id: str,
    sequence: str,
    databases: Sequence[str],
    design_params: Optional[DesignParams] = None,
    spec_params: Optional[SpecParams] = None,
    primer3_bin: Optional[str] = None,
    blastn_bin: Optional[str] = None,
    size_tolerance: int = 10,
    genome=None,
    thermo_params=None,
    thermo_gate: bool = True,
    dimer_params=None,
) -> PipelineResult:
    design_params = design_params or DesignParams()
    spec_params = spec_params or SpecParams()
    pairs, explain = design_primers(template_id, sequence, design_params, primer3_bin)

    from . import dimers as dimer_module
    for pair in pairs:
        per_db: List[Dict] = []
        for database in databases:
            per_db.append(pair_specificity(
                pair.forward,
                pair.reverse,
                database,
                designed_size=pair.product_size,
                sp=spec_params,
                blastn_bin=blastn_bin,
                size_tolerance=size_tolerance,
                genome=genome,
                thermo_params=thermo_params,
                thermo_gate=thermo_gate,
            ))
        dimer = (
            dimer_module.analyze_pair(pair.forward, pair.reverse, dimer_params)
            if dimer_module.available() else None
        )
        _score_pair(pair, per_db, dimer)

    pairs.sort(key=lambda pair: (
        pair.specificity.get("rank") == "I",
        -pair.specificity.get("score", 0.0),
        pair.specificity.get("total_off_target", 999),
        pair.penalty,
    ))

    from .design import clean_sequence
    return PipelineResult(
        template_id=template_id,
        template_len=len(clean_sequence(sequence)),
        pairs=pairs,
        primer3_explain=explain,
        databases=list(databases),
        params={
            "design": design_params.__dict__,
            "specificity": spec_params.__dict__,
            "size_tolerance": size_tolerance,
        },
    )
