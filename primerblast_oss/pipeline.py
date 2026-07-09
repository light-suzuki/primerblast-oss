"""End-to-end pipeline: design -> multi-database specificity -> scoring."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

from .design import DesignParams, PrimerPair, design_primers
from .specificity import SpecParams, pair_specificity


@dataclass
class PipelineResult:
    template_id: str
    template_len: int
    pairs: List[PrimerPair]
    primer3_explain: str
    databases: List[str]
    params: Dict = field(default_factory=dict)


def _score_pair(pair: PrimerPair, per_db: Sequence[Dict], dimer: Optional[Dict] = None) -> None:
    """Rank a pair by specificity across all databases (higher = better).

    Off-target products that co-migrate with the intended one (similar size)
    are heavily penalized; off-targets far enough in size to be resolved on a
    gel are only a minor penalty, matching how such pairs are used in practice.
    A concerning primer-dimer / hairpin (from primer3-py, when available) adds a
    penalty and caps the rank at C.
    """
    total_off = sum(d["n_off_target"] for d in per_db)
    total_on = sum(d["n_on_target"] for d in per_db)
    total_comig = sum(d.get("n_comigrating", 0) for d in per_db)
    total_distinguishable_off = total_off - total_comig
    specific_all = all(d["specific"] for d in per_db) and len(per_db) > 0
    gel_clean = all(d.get("gel_distinguishable", True) for d in per_db) and len(per_db) > 0

    score = 100.0
    score -= 25.0 * total_comig                  # co-migrating off-targets: bad
    score -= 4.0 * total_distinguishable_off     # resolvable off-targets: minor
    if total_on == 0:
        score -= 40.0                            # intended product not recovered
    score -= min(10.0, abs(pair.tm_f - pair.tm_r) * 2.0)
    for gc in (pair.gc_f, pair.gc_r):
        if gc < 30.0 or gc > 70.0:
            score -= 3.0
    score -= min(10.0, pair.self_end_th / 5.0)   # penalize strong 3' dimers
    dimer_concern = bool(dimer and dimer.get("n_concerning", 0) > 0)
    if dimer_concern:
        score -= 12.0
    score = max(0.0, min(100.0, score))

    if dimer_concern:
        rank = "C" if (total_comig == 0 and total_on > 0) else "D"
    elif specific_all and score >= 85.0:
        rank = "A"                               # single product everywhere
    elif total_comig == 0 and total_on > 0:
        rank = "B"                               # extra products, but gel-resolvable
    elif total_comig <= 1:
        rank = "C"
    else:
        rank = "D"

    pair.specificity = {
        "per_db": list(per_db),
        "total_off_target": total_off,
        "total_on_target": total_on,
        "total_comigrating": total_comig,
        "specific_all_db": specific_all,
        "gel_distinguishable": gel_clean,
        "score": round(score, 1),
        "rank": rank,
        "dimers": ({
            "worst_dg": dimer["worst_dg"], "cross_dimer_dg": dimer["cross_dimer_dg"],
            "n_concerning": dimer["n_concerning"], "ok": dimer["ok"],
            "concerning": [{"kind": s.kind, "a": s.a, "b": s.b, "tm": s.tm, "dg": s.dg}
                           for s in dimer["structures"] if s.concerning],
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

    from . import dimers as _dimers
    for pair in pairs:
        per_db: List[Dict] = []
        for db in databases:
            res = pair_specificity(
                pair.forward, pair.reverse, db,
                designed_size=pair.product_size,
                sp=spec_params, blastn_bin=blastn_bin,
                size_tolerance=size_tolerance,
                genome=genome, thermo_params=thermo_params, thermo_gate=thermo_gate,
            )
            per_db.append(res)
        dimer = (_dimers.analyze_pair(pair.forward, pair.reverse, dimer_params)
                 if _dimers.available() else None)
        _score_pair(pair, per_db, dimer)

    # best (specific + high score) first
    pairs.sort(
        key=lambda p: (
            -p.specificity.get("score", 0.0),
            p.specificity.get("total_off_target", 999),
            p.penalty,
        )
    )

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
