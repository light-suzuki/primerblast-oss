"""Experimenter-facing risk assessment."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class RiskAssessment:
    level: str
    score: float
    reasons: List[str] = field(default_factory=list)
    signals: Dict = field(default_factory=dict)


def assess_risk(
    n_comigrating_offtarget: int = 0,
    n_ff: int = 0,
    n_rr: int = 0,
    n_fr_offtarget: int = 0,
    offtarget_min_tp5: Optional[int] = None,
    snp_in_primer: bool = False,
    snp_in_primer_3prime: bool = False,
    tm_diff: float = 0.0,
    gc_f: Optional[float] = None,
    gc_r: Optional[float] = None,
    gel_distinguishable: bool = True,
    conserved_fraction: Optional[float] = None,
    repeat_overlap: bool = False,
    dimer_concern: bool = False,
    cross_dimer_dg: Optional[float] = None,
    intended_status: str = "unique",
    search_completeness: str = "complete",
) -> RiskAssessment:
    score = 100.0
    reasons: List[str] = []

    if intended_status == "missing":
        score -= 50.0
        reasons.append("designed intended amplicon was not recovered at the expected coordinates")
    elif intended_status == "ambiguous":
        score -= 50.0
        reasons.append("more than one product matches the expected intended coordinates")
    elif intended_status != "unique":
        score -= 40.0
        reasons.append("intended amplicon status is %s" % intended_status)

    if search_completeness != "complete":
        score -= 20.0
        reasons.append(
            "specificity search is %s; absence of further off-targets is not exhaustive"
            % search_completeness
        )

    if n_comigrating_offtarget > 0:
        score -= 30.0 * n_comigrating_offtarget
        reasons.append("%s co-migrating off-target product(s)" % n_comigrating_offtarget)
    if n_fr_offtarget > 0 and not (
        n_comigrating_offtarget and n_fr_offtarget <= n_comigrating_offtarget
    ):
        score -= 6.0 * n_fr_offtarget
        reasons.append("%s additional F/R off-target product(s)" % n_fr_offtarget)
    if n_ff:
        score -= 5.0 * n_ff
        reasons.append("%s forward-forward (F/F) product(s)" % n_ff)
    if n_rr:
        score -= 5.0 * n_rr
        reasons.append("%s reverse-reverse (R/R) product(s)" % n_rr)

    if offtarget_min_tp5 == 0 and (n_ff or n_rr or n_fr_offtarget):
        score -= 15.0
        reasons.append("off-target with a perfectly matched 3' end (likely to prime)")

    # Historical parameter names are retained for API compatibility; these
    # signals may represent SNPs, MNPs, indels or complex VCF records.
    if snp_in_primer_3prime:
        score -= 25.0
        reasons.append("variant within a primer's 3' end (allele-specific / may fail)")
    elif snp_in_primer:
        score -= 8.0
        reasons.append("variant under a primer (away from 3' end)")

    if tm_diff > 3.0:
        score -= min(10.0, (tm_diff - 3.0) * 3.0)
        reasons.append("Tm difference %.1f degC" % tm_diff)
    if any(gc is not None and (gc < 30.0 or gc > 70.0) for gc in (gc_f, gc_r)):
        score -= 4.0
        reasons.append("primer GC outside 30-70%")

    if not gel_distinguishable:
        reasons.append("off-target not resolvable from the intended band by size")
    if conserved_fraction is not None and conserved_fraction < 1.0:
        score -= 20.0 * (1.0 - conserved_fraction)
        reasons.append("amplicon conserved in only %.0f%% of references" % (conserved_fraction * 100.0))
    if repeat_overlap:
        score -= 10.0
        reasons.append("primer overlaps an annotated repeat region")
    if dimer_concern:
        score -= 12.0
        dg_text = " (ΔG %s kcal/mol)" % cross_dimer_dg if cross_dimer_dg is not None else ""
        reasons.append("primer-dimer / hairpin likely to form%s" % dg_text)

    primes_offtarget = (
        offtarget_min_tp5 == 0 and (n_ff or n_rr or n_fr_offtarget)
    )
    score = max(0.0, min(100.0, score))
    if intended_status != "unique":
        level = "high"
    elif search_completeness != "complete":
        # Incomplete evidence may still be useful, but it can never support a
        # green/low-risk order recommendation.
        level = "medium" if score >= 55.0 and n_comigrating_offtarget == 0 else "high"
    elif (score >= 80.0 and n_comigrating_offtarget == 0
          and not snp_in_primer_3prime and not primes_offtarget
          and not dimer_concern):
        level = "low"
    elif score >= 55.0 and n_comigrating_offtarget == 0:
        level = "medium"
    else:
        level = "high"

    return RiskAssessment(
        level=level,
        score=round(score, 1),
        reasons=reasons or ["no concerns detected"],
        signals={
            "intended_status": intended_status,
            "search_completeness": search_completeness,
            "n_comigrating_offtarget": n_comigrating_offtarget,
            "n_ff": n_ff,
            "n_rr": n_rr,
            "n_fr_offtarget": n_fr_offtarget,
            "offtarget_min_tp5": offtarget_min_tp5,
            "snp_in_primer": snp_in_primer,
            "snp_in_primer_3prime": snp_in_primer_3prime,
            "gel_distinguishable": gel_distinguishable,
            "conserved_fraction": conserved_fraction,
        },
    )
