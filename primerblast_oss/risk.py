"""Experimenter-facing risk assessment.

Rolls up every signal the pipeline produces -- off-target products (F/R, F/F,
R/R), 3'-end mismatch of off-targets, SNPs under primers, Tm/GC quality,
gel-resolvability, cross-reference conservation -- into a single
high / medium / low call with the reasons spelled out, so a wet-lab user can
decide at a glance whether a pair is worth ordering.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class RiskAssessment:
    level: str                       # "low" | "medium" | "high"
    score: float                     # 0 (worst) .. 100 (best)
    reasons: List[str] = field(default_factory=list)
    signals: Dict = field(default_factory=dict)


def assess_risk(
    n_comigrating_offtarget: int = 0,
    n_ff: int = 0,
    n_rr: int = 0,
    n_fr_offtarget: int = 0,
    offtarget_min_tp5: Optional[int] = None,  # min 3'-5bp mismatch among off-targets
    snp_in_primer: bool = False,
    snp_in_primer_3prime: bool = False,
    tm_diff: float = 0.0,
    gc_f: Optional[float] = None,
    gc_r: Optional[float] = None,
    gel_distinguishable: bool = True,
    conserved_fraction: Optional[float] = None,  # 0..1 across refs, if requested
    repeat_overlap: bool = False,
    dimer_concern: bool = False,     # concerning self/cross-dimer or hairpin
    cross_dimer_dg: Optional[float] = None,
) -> RiskAssessment:
    score = 100.0
    reasons: List[str] = []

    # off-target products
    if n_comigrating_offtarget > 0:
        score -= 30.0 * n_comigrating_offtarget
        reasons.append(f"{n_comigrating_offtarget} co-migrating off-target product(s)")
    if n_fr_offtarget > 0 and not (n_comigrating_offtarget and n_fr_offtarget <= n_comigrating_offtarget):
        score -= 6.0 * n_fr_offtarget
        reasons.append(f"{n_fr_offtarget} additional F/R off-target product(s)")
    if n_ff:
        score -= 5.0 * n_ff
        reasons.append(f"{n_ff} forward-forward (F/F) product(s)")
    if n_rr:
        score -= 5.0 * n_rr
        reasons.append(f"{n_rr} reverse-reverse (R/R) product(s)")

    # an off-target whose 3' end matches perfectly is the dangerous kind
    if offtarget_min_tp5 is not None and offtarget_min_tp5 == 0 and (n_ff or n_rr or n_fr_offtarget):
        score -= 15.0
        reasons.append("off-target with a perfectly matched 3' end (likely to prime)")

    # SNPs under primers
    if snp_in_primer_3prime:
        score -= 25.0
        reasons.append("SNP within a primer's 3' end (allele-specific / may fail)")
    elif snp_in_primer:
        score -= 8.0
        reasons.append("SNP under a primer (away from 3' end)")

    # primer quality
    if tm_diff > 3.0:
        score -= min(10.0, (tm_diff - 3.0) * 3.0)
        reasons.append(f"Tm difference {tm_diff:.1f} degC")
    for gc in (gc_f, gc_r):
        if gc is not None and (gc < 30.0 or gc > 70.0):
            score -= 4.0
            reasons.append("primer GC outside 30-70%")
            break

    if not gel_distinguishable:
        reasons.append("off-target not resolvable from the intended band by size")

    if conserved_fraction is not None and conserved_fraction < 1.0:
        score -= 20.0 * (1.0 - conserved_fraction)
        reasons.append(f"amplicon conserved in only {conserved_fraction*100:.0f}% of references")

    if repeat_overlap:
        score -= 10.0
        reasons.append("primer overlaps an annotated repeat region")

    if dimer_concern:
        score -= 12.0
        dg_txt = f" (ΔG {cross_dimer_dg} kcal/mol)" if cross_dimer_dg is not None else ""
        reasons.append(f"primer-dimer / hairpin likely to form{dg_txt}")

    # an off-target that would actually prime (perfect 3' end) blocks "low"
    primes_offtarget = (offtarget_min_tp5 is not None and offtarget_min_tp5 == 0
                        and (n_ff or n_rr or n_fr_offtarget))

    score = max(0.0, min(100.0, score))
    if (score >= 80.0 and n_comigrating_offtarget == 0
            and not snp_in_primer_3prime and not primes_offtarget
            and not dimer_concern):
        level = "low"
    elif score >= 55.0 and n_comigrating_offtarget == 0:
        level = "medium"
    else:
        level = "high"

    return RiskAssessment(
        level=level, score=round(score, 1),
        reasons=reasons or ["no concerns detected"],
        signals={
            "n_comigrating_offtarget": n_comigrating_offtarget,
            "n_ff": n_ff, "n_rr": n_rr, "n_fr_offtarget": n_fr_offtarget,
            "offtarget_min_tp5": offtarget_min_tp5,
            "snp_in_primer": snp_in_primer,
            "snp_in_primer_3prime": snp_in_primer_3prime,
            "gel_distinguishable": gel_distinguishable,
            "conserved_fraction": conserved_fraction,
        },
    )
