"""High-level assay design: the full breeding-oriented pipeline.

target (gene / interval / SNP / coords)
  -> extract template from a local genome (+flank, strand-aware)
  -> Primer3 design
  -> specificity across one or more reference genomes (intended vs F/F, R/R, F/R
     off-targets, with genomic coordinates and 3'-end mismatch)
  -> SNPs under primers (VCF) and amplicon conservation across references
  -> optional CAPS/dCAPS enzyme scan for a SNP
  -> experimenter risk (low/medium/high)
Returns plain dicts ready for the CSV / BED / HTML / order-table writers.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

from .design import DesignParams
from .specificity import SpecParams
from .pipeline import run_pipeline
from .genome import Genome, revcomp
from .regions import GenomicRegion, Template, extract_template
from .variants import (
    footprints_from_amplicon, snps_under_primers, amplicon_variants,
    conservation_from_per_db,
)
from .risk import assess_risk


def _orientation_kind(a) -> str:
    if {a.fwd_primer, a.rev_primer} == {"F", "R"}:
        return "FR"
    if a.fwd_primer == a.rev_primer == "F":
        return "FF"
    if a.fwd_primer == a.rev_primer == "R":
        return "RR"
    return "other"


def _amp_dict(a) -> Dict:
    return {
        "subject": a.subject, "start": a.start, "end": a.end, "size": a.size,
        "orientation": f"{a.fwd_primer}/{a.rev_primer}", "on_target": a.on_target,
        "fwd_mismatch": a.fwd_mismatch, "rev_mismatch": a.rev_mismatch,
        "fwd_tp5": a.fwd_tp5, "rev_tp5": a.rev_tp5,
    }


def _overlaps(a, chrom: str, lo: int, hi: int) -> bool:
    if a.subject != chrom:
        return False
    a_lo, a_hi = min(a.start, a.end), max(a.start, a.end)
    return a_lo <= hi and a_hi >= lo


def reclassify_by_anchor(design_res: Dict, chrom: str, ext_start: int, ext_end: int,
                         designed_size: int, gel_min_gap: int = 50) -> Dict:
    """Re-decide intended vs off-target using the template's genomic locus.

    The generic per-db result marks *every* perfect-size proper product as
    on_target; but only the product at the template's own locus is truly
    intended. Perfect products elsewhere are paralog/duplication off-targets --
    exactly the case the generic heuristic hides."""
    all_products = list(design_res.get("on_target", [])) + list(design_res.get("off_target", []))
    on, off = [], []
    for a in all_products:
        proper = {a.fwd_primer, a.rev_primer} == {"F", "R"}
        intended = (proper and _overlaps(a, chrom, ext_start, ext_end)
                    and a.fwd_mismatch == 0 and a.rev_mismatch == 0)
        a.on_target = intended
        (on if intended else off).append(a)
    ref_size = designed_size
    comig = [a for a in off if abs(a.size - ref_size) < gel_min_gap]
    return {
        **design_res, "on_target": on, "off_target": off,
        "n_on_target": len(on), "n_off_target": len(off),
        "n_comigrating": len(comig),
        "gel_distinguishable": len(comig) == 0,
        "specific": len(off) == 0 and len(on) >= 1,
    }


def analyze_pair(pair, per_db: Sequence[Dict], design_db: str,
                 template: Optional[Template], variants: Sequence,
                 caps_info: Optional[Dict], gel_min_gap: int = 50) -> Dict:
    """Build a full summary dict for one designed pair."""
    len_f, len_r = len(pair.forward), len(pair.reverse)
    design_res = next((d for d in per_db if d["db"] == design_db), per_db[0])
    if template is not None:
        design_res = reclassify_by_anchor(
            design_res, template.region.chrom, template.ext_start, template.ext_end,
            pair.product_size, gel_min_gap)

    per_db_products = []
    for db_res in per_db:
        view = design_res if db_res["db"] == design_res["db"] else db_res
        per_db_products.append({
            "db": view["db"],
            "n_products": view.get("n_products", 0),
            "n_on_target": view.get("n_on_target", 0),
            "n_off_target": view.get("n_off_target", 0),
            "n_comigrating": view.get("n_comigrating", 0),
            "specific": view.get("specific", False),
            "gel_distinguishable": view.get("gel_distinguishable", True),
            "nearest_offtarget_gap": view.get("nearest_offtarget_gap"),
            "products": [_amp_dict(a) for a in (
                list(view.get("on_target", [])) + list(view.get("off_target", []))
            )],
        })

    on = list(design_res.get("on_target", []))
    off = list(design_res.get("off_target", []))
    n_ff = sum(1 for a in off if _orientation_kind(a) == "FF")
    n_rr = sum(1 for a in off if _orientation_kind(a) == "RR")
    n_fr = sum(1 for a in off if _orientation_kind(a) == "FR")
    offtarget_min_tp5 = min((min(a.fwd_tp5, a.rev_tp5) for a in off), default=None)

    # genomic footprints + SNPs under primers (needs an intended amplicon)
    site_snps: List = []
    snp_in_primer = snp_in_primer_3prime = False
    amp_span_variants: List = []
    intended = on[0] if on else None
    if intended is not None and variants:
        fps = footprints_from_amplicon(intended, len_f, len_r)
        site_snps = snps_under_primers(fps, variants)
        snp_in_primer = len(site_snps) > 0
        snp_in_primer_3prime = any(s.in_3prime_5bp for s in site_snps)
        amp_span_variants = amplicon_variants(intended.subject, intended.start,
                                              intended.end, variants)

    conservation = conservation_from_per_db(per_db, pair.product_size)

    risk = assess_risk(
        n_comigrating_offtarget=design_res.get("n_comigrating", 0),
        n_ff=n_ff, n_rr=n_rr, n_fr_offtarget=n_fr,
        offtarget_min_tp5=offtarget_min_tp5,
        snp_in_primer=snp_in_primer, snp_in_primer_3prime=snp_in_primer_3prime,
        tm_diff=abs(pair.tm_f - pair.tm_r), gc_f=pair.gc_f, gc_r=pair.gc_r,
        gel_distinguishable=design_res.get("gel_distinguishable", True),
        conserved_fraction=(conservation["n_conserved"] / conservation["n_refs"]
                            if conservation["n_refs"] else None),
    )

    left_pos = [intended.start, intended.start + len_f - 1] if intended else \
        [pair.left_start + 1, pair.left_3p + 1]
    right_pos = [intended.end - len_r + 1, intended.end] if intended else \
        [pair.right_3p + 1, pair.right_start + 1]

    return {
        "name": f"{pair.template_id}_P{pair.index+1}",
        "forward": pair.forward, "reverse": pair.reverse,
        "product_size": pair.product_size,
        "tm_f": round(pair.tm_f, 1), "tm_r": round(pair.tm_r, 1),
        "gc_f": round(pair.gc_f, 1), "gc_r": round(pair.gc_r, 1),
        "left_pos": left_pos, "right_pos": right_pos,
        "risk": risk.level, "risk_score": risk.score, "risk_reasons": risk.reasons,
        "n_off_target": len(off), "n_ff": n_ff, "n_rr": n_rr, "n_fr_offtarget": n_fr,
        "tp5_mismatch_min": offtarget_min_tp5,
        "snp_in_primer": snp_in_primer, "snp_in_primer_3prime": snp_in_primer_3prime,
        "amplicon_snps": [{"pos": v.pos, "ref": v.ref, "alt": list(v.alt)}
                          for v in amp_span_variants],
        "conserved_refs": conservation["conserved_in"],
        "conservation": conservation,
        "caps_enzyme": (caps_info or {}).get("best_enzyme"),
        "caps": caps_info,
        "gel_distinguishable": design_res.get("gel_distinguishable", True),
        "per_db_products": per_db_products,
        "products": [_amp_dict(a) for a in on + off],
    }


def build_caps(template: Template, pair, snp_local_index: int,
               alt_base: str, gel_min_gap: int = 25) -> Optional[Dict]:
    """CAPS/dCAPS scan: build the two allele amplicons (ref vs alt at the SNP)
    and find enzymes that digest them differently."""
    from .caps import caps_scan, enzymes_gained_lost
    seq = template.seq
    lo = pair.left_start
    hi = pair.right_start
    if not (lo <= snp_local_index <= hi):
        return None
    amp_ref = seq[lo:hi + 1]
    rel = snp_local_index - lo
    amp_alt = amp_ref[:rel] + alt_base.upper() + amp_ref[rel + 1:]
    results = caps_scan(amp_ref, amp_alt, gel_min_gap=gel_min_gap)
    gl = enzymes_gained_lost(amp_ref, amp_alt)
    best = None
    for r in results:
        if r.distinguishable:
            best = r
            break
    return {
        "best_enzyme": best.enzyme if best else None,
        "best_distinguishable": bool(best),
        "allele_ref_fragments": best.allele_a_fragments if best else None,
        "allele_alt_fragments": best.allele_b_fragments if best else None,
        "min_gel_gap": best.min_gel_gap if best else None,
        "gained": gl.get("gained", []), "lost": gl.get("lost", []),
        "n_candidate_enzymes": sum(1 for r in results if r.distinguishable),
    }


def run_assay(
    region: GenomicRegion,
    genome: Genome,
    databases: Sequence[str],
    flank: int = 200,
    design_params: Optional[DesignParams] = None,
    spec_params: Optional[SpecParams] = None,
    variants: Optional[Sequence] = None,
    caps_snp: Optional[Dict] = None,   # {"genomic_pos":int, "alt":str}
    primer3_bin: Optional[str] = None,
    blastn_bin: Optional[str] = None,
    thermo_params=None,
    thermo_gate: bool = True,
) -> Dict:
    """Design and fully evaluate primers for one target region."""
    template = extract_template(genome, region, flank=flank)
    design_db = databases[0]

    # for a CAPS/dCAPS assay the amplicon must span the SNP: target it so every
    # designed product flanks it.
    dp = design_params or DesignParams()
    if caps_snp is not None:
        snp_local0 = _genomic_to_local(template, caps_snp["genomic_pos"])
        if snp_local0 is not None:
            from dataclasses import replace
            dp = replace(dp, target=(snp_local0, 1))

    result = run_pipeline(
        template.id, template.seq, databases,
        design_params=dp, spec_params=spec_params,
        primer3_bin=primer3_bin, blastn_bin=blastn_bin,
        genome=genome, thermo_params=thermo_params, thermo_gate=thermo_gate,
    )

    variants = variants or []
    gel_min_gap = spec_params.gel_min_gap_bp if spec_params else 50
    pair_dicts: List[Dict] = []
    for pair in result.pairs:
        per_db = pair.specificity["per_db"]
        caps_info = None
        if caps_snp is not None:
            snp_local = _genomic_to_local(template, caps_snp["genomic_pos"])
            if snp_local is not None:
                caps_info = build_caps(template, pair, snp_local, caps_snp["alt"])
        pair_dicts.append(
            analyze_pair(pair, per_db, design_db, template, variants, caps_info,
                         gel_min_gap=gel_min_gap))

    # order by risk then specificity score
    order = {"low": 0, "medium": 1, "high": 2}
    pair_dicts.sort(key=lambda p: (order.get(p["risk"], 3), -p.get("risk_score", 0)))

    return {
        "target": {
            "name": region.name, "chrom": region.chrom,
            "start": region.start, "end": region.end, "strand": region.strand,
            "source": region.source, "flank": flank,
        },
        "template_len": len(template.seq),
        "databases": list(databases),
        "n_pairs": len(pair_dicts),
        "pairs": pair_dicts,
    }


def _genomic_to_local(template: Template, genomic_pos: int) -> Optional[int]:
    """Inverse of Template.to_genomic: genomic coord -> 0-based template index."""
    if template.anchor_strand == "-":
        idx = template.anchor_coord - genomic_pos
    else:
        idx = genomic_pos - template.anchor_coord
    if 0 <= idx < len(template.seq):
        return idx
    return None


def run_batch(regions: Sequence[GenomicRegion], genome: Genome,
              databases: Sequence[str], **kwargs) -> List[Dict]:
    """Run run_assay over many targets (genes, intervals, BED). One dict each."""
    out: List[Dict] = []
    for region in regions:
        try:
            out.append(run_assay(region, genome, databases, **kwargs))
        except Exception as e:  # noqa: BLE001 -- keep the batch going, record failure
            out.append({"target": {"name": region.name, "chrom": region.chrom,
                                    "start": region.start, "end": region.end},
                        "error": str(e), "n_pairs": 0, "pairs": []})
    return out


def design_qtl_markers(interval: GenomicRegion, genome: Genome,
                       databases: Sequence[str], n_markers: int = 0,
                       spacing: int = 0, marker_flank: int = 300,
                       best_only: bool = True, **kwargs) -> List[Dict]:
    """Place evenly spaced markers across a QTL interval and design a primer
    pair at each. With best_only, keep just the lowest-risk pair per marker."""
    from .regions import tile_interval
    points = tile_interval(interval, n_markers=n_markers, spacing=spacing)
    results: List[Dict] = []
    for pt in points:
        region = GenomicRegion(pt.chrom, pt.start - marker_flank,
                               pt.start + marker_flank, "+", pt.name, "qtl-marker")
        try:
            res = run_assay(region, genome, databases, flank=0, **kwargs)
        except Exception as e:  # noqa: BLE001
            results.append({"marker": pt.name, "anchor": pt.start, "error": str(e),
                            "pairs": []})
            continue
        if best_only and res["pairs"]:
            res = {**res, "pairs": res["pairs"][:1]}
        res["marker"] = pt.name
        res["anchor"] = pt.start
        results.append(res)
    return results
