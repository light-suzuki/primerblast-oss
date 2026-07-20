"""High-level breeding assay design pipeline."""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

from .design import DesignParams
from .specificity import SpecParams
from .pipeline import run_pipeline
from .genome import Genome
from .regions import GenomicRegion, Template, extract_template
from .variants import (
    amplicon_variants,
    conservation_from_per_db,
    footprints_from_amplicon,
    snps_under_primers,
    variant_kind,
)
from .risk import assess_risk


def _orientation_kind(amplicon) -> str:
    if {amplicon.fwd_primer, amplicon.rev_primer} == {"F", "R"}:
        return "FR"
    if amplicon.fwd_primer == amplicon.rev_primer == "F":
        return "FF"
    if amplicon.fwd_primer == amplicon.rev_primer == "R":
        return "RR"
    return "other"


def _amp_dict(amplicon) -> Dict:
    return {
        "subject": amplicon.subject,
        "start": amplicon.start,
        "end": amplicon.end,
        "size": amplicon.size,
        "orientation": f"{amplicon.fwd_primer}/{amplicon.rev_primer}",
        "on_target": amplicon.on_target,
        "fwd_mismatch": amplicon.fwd_mismatch,
        "rev_mismatch": amplicon.rev_mismatch,
        "fwd_tp5": amplicon.fwd_tp5,
        "rev_tp5": amplicon.rev_tp5,
    }


def _overlaps(amplicon, chrom: str, lo: int, hi: int) -> bool:
    """Legacy interval-overlap helper retained for library compatibility."""
    if amplicon.subject != chrom:
        return False
    amp_lo, amp_hi = min(amplicon.start, amplicon.end), max(amplicon.start, amplicon.end)
    return amp_lo <= hi and amp_hi >= lo


def expected_amplicon_from_design(pair, template: Template) -> Dict:
    """Map Primer3 left/right 5' coordinates to the design genome."""
    left_5p = template.to_genomic(pair.left_start)
    right_5p = template.to_genomic(pair.right_start)
    if template.anchor_strand == "+":
        start, end = left_5p, right_5p
        left_name, right_name = "F", "R"
    else:
        # A template extracted from the minus strand is reverse-complemented.
        # Primer3's right primer therefore lies on genomic plus (left side), and
        # Primer3's left primer lies on genomic minus (right side).
        start, end = right_5p, left_5p
        left_name, right_name = "R", "F"
    if start > end:
        start, end = end, start
        left_name, right_name = right_name, left_name
    return {
        "subject": template.region.chrom,
        "start": start,
        "end": end,
        "size": end - start + 1,
        "orientation": f"{left_name}/{right_name}",
        "fwd_primer": left_name,
        "rev_primer": right_name,
    }


def _positions_from_expected(expected: Dict, len_f: int, len_r: int) -> Dict[str, List[int]]:
    lengths = {"F": len_f, "R": len_r}
    left_name = expected["fwd_primer"]
    right_name = expected["rev_primer"]
    left = [expected["start"], expected["start"] + lengths[left_name] - 1]
    right = [expected["end"] - lengths[right_name] + 1, expected["end"]]
    by_name = {left_name: left, right_name: right}
    return {
        "left": left,
        "right": right,
        "forward": by_name["F"],
        "reverse": by_name["R"],
    }


def reclassify_by_anchor(
    design_res: Dict,
    chrom: Optional[str] = None,
    ext_start: Optional[int] = None,
    ext_end: Optional[int] = None,
    designed_size: Optional[int] = None,
    gel_min_gap: int = 50,
    *,
    template: Optional[Template] = None,
    pair=None,
    coordinate_tolerance: int = 0,
    size_tolerance: int = 0,
) -> Dict:
    """Reclassify the design-reference products using a genomic anchor.

    The high-level assay path passes ``template`` and ``pair`` and therefore
    requires one exact coordinate/orientation match. The older interval arguments
    remain as a conservative compatibility fallback for direct library callers.
    """
    all_products = list(design_res.get("on_target", [])) + list(design_res.get("off_target", []))
    for amplicon in all_products:
        amplicon.on_target = False

    expected = None
    candidates = []
    if template is not None and pair is not None:
        expected = expected_amplicon_from_design(pair, template)
        reference_size = expected["size"]
        for amplicon in all_products:
            exact = (
                amplicon.subject == expected["subject"]
                and abs(amplicon.start - expected["start"]) <= coordinate_tolerance
                and abs(amplicon.end - expected["end"]) <= coordinate_tolerance
                and amplicon.fwd_primer == expected["fwd_primer"]
                and amplicon.rev_primer == expected["rev_primer"]
                and abs(amplicon.size - expected["size"]) <= size_tolerance
                and amplicon.fwd_mismatch == 0
                and amplicon.rev_mismatch == 0
            )
            if exact:
                candidates.append(amplicon)
    else:
        # Compatibility fallback: no exact Primer3/template mapping is available.
        if chrom is None or ext_start is None or ext_end is None:
            raise ValueError("reclassify_by_anchor needs template+pair or legacy interval arguments")
        reference_size = designed_size
        for amplicon in all_products:
            proper = {amplicon.fwd_primer, amplicon.rev_primer} == {"F", "R"}
            if (proper and _overlaps(amplicon, chrom, ext_start, ext_end)
                    and amplicon.fwd_mismatch == 0 and amplicon.rev_mismatch == 0):
                candidates.append(amplicon)

    if len(candidates) == 1:
        intended_status = "unique"
        candidates[0].on_target = True
        on = [candidates[0]]
    elif not candidates:
        intended_status = "missing"
        on = []
    else:
        intended_status = "ambiguous"
        on = []
        for candidate in candidates:
            candidate.__dict__["ambiguous_intended_candidate"] = True

    off = [amplicon for amplicon in all_products if not amplicon.on_target]
    comigrating = ([amplicon for amplicon in off
                    if reference_size is not None
                    and abs(amplicon.size - reference_size) < gel_min_gap])
    nearest_gap = (min(abs(amplicon.size - reference_size) for amplicon in off)
                   if reference_size is not None and off else None)
    return {
        **design_res,
        "n_products": len(all_products),
        "on_target": on,
        "off_target": off,
        "n_on_target": len(on),
        "n_off_target": len(off),
        "n_comigrating": len(comigrating),
        "nearest_offtarget_gap": nearest_gap,
        "gel_distinguishable": len(comigrating) == 0,
        "specific": intended_status == "unique" and len(off) == 0,
        "intended_status": intended_status,
        "expected_amplicon": expected,
    }


def _variant_record_dict(variant) -> Dict:
    end = getattr(variant, "end", int(variant.pos) + max(1, len(str(variant.ref))) - 1)
    return {
        "pos": int(variant.pos),
        "end": int(end),
        "ref": variant.ref,
        "alt": list(variant.alt),
        "kind": variant_kind(variant),
    }


def analyze_pair(pair, per_db: Sequence[Dict], design_db: str,
                 template: Optional[Template], variants: Sequence,
                 caps_info: Optional[Dict], gel_min_gap: int = 50,
                 dimer_params=None) -> Dict:
    """Build a full experimenter-facing summary for one designed pair."""
    len_f, len_r = len(pair.forward), len(pair.reverse)
    design_res = next((result for result in per_db if result["db"] == design_db), per_db[0])
    if template is not None:
        design_res = reclassify_by_anchor(
            design_res,
            template=template,
            pair=pair,
            gel_min_gap=gel_min_gap,
        )

    per_db_views = [design_res if result["db"] == design_res["db"] else result
                    for result in per_db]
    per_db_products = []
    for view in per_db_views:
        per_db_products.append({
            "db": view["db"],
            "n_products": view.get("n_products", 0),
            "n_on_target": view.get("n_on_target", 0),
            "n_off_target": view.get("n_off_target", 0),
            "n_comigrating": view.get("n_comigrating", 0),
            "specific": view.get("specific", False),
            "intended_status": view.get("intended_status"),
            "expected_amplicon": view.get("expected_amplicon"),
            "gel_distinguishable": view.get("gel_distinguishable", True),
            "nearest_offtarget_gap": view.get("nearest_offtarget_gap"),
            "products": [_amp_dict(amplicon) for amplicon in (
                list(view.get("on_target", [])) + list(view.get("off_target", []))
            )],
        })

    on = list(design_res.get("on_target", []))
    off = list(design_res.get("off_target", []))
    intended_status = design_res.get("intended_status", "unique" if len(on) == 1 else "missing")
    n_ff = sum(_orientation_kind(amplicon) == "FF" for amplicon in off)
    n_rr = sum(_orientation_kind(amplicon) == "RR" for amplicon in off)
    n_fr = sum(_orientation_kind(amplicon) == "FR" for amplicon in off)
    offtarget_min_tp5 = min((min(a.fwd_tp5, a.rev_tp5) for a in off), default=None)

    intended = on[0] if len(on) == 1 else None
    site_variants: List = []
    amp_span_variants: List = []
    footprints = footprints_from_amplicon(intended, len_f, len_r) if intended is not None else []
    if intended is not None and variants:
        site_variants = snps_under_primers(footprints, variants)
        amp_span_variants = amplicon_variants(intended.subject, intended.start, intended.end, variants)
    variant_in_primer = bool(site_variants)
    variant_in_primer_3prime = any(site.in_3prime_5bp for site in site_variants)

    conservation = conservation_from_per_db(per_db_views, pair.product_size)

    from . import dimers as _dimers
    dimer = (_dimers.analyze_pair(pair.forward, pair.reverse, dimer_params)
             if _dimers.available() else None)

    risk = assess_risk(
        intended_status=intended_status,
        n_comigrating_offtarget=design_res.get("n_comigrating", 0),
        n_ff=n_ff,
        n_rr=n_rr,
        n_fr_offtarget=n_fr,
        offtarget_min_tp5=offtarget_min_tp5,
        snp_in_primer=variant_in_primer,
        snp_in_primer_3prime=variant_in_primer_3prime,
        tm_diff=abs(pair.tm_f - pair.tm_r),
        gc_f=pair.gc_f,
        gc_r=pair.gc_r,
        gel_distinguishable=design_res.get("gel_distinguishable", True),
        conserved_fraction=(conservation["n_conserved"] / conservation["n_refs"]
                            if conservation["n_refs"] else None),
        dimer_concern=(dimer["n_concerning"] > 0 if dimer else False),
        cross_dimer_dg=(dimer["cross_dimer_dg"] if dimer else None),
    )

    expected = design_res.get("expected_amplicon")
    if expected is None and template is not None:
        expected = expected_amplicon_from_design(pair, template)
    if expected is not None:
        positions = _positions_from_expected(expected, len_f, len_r)
    else:
        positions = {
            "left": [pair.left_start + 1, pair.left_3p + 1],
            "right": [pair.right_3p + 1, pair.right_start + 1],
            "forward": [pair.left_start + 1, pair.left_3p + 1],
            "reverse": [pair.right_3p + 1, pair.right_start + 1],
        }

    primer_variant_dicts = [{
        "primer": site.primer,
        "chrom": site.chrom,
        "pos": site.pos,
        "end": site.end,
        "ref": site.ref,
        "alt": list(site.alt),
        "kind": site.kind,
        "in_3prime_5bp": site.in_3prime_5bp,
        "distance_from_3prime": site.distance_from_3prime,
    } for site in site_variants]
    amplicon_variant_dicts = [_variant_record_dict(variant) for variant in amp_span_variants]

    return {
        "name": f"{pair.template_id}_P{pair.index + 1}",
        "forward": pair.forward,
        "reverse": pair.reverse,
        "product_size": pair.product_size,
        "tm_f": round(pair.tm_f, 1),
        "tm_r": round(pair.tm_r, 1),
        "gc_f": round(pair.gc_f, 1),
        "gc_r": round(pair.gc_r, 1),
        "left_pos": positions["left"],
        "right_pos": positions["right"],
        "forward_pos": positions["forward"],
        "reverse_pos": positions["reverse"],
        "intended_status": intended_status,
        "expected_amplicon": expected,
        "specific": design_res.get("specific", False),
        "risk": risk.level,
        "risk_score": risk.score,
        "risk_reasons": risk.reasons,
        "n_off_target": len(off),
        "n_ff": n_ff,
        "n_rr": n_rr,
        "n_fr_offtarget": n_fr,
        "tp5_mismatch_min": offtarget_min_tp5,
        "variant_in_primer": variant_in_primer,
        "variant_in_primer_3prime": variant_in_primer_3prime,
        # Backward-compatible keys; they now include all VCF record classes.
        "snp_in_primer": variant_in_primer,
        "snp_in_primer_3prime": variant_in_primer_3prime,
        "primer_variants": primer_variant_dicts,
        "amplicon_variants": amplicon_variant_dicts,
        "amplicon_snps": amplicon_variant_dicts,
        "conserved_refs": conservation["conserved_in"],
        "conservation": conservation,
        "caps_enzyme": (caps_info or {}).get("best_enzyme"),
        "caps": caps_info,
        "gel_distinguishable": design_res.get("gel_distinguishable", True),
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
        "per_db_products": per_db_products,
        "products": [_amp_dict(amplicon) for amplicon in on + off],
    }


def build_caps(template: Template, pair, snp_local_index: int,
               alt_base: str, gel_min_gap: int = 25) -> Optional[Dict]:
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
    gained_lost = enzymes_gained_lost(amp_ref, amp_alt)
    best = next((result for result in results if result.distinguishable), None)
    return {
        "best_enzyme": best.enzyme if best else None,
        "best_distinguishable": bool(best),
        "allele_ref_fragments": best.allele_a_fragments if best else None,
        "allele_alt_fragments": best.allele_b_fragments if best else None,
        "min_gel_gap": best.min_gel_gap if best else None,
        "gained": gained_lost.get("gained", []),
        "lost": gained_lost.get("lost", []),
        "n_candidate_enzymes": sum(result.distinguishable for result in results),
    }


def run_assay(
    region: GenomicRegion,
    genome: Genome,
    databases: Sequence[str],
    flank: int = 200,
    design_params: Optional[DesignParams] = None,
    spec_params: Optional[SpecParams] = None,
    variants: Optional[Sequence] = None,
    caps_snp: Optional[Dict] = None,
    primer3_bin: Optional[str] = None,
    blastn_bin: Optional[str] = None,
    thermo_params=None,
    thermo_gate: bool = True,
    dimer_params=None,
) -> Dict:
    template = extract_template(genome, region, flank=flank)
    design_db = databases[0]
    design_params = design_params or DesignParams()
    if caps_snp is not None:
        snp_local = _genomic_to_local(template, caps_snp["genomic_pos"])
        if snp_local is not None:
            from dataclasses import replace
            design_params = replace(design_params, target=(snp_local, 1))

    result = run_pipeline(
        template.id,
        template.seq,
        databases,
        design_params=design_params,
        spec_params=spec_params,
        primer3_bin=primer3_bin,
        blastn_bin=blastn_bin,
        genome=genome,
        thermo_params=thermo_params,
        thermo_gate=thermo_gate,
    )

    variants = variants or []
    gel_min_gap = spec_params.gel_min_gap_bp if spec_params else 50
    pair_dicts: List[Dict] = []
    for pair in result.pairs:
        caps_info = None
        if caps_snp is not None:
            snp_local = _genomic_to_local(template, caps_snp["genomic_pos"])
            if snp_local is not None:
                caps_info = build_caps(template, pair, snp_local, caps_snp["alt"])
        pair_dicts.append(analyze_pair(
            pair,
            pair.specificity["per_db"],
            design_db,
            template,
            variants,
            caps_info,
            gel_min_gap=gel_min_gap,
            dimer_params=dimer_params,
        ))

    order = {"low": 0, "medium": 1, "high": 2}
    pair_dicts.sort(key=lambda pair: (order.get(pair["risk"], 3), -pair.get("risk_score", 0)))
    return {
        "target": {
            "name": region.name,
            "chrom": region.chrom,
            "start": region.start,
            "end": region.end,
            "strand": region.strand,
            "source": region.source,
            "flank": flank,
        },
        "template_len": len(template.seq),
        "databases": list(databases),
        "n_pairs": len(pair_dicts),
        "pairs": pair_dicts,
    }


def _genomic_to_local(template: Template, genomic_pos: int) -> Optional[int]:
    if template.anchor_strand == "-":
        index = template.anchor_coord - genomic_pos
    else:
        index = genomic_pos - template.anchor_coord
    return index if 0 <= index < len(template.seq) else None


def run_batch(regions: Sequence[GenomicRegion], genome: Genome,
              databases: Sequence[str], **kwargs) -> List[Dict]:
    out: List[Dict] = []
    for region in regions:
        try:
            out.append(run_assay(region, genome, databases, **kwargs))
        except Exception as exc:  # noqa: BLE001
            out.append({
                "target": {"name": region.name, "chrom": region.chrom,
                           "start": region.start, "end": region.end},
                "error": str(exc),
                "n_pairs": 0,
                "pairs": [],
            })
    return out


def design_qtl_markers(interval: GenomicRegion, genome: Genome,
                       databases: Sequence[str], n_markers: int = 0,
                       spacing: int = 0, marker_flank: int = 300,
                       best_only: bool = True, **kwargs) -> List[Dict]:
    from .regions import tile_interval
    points = tile_interval(interval, n_markers=n_markers, spacing=spacing)
    results: List[Dict] = []
    for point in points:
        region = GenomicRegion(point.chrom, point.start - marker_flank,
                               point.start + marker_flank, "+", point.name, "qtl-marker")
        try:
            result = run_assay(region, genome, databases, flank=0, **kwargs)
        except Exception as exc:  # noqa: BLE001
            results.append({"marker": point.name, "anchor": point.start,
                            "error": str(exc), "pairs": []})
            continue
        if best_only and result["pairs"]:
            result = {**result, "pairs": result["pairs"][:1]}
        result["marker"] = point.name
        result["anchor"] = point.start
        results.append(result)
    return results
