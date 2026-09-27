"""High-level breeding assay design pipeline."""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence

from .design import DesignParams
from .specificity import (
    SEARCH_COMPLETE,
    SpecParams,
    combine_search_completeness,
)
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
        "orientation": "%s/%s" % (amplicon.fwd_primer, amplicon.rev_primer),
        "on_target": amplicon.on_target,
        "fwd_mismatch": amplicon.fwd_mismatch,
        "rev_mismatch": amplicon.rev_mismatch,
        "fwd_tp5": amplicon.fwd_tp5,
        "rev_tp5": amplicon.rev_tp5,
        "fwd_end3": amplicon.fwd_end3,
        "rev_end3": amplicon.rev_end3,
    }


def _overlaps(amplicon, chrom: str, low: int, high: int) -> bool:
    if amplicon.subject != chrom:
        return False
    amplicon_low = min(amplicon.start, amplicon.end)
    amplicon_high = max(amplicon.start, amplicon.end)
    return amplicon_low <= high and amplicon_high >= low


def expected_amplicon_from_design(pair, template: Template) -> Dict:
    left_5p = template.to_genomic(pair.left_start)
    right_5p = template.to_genomic(pair.right_start)
    if template.anchor_strand == "+":
        start, end = left_5p, right_5p
        left_name, right_name = "F", "R"
    else:
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
        "orientation": "%s/%s" % (left_name, right_name),
        "fwd_primer": left_name,
        "rev_primer": right_name,
    }


def _positions_from_expected(expected: Dict, len_f: int,
                             len_r: int) -> Dict[str, List[int]]:
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
    expected_primer_mismatches: Optional[Mapping[str, int]] = None,
) -> Dict:
    allowed_mismatches = {
        name: max(0, int(value))
        for name, value in (expected_primer_mismatches or {}).items()
    }
    all_products = (
        list(design_res.get("on_target", []))
        + list(design_res.get("off_target", []))
    )
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
                and amplicon.fwd_mismatch
                <= allowed_mismatches.get(amplicon.fwd_primer, 0)
                and amplicon.rev_mismatch
                <= allowed_mismatches.get(amplicon.rev_primer, 0)
            )
            if exact:
                candidates.append(amplicon)
    else:
        if chrom is None or ext_start is None or ext_end is None:
            raise ValueError(
                "reclassify_by_anchor needs template+pair or legacy interval arguments")
        reference_size = designed_size
        for amplicon in all_products:
            proper = {amplicon.fwd_primer, amplicon.rev_primer} == {"F", "R"}
            if (proper and _overlaps(amplicon, chrom, ext_start, ext_end)
                    and amplicon.fwd_mismatch
                    <= allowed_mismatches.get(amplicon.fwd_primer, 0)
                    and amplicon.rev_mismatch
                    <= allowed_mismatches.get(amplicon.rev_primer, 0)):
                candidates.append(amplicon)

    if len(candidates) == 1:
        intended_status = "unique"
        candidates[0].on_target = True
        on_target = [candidates[0]]
    elif not candidates:
        intended_status = "missing"
        on_target = []
    else:
        intended_status = "ambiguous"
        on_target = []
        for candidate in candidates:
            candidate.__dict__["ambiguous_intended_candidate"] = True

    off_target = [amplicon for amplicon in all_products if not amplicon.on_target]
    comigrating = [
        amplicon for amplicon in off_target
        if reference_size is not None
        and abs(amplicon.size - reference_size) < gel_min_gap
    ]
    nearest_gap = (
        min(abs(amplicon.size - reference_size) for amplicon in off_target)
        if reference_size is not None and off_target else None
    )

    search_completeness = design_res.get("search_completeness", SEARCH_COMPLETE)
    observed_specific = intended_status == "unique" and len(off_target) == 0
    if not observed_specific:
        specific = False
        specificity_status = "non_specific"
    elif search_completeness == SEARCH_COMPLETE:
        specific = True
        specificity_status = "specific"
    else:
        specific = None
        specificity_status = "indeterminate"

    return {
        **design_res,
        "n_products": len(all_products),
        "on_target": on_target,
        "off_target": off_target,
        "n_on_target": len(on_target),
        "n_off_target": len(off_target),
        "n_comigrating": len(comigrating),
        "nearest_offtarget_gap": nearest_gap,
        "gel_distinguishable": len(comigrating) == 0,
        "specific_observed": observed_specific,
        "specific": specific,
        "specificity_status": specificity_status,
        "intended_status": intended_status,
        "expected_amplicon": expected,
    }


def _variant_record_dict(variant) -> Dict:
    end = getattr(
        variant, "end", int(variant.pos) + max(1, len(str(variant.ref))) - 1)
    return {
        "pos": int(variant.pos),
        "end": int(end),
        "ref": variant.ref,
        "alt": list(variant.alt),
        "kind": variant_kind(variant),
    }


def _attach_marker_gel_context(
    caps_info: Optional[Dict],
    per_db_products: Sequence[Dict],
    *,
    intended_status: str,
    search_completeness: str,
    specific_all_db,
    genomes_by_db: Optional[Mapping[str, object]] = None,
    primers: Optional[Mapping[str, str]] = None,
) -> Optional[Dict]:
    if not caps_info:
        return caps_info
    gel_analysis = caps_info.get("gel_analysis")
    if not gel_analysis:
        best = caps_info.get("best_result") or {}
        if caps_info.get("best_marker_type") == "dCAPS":
            gel_analysis = ((best.get("digest") or {}).get("gel_analysis"))
        else:
            gel_analysis = best.get("gel_analysis")
    if not gel_analysis:
        return caps_info

    enzyme = caps_info.get("best_enzyme")
    if enzyme and primers is not None:
        from .amplicon_digest import digest_offtarget_products
        annotated = digest_offtarget_products(
            per_db_products, genomes_by_db, primers, enzyme)
        if isinstance(per_db_products, list):
            per_db_products[:] = annotated
        else:
            per_db_products = annotated

    from .gel import assess_background_across_databases
    background = assess_background_across_databases(
        gel_analysis, per_db_products)
    intrinsic = bool(
        (gel_analysis.get("genotype_discrimination") or {}).get(
            "distinguishable"))
    background_complete = background.get("all_databases_complete") is True
    background_ok = background.get("all_databases_distinguishable") is True

    if intended_status != "unique":
        verdict = "invalid_intended"
    elif search_completeness != SEARCH_COMPLETE:
        verdict = "indeterminate_search"
    elif not intrinsic:
        verdict = "ambiguous_genotype_pattern"
    elif not background_complete:
        verdict = "indeterminate_offtarget_digest"
    elif not background_ok:
        verdict = "ambiguous_with_offtarget_background"
    elif specific_all_db is True:
        verdict = "specific_clean"
    else:
        verdict = "gel_scorable_with_separated_offtargets"

    caps_info["gel_analysis"] = gel_analysis
    caps_info["background_analysis"] = background
    caps_info["gel_scorable_all_db"] = (
        intended_status == "unique"
        and search_completeness == SEARCH_COMPLETE
        and intrinsic
        and background_complete
        and background_ok
    )
    caps_info["marker_verdict"] = verdict
    return caps_info


def analyze_pair(pair, per_db: Sequence[Dict], design_db: str,
                 template: Optional[Template], variants: Sequence,
                 caps_info: Optional[Dict], gel_min_gap: int = 50,
                 dimer_params=None,
                 expected_primer_mismatches: Optional[Mapping[str, int]] = None,
                 genomes_by_db: Optional[Mapping[str, object]] = None) -> Dict:
    len_f, len_r = len(pair.forward), len(pair.reverse)
    design_res = next(
        (result for result in per_db if result["db"] == design_db), per_db[0])
    if template is not None:
        design_res = reclassify_by_anchor(
            design_res,
            template=template,
            pair=pair,
            gel_min_gap=gel_min_gap,
            expected_primer_mismatches=expected_primer_mismatches,
        )

    per_db_views = [
        design_res if result["db"] == design_res["db"] else result
        for result in per_db
    ]
    overall_completeness = combine_search_completeness([
        view.get("search_completeness", SEARCH_COMPLETE) for view in per_db_views
    ])
    db_specific_values = [view.get("specific") for view in per_db_views]
    if any(value is False for value in db_specific_values):
        specific_all_db = False
        specificity_status_all_db = "non_specific"
    elif any(value is None for value in db_specific_values):
        specific_all_db = None
        specificity_status_all_db = "indeterminate"
    elif db_specific_values and all(value is True for value in db_specific_values):
        specific_all_db = True
        specificity_status_all_db = "specific"
    else:
        specific_all_db = False
        specificity_status_all_db = "non_specific"

    per_db_products = []
    for view in per_db_views:
        per_db_products.append({
            "db": view["db"],
            "n_products": view.get("n_products", 0),
            "n_on_target": view.get("n_on_target", 0),
            "n_off_target": view.get("n_off_target", 0),
            "n_comigrating": view.get("n_comigrating", 0),
            "specific": view.get("specific"),
            "specific_observed": view.get("specific_observed"),
            "specificity_status": view.get("specificity_status"),
            "search_completeness": view.get(
                "search_completeness", SEARCH_COMPLETE),
            "primer_search_completeness": view.get(
                "primer_search_completeness", {}),
            "completeness_recommendation": view.get(
                "completeness_recommendation"),
            "thermo_status": view.get("thermo_status"),
            "thermo_evaluated": view.get("thermo_evaluated", False),
            "thermo_genome_fasta": view.get("thermo_genome_fasta"),
            "thermo_genome_association": view.get(
                "thermo_genome_association"),
            "intended_status": view.get("intended_status"),
            "expected_amplicon": view.get("expected_amplicon"),
            "gel_distinguishable": view.get("gel_distinguishable", True),
            "nearest_offtarget_gap": view.get("nearest_offtarget_gap"),
            "products": [_amp_dict(amplicon) for amplicon in (
                list(view.get("on_target", []))
                + list(view.get("off_target", []))
            )],
        })

    on_target = list(design_res.get("on_target", []))
    off_target = list(design_res.get("off_target", []))
    intended_status = design_res.get(
        "intended_status", "unique" if len(on_target) == 1 else "missing")
    n_ff = sum(_orientation_kind(amplicon) == "FF" for amplicon in off_target)
    n_rr = sum(_orientation_kind(amplicon) == "RR" for amplicon in off_target)
    n_fr = sum(_orientation_kind(amplicon) == "FR" for amplicon in off_target)
    offtarget_min_tp5 = min(
        (min(amplicon.fwd_tp5, amplicon.rev_tp5) for amplicon in off_target),
        default=None,
    )

    intended = on_target[0] if len(on_target) == 1 else None
    site_variants: List = []
    amplicon_span_variants: List = []
    footprints = (
        footprints_from_amplicon(intended, len_f, len_r)
        if intended is not None else []
    )
    if intended is not None and variants:
        site_variants = snps_under_primers(footprints, variants)
        amplicon_span_variants = amplicon_variants(
            intended.subject, intended.start, intended.end, variants)
    variant_in_primer = bool(site_variants)
    variant_in_primer_3prime = any(
        site.in_3prime_5bp for site in site_variants)

    conservation = conservation_from_per_db(per_db_views, pair.product_size)

    from . import dimers as dimer_module
    dimer = (
        dimer_module.analyze_pair(pair.forward, pair.reverse, dimer_params)
        if dimer_module.available() else None
    )

    caps_info = _attach_marker_gel_context(
        caps_info,
        per_db_products,
        intended_status=intended_status,
        search_completeness=overall_completeness,
        specific_all_db=specific_all_db,
        genomes_by_db=genomes_by_db,
        primers={"F": pair.forward, "R": pair.reverse},
    )

    risk = assess_risk(
        intended_status=intended_status,
        search_completeness=overall_completeness,
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
        conserved_fraction=(
            conservation["n_conserved"] / conservation["n_refs"]
            if conservation["n_refs"] else None
        ),
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
    amplicon_variant_dicts = [
        _variant_record_dict(variant) for variant in amplicon_span_variants
    ]
    thermo_status_by_db = {
        view["db"]: view.get("thermo_status") for view in per_db_views
    }
    thermo_genomes_by_db = {
        view["db"]: view.get("thermo_genome_fasta") for view in per_db_views
    }

    return {
        "name": "%s_P%s" % (pair.template_id, pair.index + 1),
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
        "specific": design_res.get("specific"),
        "specific_observed": design_res.get("specific_observed"),
        "specificity_status": design_res.get("specificity_status"),
        "specific_all_db": specific_all_db,
        "specificity_status_all_db": specificity_status_all_db,
        "search_completeness": overall_completeness,
        "search_complete_all_db": overall_completeness == SEARCH_COMPLETE,
        "incomplete_databases": [
            view["db"] for view in per_db_views
            if view.get("search_completeness", SEARCH_COMPLETE) != SEARCH_COMPLETE
        ],
        "thermo_status_by_db": thermo_status_by_db,
        "thermo_genomes_by_db": thermo_genomes_by_db,
        "risk": risk.level,
        "risk_score": risk.score,
        "risk_reasons": risk.reasons,
        "n_off_target": len(off_target),
        "n_ff": n_ff,
        "n_rr": n_rr,
        "n_fr_offtarget": n_fr,
        "tp5_mismatch_min": offtarget_min_tp5,
        "variant_in_primer": variant_in_primer,
        "variant_in_primer_3prime": variant_in_primer_3prime,
        "snp_in_primer": variant_in_primer,
        "snp_in_primer_3prime": variant_in_primer_3prime,
        "primer_variants": primer_variant_dicts,
        "amplicon_variants": amplicon_variant_dicts,
        "amplicon_snps": amplicon_variant_dicts,
        "conserved_refs": conservation["conserved_in"],
        "conservation": conservation,
        "caps_enzyme": (caps_info or {}).get("best_enzyme"),
        "marker_type": (caps_info or {}).get("best_marker_type"),
        "marker_verdict": (caps_info or {}).get("marker_verdict"),
        "gel_scorable_all_db": (caps_info or {}).get("gel_scorable_all_db"),
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
        "products": [
            _amp_dict(amplicon) for amplicon in on_target + off_target
        ],
    }


def build_caps(template: Template, pair, snp_local_index: int,
               alt_base: str, gel_min_gap: int = 25,
               gel_ladder: str = "auto",
               custom_ladder_bands: Optional[Sequence[int]] = None,
               gel_percent: Optional[float] = None) -> Optional[Dict]:
    """Build exact natural-CAPS digest results for the designed amplicon."""
    from .caps import caps_scan, enzymes_gained_lost, result_to_dict

    sequence = template.seq
    low = pair.left_start
    high = pair.right_start
    if not (low <= snp_local_index <= high):
        return None
    amplicon_ref = sequence[low:high + 1].upper()
    relative_index = snp_local_index - low
    ref_base = amplicon_ref[relative_index]
    amplicon_alt = (
        amplicon_ref[:relative_index]
        + alt_base.upper()
        + amplicon_ref[relative_index + 1:]
    )
    results = caps_scan(
        amplicon_ref,
        amplicon_alt,
        gel_min_gap=gel_min_gap,
        ladder=gel_ladder,
        custom_ladder_bands=custom_ladder_bands,
        gel_percent=gel_percent,
    )
    gained_lost = enzymes_gained_lost(amplicon_ref, amplicon_alt)
    best = next(
        (result for result in results if result.distinguishable), None)
    return {
        "best_marker_type": "CAPS" if best else None,
        "best_enzyme": best.enzyme if best else None,
        "best_distinguishable": bool(best),
        "allele_ref_fragments": best.allele_a_fragments if best else None,
        "allele_alt_fragments": best.allele_b_fragments if best else None,
        "min_gel_gap": best.min_gel_gap if best else None,
        "gel_analysis": best.gel_analysis if best else None,
        "best_result": result_to_dict(best) if best else None,
        "natural_candidates": [result_to_dict(result) for result in results],
        "gained": gained_lost.get("gained", []),
        "lost": gained_lost.get("lost", []),
        "n_candidate_enzymes": sum(
            result.distinguishable for result in results),
        "snp_local_index": snp_local_index,
        "snp_amplicon_index": relative_index,
        "ref_base": ref_base,
        "alt_base": alt_base.upper(),
        "allele_a_sequence": amplicon_ref,
        "allele_b_sequence": amplicon_alt,
        "dcaps": None,
    }


def _attach_dcaps(
    caps_info: Dict,
    template: Template,
    pair,
    snp_local: int,
    alt_base: str,
    databases: Sequence[str],
    spec_params: SpecParams,
    blastn_bin: Optional[str],
    genomes_by_db: Mapping[str, object],
    thermo_params,
    thermo_gate: bool,
    dimer_params,
    variants: Sequence,
    max_candidates: int,
    gel_ladder: str = "auto",
    custom_ladder_bands: Optional[Sequence[int]] = None,
    gel_percent: Optional[float] = None,
) -> Dict:
    from .dcaps_workflow import evaluate_dcaps_candidates

    dcaps = evaluate_dcaps_candidates(
        template,
        pair,
        snp_local,
        alt_base,
        databases,
        spec_params=spec_params,
        blastn_bin=blastn_bin,
        genomes_by_db=genomes_by_db,
        thermo_params=thermo_params,
        thermo_gate=thermo_gate,
        dimer_params=dimer_params,
        variants=variants,
        max_candidates_to_screen=max_candidates,
        gel_ladder=gel_ladder,
        custom_ladder_bands=custom_ladder_bands,
        gel_percent=gel_percent,
    )
    caps_info["dcaps"] = dcaps
    best = dcaps.get("best")
    if best and best.get("orderable"):
        caps_info["best_marker_type"] = "dCAPS"
        caps_info["best_enzyme"] = best["enzyme"]
        caps_info["best_distinguishable"] = True
        caps_info["allele_ref_fragments"] = best["digest"][
            "allele_a_fragments"]
        caps_info["allele_alt_fragments"] = best["digest"][
            "allele_b_fragments"]
        caps_info["min_gel_gap"] = best["digest"]["min_gel_gap"]
        caps_info["gel_analysis"] = best["digest"].get("gel_analysis")
        caps_info["best_result"] = best
    return caps_info


def _template_allele_base(template: Template, genomic_base: str) -> str:
    """Convert a genomic allele base into the oriented template base."""
    base = genomic_base.upper()
    if template.anchor_strand == "-":
        from .genome import revcomp
        return revcomp(base)
    return base


def _preferred_genotyping_mode(summary: Dict) -> Optional[str]:
    caps = summary.get("caps") or {}
    marker_verdict = summary.get("marker_verdict")
    if caps.get("best_marker_type") and marker_verdict in (
            "specific_clean", "gel_scorable_with_separated_offtargets"):
        return caps.get("best_marker_type")
    aspcr = summary.get("aspcr") or {}
    screens = (aspcr.get("specificity_screen") or {}).get("sets") or {}
    def usable_screen(name):
        screen = screens.get(name) or {}
        structures = screen.get("primer_structures") or {}
        return (screen.get("genome_screen_acceptable") is True
                and structures.get("status") == "evaluated"
                and structures.get("n_concerning") == 0)

    tetra = aspcr.get("best_tetra")
    if (tetra and tetra.get("gel_scorable")
            and usable_screen("tetra_arms")):
        return "tetra-ARMS"
    if (aspcr.get("best_classical_ref") and aspcr.get("best_classical_alt")
            and all(usable_screen("classical_" + allele) for allele in ("ref", "alt"))):
        return "AS-PCR"
    return None


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
    genomes_by_db: Optional[Mapping[str, object]] = None,
    thermo_params=None,
    thermo_gate: bool = True,
    dimer_params=None,
    gel_ladder: str = "auto",
    custom_ladder_bands: Optional[Sequence[int]] = None,
    gel_percent: Optional[float] = None,
    dcaps_pairs_to_screen: int = 3,
    dcaps_candidates_per_pair: int = 6,
    aspcr_candidates_per_allele: int = 4,
    aspcr_tetra_candidates: int = 6,
    aspcr_pairs_to_screen: int = 2,
) -> Dict:
    template = extract_template(genome, region, flank=flank)
    design_db = databases[0]
    associated_genomes = dict(genomes_by_db or {})
    associated_genomes.setdefault(design_db, genome)

    design_params = design_params or DesignParams()
    specificity = spec_params or SpecParams()
    snp_local = None
    local_alt_base = None
    if caps_snp is not None:
        snp_local = _genomic_to_local(template, caps_snp["genomic_pos"])
        local_alt_base = _template_allele_base(template, caps_snp["alt"])
        if snp_local is not None:
            from dataclasses import replace
            design_params = replace(design_params, target=(snp_local, 1))

    result = run_pipeline(
        template.id,
        template.seq,
        databases,
        design_params=design_params,
        spec_params=specificity,
        primer3_bin=primer3_bin,
        blastn_bin=blastn_bin,
        genomes_by_db=associated_genomes,
        thermo_params=thermo_params,
        thermo_gate=thermo_gate,
        dimer_params=dimer_params,
    )

    variants = variants or []
    gel_min_gap = specificity.gel_min_gap_bp
    pair_dicts: List[Dict] = []
    for pair_index, pair in enumerate(result.pairs):
        caps_info = None
        aspcr_info = None
        if caps_snp is not None and snp_local is not None and local_alt_base is not None:
            from .aspcr import build_aspcr
            aspcr_info = build_aspcr(
                template.seq,
                pair,
                snp_local,
                local_alt_base,
                max_classical_per_allele=aspcr_candidates_per_allele,
                max_tetra=aspcr_tetra_candidates,
                min_length=design_params.min_size,
                max_length=design_params.max_size,
                opt_length=design_params.opt_size,
            )
            if pair_index < aspcr_pairs_to_screen:
                from .aspcr import screen_aspcr
                aspcr_info = screen_aspcr(
                    aspcr_info,
                    databases,
                    spec_params=specificity,
                    blastn_bin=blastn_bin,
                    genomes_by_db=associated_genomes,
                    thermo_params=thermo_params,
                    thermo_gate=thermo_gate,
                    size_tolerance=10,
                    template=template,
                    design_db=design_db,
                    dimer_params=dimer_params,
                )
            else:
                aspcr_info["specificity_screen"] = {
                    "status": "skipped_pair_limit",
                    "reason": (
                        "AS-PCR/tetra-ARMS BLAST screening is limited to the "
                        "top %s parent primer pairs" % aspcr_pairs_to_screen),
                    "sets": {},
                    "search_complete_all_sets": False,
                }
            caps_info = build_caps(
                template,
                pair,
                snp_local,
                local_alt_base,
                gel_min_gap=gel_min_gap,
                gel_ladder=gel_ladder,
                custom_ladder_bands=custom_ladder_bands,
                gel_percent=gel_percent,
            )
            if (caps_info is not None
                    and not caps_info.get("best_distinguishable")):
                if pair_index < dcaps_pairs_to_screen:
                    caps_info = _attach_dcaps(
                        caps_info,
                        template,
                        pair,
                        snp_local,
                        local_alt_base,
                        databases,
                        specificity,
                        blastn_bin,
                        associated_genomes,
                        thermo_params,
                        thermo_gate,
                        dimer_params,
                        variants,
                        dcaps_candidates_per_pair,
                        gel_ladder=gel_ladder,
                        custom_ladder_bands=custom_ladder_bands,
                        gel_percent=gel_percent,
                    )
                else:
                    caps_info["dcaps"] = {
                        "status": "skipped_pair_limit",
                        "reason": (
                            "dCAPS specificity screening is limited to the "
                            "top %s parent primer pairs" % dcaps_pairs_to_screen),
                        "candidates": [],
                        "n_orderable": 0,
                    }
        pair_summary = analyze_pair(
            pair,
            pair.specificity["per_db"],
            design_db,
            template,
            variants,
            caps_info,
            gel_min_gap=gel_min_gap,
            dimer_params=dimer_params,
            genomes_by_db=associated_genomes,
        )
        pair_summary["aspcr"] = aspcr_info
        pair_summary["preferred_genotyping_mode"] = _preferred_genotyping_mode(
            pair_summary)
        pair_summary["available_genotyping_modes"] = [
            mode for mode, available in (
                ("CAPS/dCAPS", bool(
                    (pair_summary.get("caps") or {}).get("best_marker_type"))),
                ("AS-PCR", bool(
                    (aspcr_info or {}).get("best_classical_ref")
                    and (aspcr_info or {}).get("best_classical_alt"))),
                ("tetra-ARMS", bool(
                    (aspcr_info or {}).get("best_tetra"))),
            ) if available
        ]
        pair_dicts.append(pair_summary)

    risk_order = {"low": 0, "medium": 1, "high": 2}
    marker_order = {
        "specific_clean": 0,
        "gel_scorable_with_separated_offtargets": 1,
        "indeterminate_search": 3,
        "indeterminate_offtarget_digest": 3,
        "ambiguous_genotype_pattern": 4,
        "ambiguous_with_offtarget_background": 4,
        "invalid_intended": 5,
        None: 6,
    }
    if caps_snp is not None:
        mode_order = {"CAPS": 0, "dCAPS": 0, "tetra-ARMS": 1, "AS-PCR": 2, None: 6}
        pair_dicts.sort(key=lambda pair_dict: (
            mode_order.get(pair_dict.get("preferred_genotyping_mode"), 5),
            marker_order.get(pair_dict.get("marker_verdict"), 2)
            if pair_dict.get("preferred_genotyping_mode") in ("CAPS", "dCAPS")
            else 0,
            -float(((pair_dict.get("caps") or {}).get("gel_analysis") or {}).get(
                "score", 0)),
            -float((((pair_dict.get("aspcr") or {}).get("best_tetra") or {}).get(
                "score", 0))),
            risk_order.get(pair_dict["risk"], 3),
            -pair_dict.get("risk_score", 0),
        ))
    else:
        pair_dicts.sort(key=lambda pair_dict: (
            risk_order.get(pair_dict["risk"], 3),
            -pair_dict.get("risk_score", 0),
        ))
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
        "template": {
            "sequence": template.seq, "chrom": region.chrom,
            "start": template.ext_start, "end": template.ext_end,
            "anchor": template.anchor_coord, "strand": template.anchor_strand,
        },
        "databases": list(databases),
        "thermo_genomes": {
            database: getattr(associated_genomes.get(database), "fasta", None)
            for database in databases
        },
        "caps_metadata": {
            "cut_coordinates": "top and bottom strand boundaries",
            "default_recommendation_policy": (
                "verified cleavage metadata and ordinary-PCR-compatible substrate"),
            "dcaps_pairs_screened": min(
                dcaps_pairs_to_screen, len(result.pairs)) if caps_snp else 0,
            "dcaps_candidates_per_pair": (
                dcaps_candidates_per_pair if caps_snp else 0),
            "aspcr_candidates_per_allele": (
                aspcr_candidates_per_allele if caps_snp else 0),
            "aspcr_tetra_candidates": (
                aspcr_tetra_candidates if caps_snp else 0),
            "aspcr_pairs_screened": min(
                aspcr_pairs_to_screen, len(result.pairs)) if caps_snp else 0,
            "gel_ladder": gel_ladder,
            "custom_ladder_bands": (
                list(custom_ladder_bands) if custom_ladder_bands else None),
            "gel_percent": gel_percent,
        },
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
    output: List[Dict] = []
    for region in regions:
        try:
            output.append(run_assay(region, genome, databases, **kwargs))
        except Exception as error:  # noqa: BLE001
            output.append({
                "target": {
                    "name": region.name,
                    "chrom": region.chrom,
                    "start": region.start,
                    "end": region.end,
                },
                "error": str(error),
                "n_pairs": 0,
                "pairs": [],
            })
    return output


def design_qtl_markers(interval: GenomicRegion, genome: Genome,
                       databases: Sequence[str], n_markers: int = 0,
                       spacing: int = 0, marker_flank: int = 300,
                       best_only: bool = True, **kwargs) -> List[Dict]:
    from .regions import tile_interval
    points = tile_interval(interval, n_markers=n_markers, spacing=spacing)
    results: List[Dict] = []
    for point in points:
        region = GenomicRegion(
            point.chrom,
            point.start - marker_flank,
            point.start + marker_flank,
            "+",
            point.name,
            "qtl-marker",
        )
        try:
            result = run_assay(region, genome, databases, flank=0, **kwargs)
        except Exception as error:  # noqa: BLE001
            results.append({
                "marker": point.name,
                "anchor": point.start,
                "error": str(error),
                "pairs": [],
            })
            continue
        if best_only and result["pairs"]:
            result = {**result, "pairs": result["pairs"][:1]}
        result["marker"] = point.name
        result["anchor"] = point.start
        results.append(result)
    return results
