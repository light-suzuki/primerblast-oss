"""Render results as JSON, TSV, or human-readable text."""
from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from typing import Dict, List

from .pipeline import PipelineResult


def _amp_dict(amplicon) -> dict:
    data = asdict(amplicon) if is_dataclass(amplicon) else dict(amplicon)
    if hasattr(amplicon, "__dict__"):
        for key, value in amplicon.__dict__.items():
            data.setdefault(key, value)
    data["orientation"] = "%s/%s" % (
        data["fwd_primer"], data["rev_primer"])
    return data


def _specificity_to_dict(spec: dict) -> dict:
    data = dict(spec)
    per_database = []
    for database in spec.get("per_db", []):
        database_dict = dict(database)
        for key in ("on_target", "off_target"):
            database_dict[key] = [
                _amp_dict(amplicon) for amplicon in database.get(key, [])
            ]
        per_database.append(database_dict)
    data["per_db"] = per_database
    return data


def to_dict(result: PipelineResult) -> dict:
    pairs = []
    for pair in result.pairs:
        data = asdict(pair)
        data["specificity"] = _specificity_to_dict(pair.specificity)
        pairs.append(data)
    return {
        "mode": "design",
        "template_id": result.template_id,
        "template_len": result.template_len,
        "databases": result.databases,
        "params": result.params,
        "primer3_explain": result.primer3_explain,
        "n_pairs": len(result.pairs),
        "pairs": pairs,
    }


def to_json(result: PipelineResult, indent: int = 2) -> str:
    return json.dumps(to_dict(result), indent=indent, default=str)


def _thermo_summary(per_database: List[Dict]) -> str:
    return ";".join(
        "%s=%s:%s" % (
            database.get("db", "").split("/")[-1],
            database.get("thermo_status", "unknown"),
            database.get("thermo_genome_fasta") or "-",
        )
        for database in per_database
    )


def to_tsv(result: PipelineResult) -> str:
    columns = [
        "rank", "score", "specificity_status", "specific_all_db",
        "specific_observed_all_db", "search_completeness",
        "search_complete_all_db", "incomplete_databases", "thermo_by_db",
        "gel_distinguishable", "forward", "reverse", "product_size",
        "tm_f", "tm_r", "gc_f", "gc_r", "total_on_target",
        "total_off_target", "total_comigrating", "dimer_ok",
        "cross_dimer_dg",
    ]
    lines = ["\t".join(columns)]
    for pair in result.pairs:
        specificity = pair.specificity
        dimers = specificity.get("dimers") or {}
        lines.append("\t".join(str(value) for value in [
            specificity.get("rank"),
            specificity.get("score"),
            specificity.get("specificity_status"),
            specificity.get("specific_all_db"),
            specificity.get("specific_observed_all_db"),
            specificity.get("search_completeness"),
            specificity.get("search_complete_all_db"),
            ";".join(specificity.get("incomplete_databases", [])),
            _thermo_summary(specificity.get("per_db", [])),
            specificity.get("gel_distinguishable"),
            pair.forward,
            pair.reverse,
            pair.product_size,
            round(pair.tm_f, 1),
            round(pair.tm_r, 1),
            round(pair.gc_f, 1),
            round(pair.gc_r, 1),
            specificity.get("total_on_target"),
            specificity.get("total_off_target"),
            specificity.get("total_comigrating"),
            dimers.get("ok"),
            dimers.get("cross_dimer_dg"),
        ]))
    return "\n".join(lines)


def _completeness_warning(result: Dict, indent: str = "") -> List[str]:
    completeness = result.get("search_completeness", "complete")
    if completeness == "complete":
        return []
    hit_counts = result.get("raw_hits_per_primer", {})
    subject_counts = result.get("unique_subjects_per_primer", {})
    primer_states = result.get("primer_search_completeness", {})
    details = ", ".join(
        "%s=%s (subjects=%s, HSPs=%s)" % (
            primer, state, subject_counts.get(primer, "?"),
            hit_counts.get(primer, "?"))
        for primer, state in sorted(primer_states.items())
        if state != "complete"
    )
    recommendation = result.get("completeness_recommendation") or (
        "rerun with --exhaustive or a larger --max-target-seqs")
    return [
        "%sWARNING: specificity evidence is %s%s" % (
            indent, completeness, " (%s)" % details if details else ""),
        "%s         %s" % (indent, recommendation),
    ]


def _thermo_line(result: Dict, indent: str = "") -> str:
    status = result.get("thermo_status")
    if not status:
        status = "evaluated" if result.get("thermo_evaluated") else "not_evaluated"
    fasta = result.get("thermo_genome_fasta")
    association = result.get("thermo_genome_association")
    detail = ""
    if fasta:
        detail += "  FASTA=%s" % fasta
    if association:
        detail += "  association=%s" % association
    return "%sthermo: %s%s" % (indent, status, detail)


def to_text(result: PipelineResult, max_offtarget_rows: int = 6) -> str:
    output: List[str] = [
        "=" * 72,
        "primerblast-oss (design)  |  template: %s (%s bp)" % (
            result.template_id, result.template_len),
        "databases: %s" % ", ".join(
            database.split("/")[-1] for database in result.databases),
        "primer pairs returned: %s" % len(result.pairs),
        "=" * 72,
    ]
    for rank_number, pair in enumerate(result.pairs, 1):
        specificity = pair.specificity
        status = specificity.get("specificity_status")
        if status == "specific":
            verdict = "SPECIFIC"
        elif status == "indeterminate":
            verdict = "INDETERMINATE (clean only in returned hits)"
        elif specificity.get("gel_distinguishable"):
            verdict = "GEL-RESOLVABLE"
        else:
            verdict = "AMBIGUOUS"
        output.extend([
            "",
            "[%s] rank %s  score %s  %s" % (
                rank_number, specificity.get("rank"),
                specificity.get("score"), verdict),
            "    F  5'-%s-3'   Tm %.1f  GC %.0f%%  pos %s..%s" % (
                pair.forward, pair.tm_f, pair.gc_f,
                pair.left_start + 1, pair.left_3p + 1),
            "    R  5'-%s-3'   Tm %.1f  GC %.0f%%  pos %s..%s" % (
                pair.reverse, pair.tm_r, pair.gc_r,
                pair.right_3p + 1, pair.right_start + 1),
            "    product %s bp   on-target %s  off-target %s (co-migrating %s)" % (
                pair.product_size, specificity.get("total_on_target"),
                specificity.get("total_off_target"),
                specificity.get("total_comigrating")),
        ])
        dimers = specificity.get("dimers")
        if dimers:
            tag = "OK" if dimers.get("ok") else "%s concerning" % dimers.get("n_concerning")
            output.append(
                "    primer-dimer/hairpin: %s  worst ΔG %s  F×R ΔG %s kcal/mol"
                % (tag, dimers.get("worst_dg"), dimers.get("cross_dimer_dg")))
        for database in specificity.get("per_db", []):
            if database.get("specific") is True:
                tag = "OK (single exhaustive product)"
            elif database.get("specific") is None:
                tag = "INDETERMINATE (single observed product; search incomplete)"
            else:
                tag = "%s products; nearest off-target size gap %s bp" % (
                    database.get("n_products"),
                    database.get("nearest_offtarget_gap"))
            output.append("      %s: %s" % (
                database["db"].split("/")[-1], tag))
            output.append(_thermo_line(database, indent="         "))
            output.extend(_completeness_warning(database, indent="         "))
            if database.get("high_copy_primers"):
                output.append("         repeat-prone primer(s): %s" % ", ".join(
                    database["high_copy_primers"]))
            for amplicon in database.get("off_target", [])[:max_offtarget_rows]:
                tm_text = ""
                if getattr(amplicon, "fwd_tm", None) is not None:
                    tm_text = "  Tm %s/%s" % (
                        amplicon.fwd_tm, amplicon.rev_tm)
                output.append(
                    "         off: %s:%s-%s %sbp  %s  mm %s+%s%s" % (
                        amplicon.subject, amplicon.start, amplicon.end,
                        amplicon.size, amplicon.orientation,
                        amplicon.fwd_mismatch, amplicon.rev_mismatch, tm_text))
            extra = len(database.get("off_target", [])) - max_offtarget_rows
            if extra > 0:
                output.append("         ... %s more off-target products" % extra)
    return "\n".join(output)


def insilico_to_dict(results: List[Dict], primers: Dict[str, str]) -> dict:
    databases = []
    for result in results:
        databases.append({
            "db": result["db"],
            "sites_per_primer": result["sites_per_primer"],
            "thermo_status": result.get("thermo_status"),
            "thermo_evaluated": result.get("thermo_evaluated"),
            "thermo_genome_fasta": result.get("thermo_genome_fasta"),
            "thermo_genome_association": result.get("thermo_genome_association"),
            "viable_sites_per_primer": result.get("viable_sites_per_primer", {}),
            "search_completeness": result.get("search_completeness", "complete"),
            "search_complete": result.get("search_complete", True),
            "primer_search_completeness": result.get(
                "primer_search_completeness", {}),
            "completeness_recommendation": result.get(
                "completeness_recommendation"),
            "raw_hits_per_primer": result.get("raw_hits_per_primer", {}),
            "unique_subjects_per_primer": result.get(
                "unique_subjects_per_primer", {}),
            "near_blast_limit": result.get("near_blast_limit", []),
            "at_blast_limit": result.get("at_blast_limit", []),
            "high_copy_primers": result.get("high_copy_primers", []),
            "blast_limits": result.get("blast_limits", {}),
            "n_products": result["n_products"],
            "products": [_amp_dict(amplicon) for amplicon in result["products"]],
        })
    return {"mode": "check", "primers": primers, "results": databases}


def insilico_to_text(results: List[Dict], primers: Dict[str, str]) -> str:
    output: List[str] = ["=" * 72, "primerblast-oss (in-silico PCR)"]
    for name, sequence in primers.items():
        output.append("    %s: 5'-%s-3'  (%s nt)" % (
            name, sequence, len(sequence)))
    output.append("=" * 72)
    for result in results:
        output.extend([
            "",
            "# %s: %s predicted product(s)   sites/primer: %s" % (
                result["db"].split("/")[-1], result["n_products"],
                result["sites_per_primer"]),
            "    search completeness: %s" % result.get(
                "search_completeness", "complete"),
            _thermo_line(result, indent="    "),
        ])
        output.extend(_completeness_warning(result, indent="    "))
        if result.get("thermo_evaluated"):
            output.append(
                "    viable priming sites/primer: %s" % result.get(
                    "viable_sites_per_primer", {}))
        if not result["products"]:
            output.append("    (no products within the size window)")
            continue
        thermo = result.get("thermo_evaluated")
        header = "    size   subject:start-end            primers    mm      Δsize"
        output.append(header + ("   Tm(F/R)" if thermo else ""))
        for amplicon in result["products"]:
            gap = amplicon.__dict__.get("nearest_gap")
            gap_text = "-" if gap is None else str(gap)
            tm_text = ""
            if thermo:
                forward_tm = getattr(amplicon, "fwd_tm", None)
                reverse_tm = getattr(amplicon, "rev_tm", None)
                tm_text = "   %s/%s" % (
                    forward_tm if forward_tm is not None else "-",
                    reverse_tm if reverse_tm is not None else "-")
            output.append(
                "    %5s  %s:%s-%-12s  %-9s  %s+%-4s  %s%s" % (
                    amplicon.size, amplicon.subject, amplicon.start,
                    amplicon.end, amplicon.orientation,
                    amplicon.fwd_mismatch, amplicon.rev_mismatch,
                    gap_text, tm_text))
    return "\n".join(output)


def tiling_to_dict(tiles: List[Dict], template_id: str, region,
                   databases) -> dict:
    output_tiles = []
    for tile in tiles:
        pair = tile["pair"]
        output_tiles.append({
            "index": tile["index"],
            "covers": tile["covers"],
            "gap_to_prev": tile.get("gap_to_prev"),
            "forward": pair.forward,
            "reverse": pair.reverse,
            "product_size": pair.product_size,
            "left_pos": [pair.left_start + 1, pair.left_3p + 1],
            "right_pos": [pair.right_3p + 1, pair.right_start + 1],
            "tm_f": round(pair.tm_f, 1),
            "tm_r": round(pair.tm_r, 1),
            "specificity": _specificity_to_dict(pair.specificity),
        })
    return {
        "mode": "tile", "template_id": template_id, "region": region,
        "databases": databases, "n_tiles": len(tiles), "tiles": output_tiles,
    }


def assay_to_text(result: Dict) -> str:
    target = result["target"]
    output: List[str] = [
        "=" * 72,
        "primerblast-oss (assay)  |  %s  %s:%s-%s (%s) [%s]" % (
            target["name"], target["chrom"], target["start"], target["end"],
            target["strand"], target["source"]),
        "template %s bp  |  references: %s" % (
            result["template_len"], ", ".join(
                database.split("/")[-1] for database in result["databases"])),
        "primer pairs: %s" % result["n_pairs"],
        "=" * 72,
    ]
    for index, pair in enumerate(result["pairs"], 1):
        output.extend([
            "",
            "[%s] risk %s (score %s)  product %s bp" % (
                index, pair["risk"].upper(), pair.get("risk_score"),
                pair["product_size"]),
            "    specificity: %s; search %s" % (
                pair.get("specificity_status"),
                pair.get("search_completeness", "complete")),
            "    F 5'-%s-3'  Tm %s GC %s%%  @ %s-%s" % (
                pair["forward"], pair["tm_f"], pair["gc_f"],
                pair["forward_pos"][0], pair["forward_pos"][1]),
            "    R 5'-%s-3'  Tm %s GC %s%%  @ %s-%s" % (
                pair["reverse"], pair["tm_r"], pair["gc_r"],
                pair["reverse_pos"][0], pair["reverse_pos"][1]),
            "    off-target %s (F/R %s, F/F %s, R/R %s)" % (
                pair["n_off_target"], pair["n_fr_offtarget"],
                pair["n_ff"], pair["n_rr"]),
        ])
        conservation = pair.get("conservation", {})
        output.append("    conserved in %s/%s refs: %s" % (
            conservation.get("n_conserved"), conservation.get("n_refs"),
            ", ".join(pair.get("conserved_refs", []))))
        if pair.get("variant_in_primer"):
            suffix = " (3prime!)" if pair.get("variant_in_primer_3prime") else ""
            output.append("    variant under primer%s" % suffix)
        if pair.get("caps"):
            caps = pair["caps"]
            if caps.get("best_enzyme"):
                output.append("    CAPS: %s  ref %s vs alt %s" % (
                    caps["best_enzyme"], caps["allele_ref_fragments"],
                    caps["allele_alt_fragments"]))
            else:
                output.append("    CAPS: no distinguishing enzyme found")
        for database in pair.get("per_db_products", []):
            output.append("    [%s] %s" % (
                database["db"].split("/")[-1],
                _thermo_line(database)))
            if database.get("search_completeness") != "complete":
                output.append("    warning [%s]: search %s; %s" % (
                    database["db"].split("/")[-1],
                    database.get("search_completeness"),
                    database.get("completeness_recommendation")
                    or "increase the BLAST hit cap"))
        output.append("    reasons: %s" % "; ".join(
            pair.get("risk_reasons", [])))
    return "\n".join(output)


def tiling_to_text(tiles: List[Dict], template_id: str, region) -> str:
    output: List[str] = [
        "=" * 72,
        "primerblast-oss (tiling)  |  template: %s  region %s..%s" % (
            template_id, region[0] + 1, region[1] + 1),
        "amplicons covering the region: %s" % len(tiles),
        "=" * 72,
    ]
    if not tiles:
        output.append("  (no amplicons could be placed)")
        return "\n".join(output)
    covered_low = tiles[0]["covers"][0]
    covered_high = tiles[-1]["covers"][1]
    for tile in tiles:
        pair = tile["pair"]
        specificity = pair.specificity
        coverage = tile["covers"]
        output.extend([
            "",
            "  amplicon %s: covers %s..%s (%s bp, rank %s)" % (
                tile["index"], coverage[0] + 1, coverage[1] + 1,
                pair.product_size, specificity.get("rank", "?")),
            "      F 5'-%s-3'  Tm %.1f" % (pair.forward, pair.tm_f),
            "      R 5'-%s-3'  Tm %.1f" % (pair.reverse, pair.tm_r),
        ])
        for database in specificity.get("per_db", []):
            output.append("      %s" % _thermo_line(database))
        if specificity.get("specificity_status") == "indeterminate":
            output.append("      WARNING: specificity search incomplete")
    output.append("")
    output.append("  region coverage: %s..%s of requested %s..%s" % (
        covered_low + 1, covered_high + 1, region[0] + 1, region[1] + 1))
    return "\n".join(output)
