"""Render results as JSON, TSV, or human-readable text.

All renderers emit stable, GUI-friendly structures: primers carry sequence /
Tm / GC / position, and every predicted product carries subject, coordinates,
size, primer orientation, mismatches, and its nearest size gap to another
product (for gel-resolvability).
"""
from __future__ import annotations

import json
from dataclasses import asdict, is_dataclass
from typing import Dict, List

from .pipeline import PipelineResult


def _amp_dict(a) -> dict:
    d = asdict(a) if is_dataclass(a) else dict(a)
    if hasattr(a, "__dict__"):
        for key, val in a.__dict__.items():
            d.setdefault(key, val)
    d["orientation"] = f"{d['fwd_primer']}/{d['rev_primer']}"
    return d


# --------------------------------------------------------------------------- #
# design pipeline
# --------------------------------------------------------------------------- #
def _enrich_amplicon(a: dict) -> dict:
    a["orientation"] = f"{a['fwd_primer']}/{a['rev_primer']}"
    return a


def _specificity_to_dict(spec: dict) -> dict:
    d = dict(spec)
    per_db = []
    for db in spec.get("per_db", []):
        dbd = dict(db)
        for key in ("on_target", "off_target"):
            dbd[key] = [_amp_dict(a) for a in db.get(key, [])]
        per_db.append(dbd)
    d["per_db"] = per_db
    return d


def to_dict(result: PipelineResult) -> dict:
    pairs = []
    for p in result.pairs:
        d = asdict(p)
        d["specificity"] = _specificity_to_dict(p.specificity)
        pairs.append(d)
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


def to_tsv(result: PipelineResult) -> str:
    cols = [
        "rank", "score", "specific_all_db", "gel_distinguishable",
        "forward", "reverse", "product_size", "tm_f", "tm_r", "gc_f", "gc_r",
        "total_on_target", "total_off_target", "total_comigrating",
    ]
    lines = ["\t".join(cols)]
    for p in result.pairs:
        s = p.specificity
        lines.append("\t".join(str(x) for x in [
            s.get("rank"), s.get("score"), s.get("specific_all_db"),
            s.get("gel_distinguishable"),
            p.forward, p.reverse, p.product_size,
            round(p.tm_f, 1), round(p.tm_r, 1), round(p.gc_f, 1), round(p.gc_r, 1),
            s.get("total_on_target"), s.get("total_off_target"),
            s.get("total_comigrating"),
        ]))
    return "\n".join(lines)


def to_text(result: PipelineResult, max_offtarget_rows: int = 6) -> str:
    out: List[str] = []
    out.append("=" * 72)
    out.append(f"primerblast-oss (design)  |  template: {result.template_id} "
               f"({result.template_len} bp)")
    out.append(f"databases: {', '.join(db.split('/')[-1] for db in result.databases)}")
    out.append(f"primer pairs returned: {len(result.pairs)}")
    out.append("=" * 72)
    for rankn, p in enumerate(result.pairs, 1):
        s = p.specificity
        verdict = ("SPECIFIC" if s.get("specific_all_db")
                   else "GEL-RESOLVABLE" if s.get("gel_distinguishable")
                   else "AMBIGUOUS")
        out.append("")
        out.append(f"[{rankn}] rank {s.get('rank')}  score {s.get('score')}  {verdict}")
        out.append(f"    F  5'-{p.forward}-3'   Tm {p.tm_f:.1f}  GC {p.gc_f:.0f}%  "
                   f"pos {p.left_start+1}..{p.left_3p+1}")
        out.append(f"    R  5'-{p.reverse}-3'   Tm {p.tm_r:.1f}  GC {p.gc_r:.0f}%  "
                   f"pos {p.right_3p+1}..{p.right_start+1}")
        out.append(f"    product {p.product_size} bp   "
                   f"on-target {s.get('total_on_target')}  "
                   f"off-target {s.get('total_off_target')} "
                   f"(co-migrating {s.get('total_comigrating')})")
        for db in s.get("per_db", []):
            if db["specific"]:
                tag = "OK (single product)"
            else:
                tag = (f"{db['n_products']} products; nearest off-target size gap "
                       f"{db.get('nearest_offtarget_gap')} bp")
            out.append(f"      {db['db'].split('/')[-1]}: {tag}")
            if db.get("near_blast_limit"):
                hit_counts = db.get("raw_hits_per_primer", {})
                subjects = db.get("unique_subjects_per_primer", {})
                capped = ", ".join(
                    f"{p} subjects={subjects.get(p, '?')} hits={hit_counts.get(p, '?')}"
                    for p in db["near_blast_limit"])
                limit = db.get("blast_limits", {}).get("max_target_seqs")
                out.append(f"         warning: BLAST target list near limit ({capped}; "
                           f"max_target_seqs={limit}); raise --max-target-seqs for exhaustive checks")
            if db.get("high_copy_primers"):
                hit_counts = db.get("raw_hits_per_primer", {})
                copies = ", ".join(f"{p}={hit_counts.get(p, '?')}" for p in db["high_copy_primers"])
                threshold = db.get("blast_limits", {}).get("high_copy_hit_threshold")
                out.append(f"         warning: high-copy primer hit list ({copies}; "
                           f"threshold={threshold}); treat specificity as repeat-sensitive")
            for a in db["off_target"][:max_offtarget_rows]:
                out.append(
                    f"         off: {a.subject}:{a.start}-{a.end} "
                    f"{a.size}bp  {a.orientation}  mm {a.fwd_mismatch}+{a.rev_mismatch}"
                )
            extra = len(db["off_target"]) - max_offtarget_rows
            if extra > 0:
                out.append(f"         ... {extra} more off-target products")
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# in-silico PCR (check mode)
# --------------------------------------------------------------------------- #
def insilico_to_dict(results: List[Dict], primers: Dict[str, str]) -> dict:
    dbs = []
    for r in results:
        dbs.append({
            "db": r["db"],
            "sites_per_primer": r["sites_per_primer"],
            "raw_hits_per_primer": r.get("raw_hits_per_primer", {}),
            "unique_subjects_per_primer": r.get("unique_subjects_per_primer", {}),
            "near_blast_limit": r.get("near_blast_limit", []),
            "high_copy_primers": r.get("high_copy_primers", []),
            "blast_limits": r.get("blast_limits", {}),
            "n_products": r["n_products"],
            "products": [_amp_dict(a) for a in r["products"]],
        })
    return {"mode": "check", "primers": primers, "results": dbs}


def insilico_to_text(results: List[Dict], primers: Dict[str, str]) -> str:
    out: List[str] = []
    out.append("=" * 72)
    out.append("primerblast-oss (in-silico PCR)")
    for name, seq in primers.items():
        out.append(f"    {name}: 5'-{seq}-3'  ({len(seq)} nt)")
    out.append("=" * 72)
    for r in results:
        name = r["db"].split("/")[-1]
        out.append("")
        out.append(f"# {name}: {r['n_products']} predicted product(s)   "
                   f"sites/primer: {r['sites_per_primer']}")
        if r.get("near_blast_limit"):
            hit_counts = r.get("raw_hits_per_primer", {})
            subjects = r.get("unique_subjects_per_primer", {})
            capped = ", ".join(
                f"{p} subjects={subjects.get(p, '?')} hits={hit_counts.get(p, '?')}"
                for p in r["near_blast_limit"])
            limit = r.get("blast_limits", {}).get("max_target_seqs")
            out.append(f"    warning: BLAST target list near limit ({capped}; "
                       f"max_target_seqs={limit}); raise --max-target-seqs for exhaustive checks")
        if r.get("high_copy_primers"):
            hit_counts = r.get("raw_hits_per_primer", {})
            copies = ", ".join(f"{p}={hit_counts.get(p, '?')}" for p in r["high_copy_primers"])
            threshold = r.get("blast_limits", {}).get("high_copy_hit_threshold")
            out.append(f"    warning: high-copy primer hit list ({copies}; "
                       f"threshold={threshold}); treat specificity as repeat-sensitive")
        if not r["products"]:
            out.append("    (no products within the size window)")
            continue
        out.append("    size   subject:start-end            primers    mm      Δsize")
        for a in r["products"]:
            gap = a.__dict__.get("nearest_gap")
            gap_s = "-" if gap is None else str(gap)
            out.append(
                f"    {a.size:>5}  {a.subject}:{a.start}-{a.end:<12}  "
                f"{a.orientation:<9}  {a.fwd_mismatch}+{a.rev_mismatch:<4}  {gap_s}"
            )
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# tiling (whole-region walk)
# --------------------------------------------------------------------------- #
def tiling_to_dict(tiles: List[Dict], template_id: str, region, databases) -> dict:
    out_tiles = []
    for t in tiles:
        pair = t["pair"]
        out_tiles.append({
            "index": t["index"],
            "covers": t["covers"],
            "gap_to_prev": t.get("gap_to_prev"),
            "forward": pair.forward,
            "reverse": pair.reverse,
            "product_size": pair.product_size,
            "left_pos": [pair.left_start + 1, pair.left_3p + 1],
            "right_pos": [pair.right_3p + 1, pair.right_start + 1],
            "tm_f": round(pair.tm_f, 1), "tm_r": round(pair.tm_r, 1),
            "specificity": _specificity_to_dict(pair.specificity),
        })
    return {
        "mode": "tile",
        "template_id": template_id,
        "region": region,
        "databases": databases,
        "n_tiles": len(tiles),
        "tiles": out_tiles,
    }


def assay_to_text(result: Dict) -> str:
    t = result["target"]
    out: List[str] = []
    out.append("=" * 72)
    out.append(f"primerblast-oss (assay)  |  {t['name']}  "
               f"{t['chrom']}:{t['start']}-{t['end']} ({t['strand']}) [{t['source']}]")
    out.append(f"template {result['template_len']} bp  |  "
               f"references: {', '.join(db.split('/')[-1] for db in result['databases'])}")
    out.append(f"primer pairs: {result['n_pairs']}")
    out.append("=" * 72)
    for i, p in enumerate(result["pairs"], 1):
        out.append("")
        out.append(f"[{i}] risk {p['risk'].upper()} (score {p.get('risk_score')})  "
                   f"product {p['product_size']} bp")
        out.append(f"    F 5'-{p['forward']}-3'  Tm {p['tm_f']} GC {p['gc_f']}%  "
                   f"@ {p['left_pos'][0]}-{p['left_pos'][1]}")
        out.append(f"    R 5'-{p['reverse']}-3'  Tm {p['tm_r']} GC {p['gc_r']}%  "
                   f"@ {p['right_pos'][0]}-{p['right_pos'][1]}")
        out.append(f"    off-target {p['n_off_target']} "
                   f"(F/R {p['n_fr_offtarget']}, F/F {p['n_ff']}, R/R {p['n_rr']})  "
                   f"3'-5bp mm min {p.get('tp5_mismatch_min')}")
        cons = p.get("conservation", {})
        out.append(f"    conserved in {cons.get('n_conserved')}/{cons.get('n_refs')} refs: "
                   f"{', '.join(p.get('conserved_refs', []))}")
        if p.get("snp_in_primer"):
            out.append(f"    SNP under primer{' (3prime!)' if p.get('snp_in_primer_3prime') else ''}")
        if p.get("caps"):
            c = p["caps"]
            if c.get("best_enzyme"):
                out.append(f"    CAPS: {c['best_enzyme']}  "
                           f"ref {c['allele_ref_fragments']} vs alt {c['allele_alt_fragments']}")
            else:
                out.append("    CAPS: no distinguishing enzyme found")
        out.append(f"    reasons: {'; '.join(p.get('risk_reasons', []))}")
    return "\n".join(out)


def tiling_to_text(tiles: List[Dict], template_id: str, region) -> str:
    out: List[str] = []
    out.append("=" * 72)
    out.append(f"primerblast-oss (tiling)  |  template: {template_id}  "
               f"region {region[0]+1}..{region[1]+1}")
    out.append(f"amplicons covering the region: {len(tiles)}")
    out.append("=" * 72)
    if not tiles:
        out.append("  (no amplicons could be placed - relax design/product-size options)")
        return "\n".join(out)
    covered_lo = tiles[0]["covers"][0]
    covered_hi = tiles[-1]["covers"][1]
    for t in tiles:
        p = t["pair"]
        s = p.specificity
        rank = s.get("rank", "?") if s else "?"
        c = t["covers"]
        ov = "" if t.get("gap_to_prev") is None else f"  overlap/gap {t['gap_to_prev']} bp"
        out.append("")
        out.append(f"  amplicon {t['index']}: covers {c[0]+1}..{c[1]+1} "
                   f"({p.product_size} bp, rank {rank}){ov}")
        out.append(f"      F 5'-{p.forward}-3'  Tm {p.tm_f:.1f}")
        out.append(f"      R 5'-{p.reverse}-3'  Tm {p.tm_r:.1f}")
    out.append("")
    out.append(f"  region coverage: {covered_lo+1}..{covered_hi+1} "
               f"of requested {region[0]+1}..{region[1]+1}")
    return "\n".join(out)
