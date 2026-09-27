"""Sequencing / primer-walking orchestration built on the tiling engine."""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

M13_FORWARD_TAIL = "TGTAAAACGACGGCCAGT"
M13_REVERSE_TAIL = "CAGGAAACAGCTATGACC"


def local_region_for_template(template) -> Tuple[int, int]:
    """Return the target region as 0-based inclusive indices on a Template."""
    if template.anchor_strand == "-":
        first = template.anchor_coord - template.region.end
        last = template.anchor_coord - template.region.start
    else:
        first = template.region.start - template.anchor_coord
        last = template.region.end - template.anchor_coord
    return (min(first, last), max(first, last))


def coverage_summary(tiles: Sequence[Dict], region: Tuple[int, int]) -> Dict:
    """Summarize exact target coverage and uncovered intervals."""
    low, high = region
    if high < low:
        raise ValueError("sequencing region end must be >= start")
    spans = []
    for tile in tiles:
        start, end = tile["covers"]
        start, end = max(low, min(start, end)), min(high, max(start, end))
        if end >= start:
            spans.append((start, end))
    spans.sort()

    merged: List[List[int]] = []
    for start, end in spans:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    gaps: List[List[int]] = []
    cursor = low
    for start, end in merged:
        if start > cursor:
            gaps.append([cursor, start - 1])
        cursor = max(cursor, end + 1)
    if cursor <= high:
        gaps.append([cursor, high])

    target_bases = high - low + 1
    covered_bases = sum(end - start + 1 for start, end in merged)
    return {
        "target_bases": target_bases,
        "covered_bases": covered_bases,
        "coverage_fraction": (
            covered_bases / target_bases if target_bases else 0.0),
        "full_coverage": covered_bases == target_bases,
        "covered_intervals": merged,
        "gaps": gaps,
    }


def build_sequence_plan(
    tiles: Sequence[Dict],
    template_id: str,
    region: Tuple[int, int],
    databases: Sequence[str],
    *,
    genomic_template=None,
    requested_overlap: int = 0,
    amplicon_range: Optional[Tuple[int, int]] = None,
    m13_tails: bool = False,
    forward_tail: str = M13_FORWARD_TAIL,
    reverse_tail: str = M13_REVERSE_TAIL,
    template_sequence: Optional[str] = None,
) -> Dict:
    """Convert tiling output into an experiment-facing sequencing plan."""
    from .report import _specificity_to_dict

    records = []
    for index, tile in enumerate(tiles):
        pair = tile["pair"]
        gap_to_prev = tile.get("gap_to_prev")
        overlap_to_prev = max(0, int(gap_to_prev or 0))
        gap_from_prev = max(0, -int(gap_to_prev or 0))
        local_start, local_end = tile["covers"]

        record = {
            "index": tile["index"],
            "covers": [local_start, local_end],
            "product_size": pair.product_size,
            "forward": pair.forward,
            "reverse": pair.reverse,
            "forward_pos": [pair.left_start, pair.left_start + len(pair.forward) - 1],
            "reverse_pos": [pair.right_start - len(pair.reverse) + 1, pair.right_start],
            "order_forward": (
                forward_tail + pair.forward if m13_tails else pair.forward),
            "order_reverse": (
                reverse_tail + pair.reverse if m13_tails else pair.reverse),
            "tm_f": round(pair.tm_f, 1),
            "tm_r": round(pair.tm_r, 1),
            "overlap_to_prev": overlap_to_prev,
            "gap_from_prev": gap_from_prev,
            "specificity": _specificity_to_dict(pair.specificity),
        }
        if index + 1 < len(tiles):
            next_start = min(tiles[index + 1]["covers"])
            current_end = max(tile["covers"])
            record["overlap_to_next"] = max(
                0, current_end - next_start + 1)
            record["gap_to_next"] = max(
                0, next_start - current_end - 1)
        else:
            record["overlap_to_next"] = 0
            record["gap_to_next"] = 0

        if genomic_template is not None:
            g1 = genomic_template.to_genomic(local_start)
            g2 = genomic_template.to_genomic(local_end)
            record["genomic"] = {
                "chrom": genomic_template.region.chrom,
                "start": min(g1, g2),
                "end": max(g1, g2),
                "strand": genomic_template.region.strand,
            }
        records.append(record)

    plan = {
        "mode": "sequence",
        "template_id": template_id,
        "region": list(region),
        "databases": list(databases),
        "requested_overlap": requested_overlap,
        "amplicon_range": (
            list(amplicon_range) if amplicon_range is not None else None),
        "m13_tails": m13_tails,
        "tails": {
            "forward": forward_tail if m13_tails else "",
            "reverse": reverse_tail if m13_tails else "",
        },
        "coverage": coverage_summary(tiles, region),
        "n_amplicons": len(records),
        "amplicons": records,
    }
    if genomic_template is not None:
        plan["template"] = {
            "sequence": genomic_template.seq, "chrom": genomic_template.region.chrom,
            "start": genomic_template.ext_start, "end": genomic_template.ext_end,
            "anchor": genomic_template.anchor_coord, "strand": genomic_template.anchor_strand,
        }
        plan["target"] = {
            "name": genomic_template.region.name,
            "chrom": genomic_template.region.chrom,
            "start": genomic_template.region.start,
            "end": genomic_template.region.end,
            "strand": genomic_template.region.strand,
            "source": genomic_template.region.source,
            "flank": genomic_template.flank,
        }
    elif template_sequence is not None:
        plan["template"] = {"sequence": template_sequence, "anchor": 1, "strand": "+"}
    return plan


def plan_to_text(plan: Dict) -> str:
    coverage = plan["coverage"]
    lines = [
        "=" * 72,
        "primerblast-oss (sequencing) | %s" % plan["template_id"],
        "amplicons: %s | requested overlap: %s bp | coverage: %.1f%%" % (
            plan["n_amplicons"], plan["requested_overlap"],
            100.0 * coverage["coverage_fraction"]),
        "=" * 72,
    ]
    if plan.get("target"):
        target = plan["target"]
        lines.append("target: %s:%s-%s (%s)" % (
            target["chrom"], target["start"], target["end"],
            target["strand"]))
    if plan.get("amplicon_range"):
        lines.append("requested amplicon size: %s-%s bp" % tuple(
            plan["amplicon_range"]))
    if plan.get("m13_tails"):
        lines.append("M13 tails: enabled for ordered oligos only")
    if coverage["gaps"]:
        lines.append("WARNING: uncovered target intervals (template coords): %s" % (
            ", ".join("%s-%s" % (a + 1, b + 1)
                      for a, b in coverage["gaps"])))
    for amplicon in plan["amplicons"]:
        lines.extend([
            "",
            "[%s] %s bp | local %s-%s | overlap prev/next %s/%s bp | rank %s" % (
                amplicon["index"], amplicon["product_size"],
                amplicon["covers"][0] + 1, amplicon["covers"][1] + 1,
                amplicon["overlap_to_prev"], amplicon["overlap_to_next"],
                amplicon["specificity"].get("rank")),
            "    F anneal 5'-%s-3'" % amplicon["forward"],
            "    R anneal 5'-%s-3'" % amplicon["reverse"],
        ])
        if plan.get("m13_tails"):
            lines.extend([
                "    F order  5'-%s-3'" % amplicon["order_forward"],
                "    R order  5'-%s-3'" % amplicon["order_reverse"],
            ])
        if amplicon.get("genomic"):
            genomic = amplicon["genomic"]
            lines.append("    genomic %s:%s-%s (%s)" % (
                genomic["chrom"], genomic["start"], genomic["end"],
                genomic["strand"]))
    return "\n".join(lines)


def order_table(plan: Dict) -> str:
    rows = [("Name", "Sequence (5'->3')", "Role", "Amplicon")]
    for amplicon in plan["amplicons"]:
        prefix = "%s_A%s" % (plan["template_id"], amplicon["index"])
        rows.append((prefix + "_F", amplicon["order_forward"], "F",
                     str(amplicon["product_size"])))
        rows.append((prefix + "_R", amplicon["order_reverse"], "R",
                     str(amplicon["product_size"])))
    widths = [max(len(str(row[i])) for row in rows) for i in range(4)]
    return "\n".join(
        "  ".join(str(row[i]).ljust(widths[i]) for i in range(4)).rstrip()
        for row in rows
    ) + "\n"
