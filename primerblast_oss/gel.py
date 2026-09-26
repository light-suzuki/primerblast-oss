"""Gel-aware scoring and virtual-gel rendering for CAPS/dCAPS assays.

The scoring model is intentionally heuristic and interpretable. It is used to
rank experimentally convenient digest patterns, not to estimate a probability
of PCR or electrophoresis success.
"""
from __future__ import annotations

import html
import math
from collections import Counter
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


LADDER_PRESETS: Dict[str, List[int]] = {
    # Generic presets. Users can supply exact manufacturer bands with
    # --ladder-bands when needed.
    "100bp": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500],
    "1kb": [500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000],
}

# Approximate agarose separation ranges (bp), curated from the Thermo Fisher
# nucleic-acid electrophoresis workflow table. The preference order favors
# common concentrations before moving to unusually dense gels.
GEL_RANGES: Dict[float, Tuple[int, int]] = {
    0.5: (2000, 50000),
    0.6: (1000, 20000),
    0.7: (800, 12000),
    0.8: (800, 10000),
    0.9: (600, 10000),
    1.0: (400, 8000),
    1.2: (300, 7000),
    1.5: (200, 3000),
    2.0: (100, 2000),
    3.0: (25, 1000),
    4.0: (10, 500),
    5.0: (10, 300),
}
_GEL_PREFERENCE = [2.0, 1.5, 1.2, 1.0, 3.0, 0.8, 0.7, 0.6, 0.5, 4.0, 5.0]


def _multiset_difference(first: Sequence[int], second: Sequence[int]) -> List[int]:
    remaining = Counter(int(value) for value in first)
    remaining.subtract(Counter(int(value) for value in second))
    output: List[int] = []
    for size, count in remaining.items():
        if count > 0:
            output.extend([size] * count)
    return sorted(output, reverse=True)


def heterozygote_bands(
    fragments_a: Sequence[int], fragments_b: Sequence[int]
) -> List[Dict[str, int]]:
    """Return distinct AB band sizes with multiplicity/intensity proxy."""
    counts = Counter(int(value) for value in fragments_a)
    counts.update(int(value) for value in fragments_b)
    return [
        {"size": size, "copies": counts[size]}
        for size in sorted(counts, reverse=True)
    ]


def recommend_ladder(
    fragments: Sequence[int],
    ladder: str = "auto",
    custom_bands: Optional[Sequence[int]] = None,
) -> Tuple[str, List[int]]:
    sizes = [int(value) for value in fragments if int(value) > 0]
    if ladder == "custom":
        bands = sorted({int(value) for value in (custom_bands or []) if int(value) > 0})
        if not bands:
            raise ValueError("custom ladder requires at least one positive band")
        return "custom", bands
    if ladder != "auto":
        if ladder not in LADDER_PRESETS:
            raise ValueError("unknown ladder %r" % ladder)
        return ladder, list(LADDER_PRESETS[ladder])
    if sizes and min(sizes) >= 100 and max(sizes) <= 1500:
        return "100bp", list(LADDER_PRESETS["100bp"])
    return "1kb", list(LADDER_PRESETS["1kb"])


def recommend_gel_percent(
    fragments: Sequence[int], requested: Optional[float] = None
) -> float:
    sizes = [int(value) for value in fragments if int(value) > 0]
    if requested is not None:
        if requested not in GEL_RANGES:
            raise ValueError(
                "unsupported gel percentage %s; choose one of %s"
                % (requested, ", ".join(str(x) for x in sorted(GEL_RANGES)))
            )
        return requested
    if not sizes:
        return 2.0
    low, high = min(sizes), max(sizes)
    for percent in _GEL_PREFERENCE:
        lower, upper = GEL_RANGES[percent]
        if lower <= low and high <= upper:
            return percent
    # No single preset spans every fragment. Choose the preset that contains
    # the largest fraction, breaking ties toward the common lower percentage.
    return max(
        _GEL_PREFERENCE,
        key=lambda percent: (
            sum(GEL_RANGES[percent][0] <= size <= GEL_RANGES[percent][1]
                for size in sizes),
            -_GEL_PREFERENCE.index(percent),
        ),
    )


def _minimum_log_gap(sizes: Sequence[int]) -> Optional[float]:
    unique = sorted({int(value) for value in sizes if int(value) > 0})
    if len(unique) < 2:
        return None
    return min(
        abs(math.log10(high) - math.log10(low))
        for low, high in zip(unique, unique[1:])
    )


def _log_gap(a: int, b: int) -> float:
    if a <= 0 or b <= 0:
        return 0.0
    return abs(math.log10(a) - math.log10(b))


def _unique_band_separation(
    source: Sequence[int], other: Sequence[int]
) -> Optional[float]:
    """Best unique-band separation of one lane relative to another lane."""
    candidates = []
    for size in sorted({int(value) for value in source if int(value) > 0}):
        if not other:
            candidates.append(float("inf"))
            continue
        candidates.append(min(_log_gap(size, int(value)) for value in other))
    return max(candidates) if candidates else None


def genotype_pattern_discrimination(
    fragments_a: Sequence[int],
    fragments_b: Sequence[int],
    *,
    background_sizes: Sequence[int] = (),
    min_log10_gap: float = 0.04,
) -> Dict:
    """Test whether AA, AB and BB retain a visibly distinct band pattern.

    Background bands are added to every genotype lane. This models a
    non-specific PCR product that is present regardless of the target allele:
    it is acceptable unless it erases the diagnostic difference between two
    genotype lanes. Band intensity/copy-number differences are not used for the
    pass/fail decision because they are less robust than size differences.
    """
    background = [int(value) for value in background_sizes if int(value) > 0]
    aa = sorted(set([int(value) for value in fragments_a] + background), reverse=True)
    bb = sorted(set([int(value) for value in fragments_b] + background), reverse=True)
    ab = sorted(set(
        [int(value) for value in fragments_a]
        + [int(value) for value in fragments_b]
        + background
    ), reverse=True)
    lanes = {"AA": aa, "AB": ab, "BB": bb}
    comparisons = {}
    pair_names = (("AA", "AB"), ("AA", "BB"), ("AB", "BB"))
    for left_name, right_name in pair_names:
        left = lanes[left_name]
        right = lanes[right_name]
        left_gap = _unique_band_separation(left, right)
        right_gap = _unique_band_separation(right, left)
        usable_gaps = [
            value for value in (left_gap, right_gap)
            if value is not None and not math.isinf(value)
        ]
        if (left_gap is not None and math.isinf(left_gap)) or (
                right_gap is not None and math.isinf(right_gap)):
            best_gap = float("inf")
        else:
            best_gap = max(usable_gaps) if usable_gaps else 0.0
        distinguishable = best_gap >= min_log10_gap
        comparisons["%s_vs_%s" % (left_name, right_name)] = {
            "distinguishable": distinguishable,
            "best_unique_log10_gap": (
                None if math.isinf(best_gap) else round(best_gap, 5)),
        }

    overall = all(
        item["distinguishable"] for item in comparisons.values())
    finite = [
        item["best_unique_log10_gap"]
        for item in comparisons.values()
        if item["best_unique_log10_gap"] is not None
    ]
    return {
        "distinguishable": overall,
        "min_required_log10_gap": min_log10_gap,
        "weakest_pair_log10_gap": min(finite) if finite else None,
        "background_sizes": sorted(set(background), reverse=True),
        "lanes": lanes,
        "comparisons": comparisons,
    }


def assess_background_across_databases(
    gel_analysis: Dict,
    per_db_products: Sequence[Dict],
    *,
    min_log10_gap: float = 0.04,
) -> Dict:
    """Overlay non-target product sizes on CAPS genotype patterns per database.

    This intentionally models off-target products by their predicted undigested
    amplicon size. It is conservative about band overlap but does not yet claim
    to simulate restriction digestion of every off-target amplicon.
    """
    lane_a = [
        int(entry["size"]) for entry in gel_analysis.get("lanes", {}).get("AA", [])
    ]
    lane_b = [
        int(entry["size"]) for entry in gel_analysis.get("lanes", {}).get("BB", [])
    ]
    per_db = []
    for database in per_db_products:
        backgrounds = sorted({
            int(product["size"])
            for product in database.get("products", [])
            if not product.get("on_target") and product.get("size") is not None
        }, reverse=True)
        pattern = genotype_pattern_discrimination(
            lane_a,
            lane_b,
            background_sizes=backgrounds,
            min_log10_gap=min_log10_gap,
        )
        per_db.append({
            "db": database.get("db"),
            "n_off_target": database.get("n_off_target", len(backgrounds)),
            "background_sizes": backgrounds,
            "distinguishable": pattern["distinguishable"],
            "pattern": pattern,
        })
    all_ok = bool(per_db) and all(item["distinguishable"] for item in per_db)
    return {
        "model": "undigested_offtarget_amplicon_sizes",
        "all_databases_distinguishable": all_ok,
        "min_log10_gap": min_log10_gap,
        "per_db": per_db,
        "note": (
            "Off-target PCR products are overlaid by predicted amplicon size. "
            "Restriction digestion of off-target products is not simulated."
        ),
    }


def analyze_digest_patterns(
    fragments_a: Sequence[int],
    fragments_b: Sequence[int],
    *,
    ladder: str = "auto",
    custom_ladder_bands: Optional[Sequence[int]] = None,
    gel_percent: Optional[float] = None,
    min_visible_bp: int = 50,
) -> Dict:
    """Score genotype band patterns using fragment size, ladder and gel range."""
    a = [int(value) for value in fragments_a]
    b = [int(value) for value in fragments_b]
    diagnostic_a = _multiset_difference(a, b)
    diagnostic_b = _multiset_difference(b, a)
    ab = heterozygote_bands(a, b)
    combined_sizes = [entry["size"] for entry in ab]
    ladder_name, ladder_bands = recommend_ladder(
        diagnostic_a + diagnostic_b or combined_sizes,
        ladder=ladder,
        custom_bands=custom_ladder_bands,
    )
    percent = recommend_gel_percent(combined_sizes, requested=gel_percent)
    gel_low, gel_high = GEL_RANGES[percent]

    diagnostic = diagnostic_a + diagnostic_b
    genotype_discrimination = genotype_pattern_discrimination(
        a, b, min_log10_gap=0.04)
    log_gap = _minimum_log_gap(combined_sizes)
    separation_score = 100.0 if log_gap is None else min(100.0, 1000.0 * log_gap)
    ladder_score = (
        100.0 * sum(min(ladder_bands) <= size <= max(ladder_bands)
                    for size in diagnostic) / len(diagnostic)
        if diagnostic else 0.0
    )
    gel_range_score = (
        100.0 * sum(gel_low <= size <= gel_high for size in diagnostic)
        / len(diagnostic)
        if diagnostic else 0.0
    )
    visibility_score = (
        100.0 * sum(size >= min_visible_bp for size in diagnostic)
        / len(diagnostic)
        if diagnostic else 0.0
    )
    complexity_penalty = max(0, len(set(combined_sizes)) - 6) * 3.0
    score = max(0.0, min(
        100.0,
        0.55 * separation_score
        + 0.20 * ladder_score
        + 0.20 * gel_range_score
        + 0.05 * visibility_score
        - complexity_penalty,
    ))

    if score >= 85:
        rating = "excellent"
    elif score >= 70:
        rating = "good"
    elif score >= 50:
        rating = "marginal"
    else:
        rating = "poor"

    reasons = [
        "minimum log10 band gap=%s" % (
            "n/a" if log_gap is None else "%.4f" % log_gap),
        "ladder coverage %.0f%%" % ladder_score,
        "gel-range coverage %.0f%%" % gel_range_score,
    ]
    tiny = [size for size in diagnostic if size < min_visible_bp]
    outside_gel = [
        size for size in diagnostic if not (gel_low <= size <= gel_high)]
    outside_ladder = [
        size for size in diagnostic
        if not (min(ladder_bands) <= size <= max(ladder_bands))]
    if tiny:
        reasons.append(
            "diagnostic fragment(s) below %s bp: %s"
            % (min_visible_bp, ",".join(str(value) for value in sorted(tiny))))
    if outside_gel:
        reasons.append(
            "diagnostic fragment(s) outside recommended gel range: %s"
            % ",".join(str(value) for value in sorted(outside_gel)))
    if outside_ladder:
        reasons.append(
            "diagnostic fragment(s) outside ladder span: %s"
            % ",".join(str(value) for value in sorted(outside_ladder)))
    if complexity_penalty:
        reasons.append("complex multi-band pattern penalty %.0f" % complexity_penalty)

    lanes = {
        "M": [{"size": size, "copies": 1} for size in sorted(ladder_bands, reverse=True)],
        "AA": [{"size": size, "copies": count}
               for size, count in sorted(Counter(a).items(), reverse=True)],
        "AB": ab,
        "BB": [{"size": size, "copies": count}
               for size, count in sorted(Counter(b).items(), reverse=True)],
    }
    return {
        "score": round(score, 1),
        "rating": rating,
        "diagnostic_a": diagnostic_a,
        "diagnostic_b": diagnostic_b,
        "heterozygote_bands": ab,
        "genotype_discrimination": genotype_discrimination,
        "range_policy": "soft_warning_only",
        "tiny_diagnostic_fragments": tiny,
        "outside_gel_range_fragments": outside_gel,
        "outside_ladder_span_fragments": outside_ladder,
        "minimum_log10_gap": (
            None if log_gap is None else round(log_gap, 5)),
        "separation_score": round(separation_score, 1),
        "ladder": ladder_name,
        "ladder_bands": ladder_bands,
        "ladder_coverage_score": round(ladder_score, 1),
        "gel_percent": percent,
        "gel_range_bp": [gel_low, gel_high],
        "gel_range_score": round(gel_range_score, 1),
        "visibility_score": round(visibility_score, 1),
        "reasons": reasons,
        "lanes": lanes,
    }


def virtual_gel_svg(
    analysis: Dict,
    width: int = 520,
    height: int = 440,
) -> str:
    """Render deterministic M/AA/AB/BB lanes as a dependency-free SVG."""
    lanes = analysis["lanes"]
    all_sizes = [
        int(band["size"])
        for bands in lanes.values()
        for band in bands
        if int(band["size"]) > 0
    ]
    if not all_sizes:
        all_sizes = [100, 1000]
    min_size = max(1, min(all_sizes))
    max_size = max(all_sizes)
    if max_size == min_size:
        max_size = min_size * 10

    margin_top, margin_bottom = 35, 35
    plot_h = height - margin_top - margin_bottom
    labels = list(lanes.keys())
    lane_w = (width - 80) / len(labels)

    def y_for(size: int) -> float:
        top = math.log10(max_size)
        bottom = math.log10(min_size)
        fraction = (top - math.log10(max(1, size))) / max(1e-9, top - bottom)
        return margin_top + fraction * plot_h

    chunks = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" '
        'viewBox="0 0 %s %s">' % (width, height, width, height),
        '<rect width="100%" height="100%" fill="white"/>',
        '<rect x="55" y="%s" width="%s" height="%s" fill="#eceff1" '
        'stroke="#9aa0a6"/>' % (margin_top, width - 75, plot_h),
    ]
    for index, label in enumerate(labels):
        center = 70 + lane_w * (index + 0.5)
        chunks.append(
            '<text x="%.1f" y="20" text-anchor="middle" '
            'font-family="sans-serif" font-size="13">%s</text>'
            % (center, html.escape(label)))
        for band in lanes[label]:
            size = int(band["size"])
            copies = max(1, int(band.get("copies", 1)))
            y = y_for(size)
            opacity = min(1.0, 0.55 + 0.15 * (copies - 1))
            chunks.append(
                '<rect x="%.1f" y="%.1f" width="%.1f" height="4" '
                'rx="1" fill="#263238" opacity="%.2f"/>'
                % (center - lane_w * 0.30, y - 2, lane_w * 0.60, opacity))
            if label == "M":
                chunks.append(
                    '<text x="48" y="%.1f" text-anchor="end" '
                    'font-family="sans-serif" font-size="9">%s</text>'
                    % (y + 3, size))
    chunks.append(
        '<text x="%s" y="%s" text-anchor="middle" font-family="sans-serif" '
        'font-size="11">%.1f%% agarose · %s ladder · gel score %.1f (%s)</text>'
        % (width / 2, height - 8, analysis["gel_percent"],
           html.escape(str(analysis["ladder"])), analysis["score"],
           html.escape(str(analysis["rating"]))))
    chunks.append("</svg>")
    return "\n".join(chunks)


def best_analysis_from_assay(result: Dict) -> Optional[Dict]:
    """Find the first ranked CAPS/dCAPS gel analysis in an assay result."""
    for pair in result.get("pairs", []):
        caps = pair.get("caps") or {}
        if caps.get("best_marker_type") == "dCAPS":
            best = ((caps.get("dcaps") or {}).get("best") or {})
            digest = best.get("digest") or {}
            if digest.get("gel_analysis"):
                return digest["gel_analysis"]
        best_result = caps.get("best_result") or {}
        if best_result.get("gel_analysis"):
            return best_result["gel_analysis"]
    return None
