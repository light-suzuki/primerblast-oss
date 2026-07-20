"""Static CSV, order-sheet, BED, ASCII and HTML outputs."""
from __future__ import annotations

import csv
import html
import io
import json
from typing import Any, Dict, List, Sequence


def _get(mapping: Dict[str, Any], key: str, default: Any = "") -> Any:
    value = mapping.get(key, default)
    return default if value is None else value


def _num(value: Any) -> str:
    if value is None or value == "":
        return ""
    if isinstance(value, float):
        return "%.1f" % value
    return str(value)


def _bool_str(value: Any) -> str:
    if value is None or value == "":
        return ""
    return "yes" if bool(value) else "no"


def _tm_diff(pair: Dict[str, Any]) -> str:
    forward = pair.get("tm_f")
    reverse = pair.get("tm_r")
    if isinstance(forward, (int, float)) and isinstance(reverse, (int, float)):
        return "%.1f" % abs(float(forward) - float(reverse))
    return ""


def _synthesis_scale(length: int) -> str:
    if length <= 30:
        return "25 nmol (std desalt)"
    if length <= 45:
        return "100 nmol (std desalt)"
    return "250 nmol (PAGE)"


def _best_dcaps(pair: Dict[str, Any]) -> Dict[str, Any]:
    caps = pair.get("caps") or {}
    dcaps = caps.get("dcaps") or {}
    best = dcaps.get("best") or {}
    if caps.get("best_marker_type") == "dCAPS" and best.get("orderable"):
        return best
    return {}


def _orderable_pair(pair: Dict[str, Any]) -> Dict[str, Any]:
    """Return the actual oligos recommended for ordering.

    Natural CAPS and generic assays use the parent pair. A validated dCAPS result
    replaces the modified side with its engineered sequence and keeps the other
    parent primer.
    """
    derived = _best_dcaps(pair)
    if not derived:
        return {
            "name": _get(pair, "name", "pair"),
            "marker_type": _get(pair, "marker_type", "generic"),
            "enzyme": _get(pair, "caps_enzyme"),
            "forward": _get(pair, "forward"),
            "reverse": _get(pair, "reverse"),
            "tm_f": pair.get("tm_f"),
            "tm_r": pair.get("tm_r"),
            "gc_f": pair.get("gc_f"),
            "gc_r": pair.get("gc_r"),
            "product_size": pair.get("product_size"),
            "engineered_role": "",
            "engineered_mismatches": 0,
        }
    return {
        "name": "%s_dCAPS_%s" % (
            _get(pair, "name", "pair"), derived.get("enzyme", "enzyme")),
        "marker_type": "dCAPS",
        "enzyme": derived.get("enzyme"),
        "forward": derived.get("forward"),
        "reverse": derived.get("reverse"),
        "tm_f": derived.get("tm_f"),
        "tm_r": derived.get("tm_r"),
        "gc_f": derived.get("gc_f"),
        "gc_r": derived.get("gc_r"),
        "product_size": derived.get("product_size"),
        "engineered_role": derived.get("modified_primer_role"),
        "engineered_mismatches": derived.get("engineered_mismatches"),
    }


def to_csv_generic(rows: List[dict], columns: List[str]) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(columns)
    for row in rows:
        writer.writerow([
            "" if row.get(column) is None else row.get(column, "")
            for column in columns
        ])
    return buffer.getvalue()


_CSV_COLUMNS = [
    # Historical columns remain in their original order for compatibility.
    "name", "forward", "reverse", "product_size",
    "tm_f", "tm_r", "tm_diff", "gc_f", "gc_r",
    "n_off_target", "n_ff", "n_rr", "n_fr_offtarget",
    "tp5_mismatch_min", "snp_in_primer", "conserved_refs",
    "caps_enzyme", "gel_distinguishable", "cross_dimer_dg", "risk",
    # New fields are append-only.
    "marker_type", "enzyme", "engineered_role", "engineered_mismatches",
    "specificity_status", "specificity_status_all_db",
    "search_completeness", "variant_in_primer",
]


def pairs_to_csv(pairs: List[dict]) -> str:
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(_CSV_COLUMNS)
    for pair in pairs:
        orderable = _orderable_pair(pair)
        conserved = pair.get("conserved_refs", [])
        conserved_text = (
            ";".join(str(value) for value in conserved)
            if isinstance(conserved, (list, tuple)) else str(conserved)
        )
        variant_flag = pair.get(
            "variant_in_primer", pair.get("snp_in_primer"))
        writer.writerow([
            orderable["name"], orderable["forward"], orderable["reverse"],
            _num(orderable["product_size"]), _num(orderable["tm_f"]),
            _num(orderable["tm_r"]), _tm_diff(orderable),
            _num(orderable["gc_f"]), _num(orderable["gc_r"]),
            _num(pair.get("n_off_target")), _num(pair.get("n_ff")),
            _num(pair.get("n_rr")), _num(pair.get("n_fr_offtarget")),
            _num(pair.get("tp5_mismatch_min")),
            _bool_str(pair.get("snp_in_primer", variant_flag)),
            conserved_text,
            orderable["enzyme"] or pair.get("caps_enzyme"),
            _bool_str(pair.get("gel_distinguishable")),
            _num((pair.get("dimers") or {}).get("cross_dimer_dg")),
            pair.get("risk"), orderable["marker_type"], orderable["enzyme"],
            orderable["engineered_role"], orderable["engineered_mismatches"],
            pair.get("specificity_status"),
            pair.get("specificity_status_all_db"),
            pair.get("search_completeness"), _bool_str(variant_flag),
        ])
    return buffer.getvalue()


def order_table(pairs: List[dict]) -> str:
    header = (
        "Oligo name", "Sequence (5'->3')", "Len", "Tm", "Scale",
        "Marker", "Enzyme", "Engineered",
    )
    rows: List[Sequence[str]] = []
    for pair in pairs:
        orderable = _orderable_pair(pair)
        engineered_role = orderable.get("engineered_role") or ""
        mismatch_count = orderable.get("engineered_mismatches") or 0
        for role, sequence_key, tm_key in (
            ("F", "forward", "tm_f"), ("R", "reverse", "tm_r")
        ):
            sequence = str(orderable.get(sequence_key) or "")
            engineered = (
                "%s mismatch(es)" % mismatch_count
                if engineered_role == role else ""
            )
            rows.append((
                "%s_%s" % (orderable["name"], role),
                sequence,
                str(len(sequence)) if sequence else "",
                _num(orderable.get(tm_key)),
                _synthesis_scale(len(sequence)),
                str(orderable.get("marker_type") or ""),
                str(orderable.get("enzyme") or ""),
                engineered,
            ))
    widths = [len(column) for column in header]
    for row in rows:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], len(str(cell)))

    def format_row(cells: Sequence[str]) -> str:
        return "  ".join(
            str(cells[index]).ljust(widths[index])
            for index in range(len(header))
        ).rstrip()

    lines = [
        "Primer synthesis order sheet",
        "",
        format_row(header),
        "  ".join("-" * width for width in widths),
    ]
    lines.extend(format_row(row) for row in rows)
    lines.extend([
        "",
        "dCAPS rows contain the validated engineered primer sequence. "
        "Do not substitute the unmodified parent primer.",
        "Scales are suggestions; adjust to the supplier and application.",
    ])
    return "\n".join(lines) + "\n"


def products_to_bed(pairs: List[dict],
                    track_name: str = "primerblast_oss") -> str:
    lines = [
        'track name=%s description="predicted products"' % track_name]
    for pair in pairs:
        pair_name = str(_get(pair, "name", "?"))
        for product in pair.get("products", []) or []:
            chromosome = str(_get(product, "subject", ""))
            start, end = product.get("start"), product.get("end")
            if not chromosome or start is None or end is None:
                continue
            try:
                start_zero = int(start) - 1
                end_zero = int(end)
            except (TypeError, ValueError):
                continue
            orientation = str(_get(product, "orientation", ""))
            size = product.get("size")
            label = "%s:%s:%sbp" % (
                pair_name, orientation, size if size is not None else "")
            lines.append("%s\t%s\t%s\t%s\t.\t+" % (
                chromosome, start_zero, end_zero, label.replace(" ", "_")))
    return "\n".join(lines) + "\n"


def ascii_offtarget_map(pair: dict, width: int = 56) -> str:
    products = list(pair.get("products", []) or [])
    name = str(_get(pair, "name", "?"))
    if not products:
        return "%s: (no predicted products)\n" % name
    grouped: Dict[str, List[dict]] = {}
    for product in products:
        grouped.setdefault(str(_get(product, "subject", "?")), []).append(product)
    lines = ["Off-target map for %s" % name]
    for subject, subject_products in grouped.items():
        numeric_sizes = [
            product.get("size") for product in subject_products
            if isinstance(product.get("size"), (int, float))
        ]
        maximum = max(numeric_sizes) if numeric_sizes else 1
        for product in subject_products:
            orientation = str(_get(product, "orientation", "F/R")).upper()
            size = product.get("size")
            target = "intended" if product.get("on_target") else "off-target"
            span = int(round(float(size) / maximum * width)) if isinstance(
                size, (int, float)) else width
            span = max(4, span)
            parts = orientation.split("/")
            left = parts[0] if parts else "F"
            right = parts[1] if len(parts) > 1 else "R"
            left_head = ">" if left == "F" else "<"
            right_head = "<" if right == "R" else ">"
            lines.extend([
                "%s  (%s, %s)" % (subject, target, orientation),
                "  %s %s%s%s %s    %s bp" % (
                    left, left_head, "-" * span, right_head, right,
                    size if size is not None else ""),
            ])
    return "\n".join(lines) + "\n"


_RISK_COLORS = {"low": "#1a7f37", "medium": "#bf8700", "high": "#cf222e"}
_HTML_CSS = """
*{box-sizing:border-box}body{font-family:system-ui,sans-serif;color:#1f2328;
margin:0;padding:24px 32px;line-height:1.45;font-size:14px}h1{font-size:22px}
h2{font-size:17px;margin-top:28px;border-bottom:2px solid #d0d7de}
table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid
#d0d7de;padding:5px 8px;text-align:left;vertical-align:top}th{background:#f6f8fa}
code,.seq{font-family:ui-monospace,monospace;word-break:break-all}.pair{border:1px
solid #d0d7de;border-radius:8px;padding:12px 16px;margin:12px 0}.badge{display:inline-block;
padding:1px 9px;border-radius:10px;color:white;font-weight:600}.warning{background:#fff8c5;
border:1px solid #d4a72c;padding:8px}.ok{color:#1a7f37}.bad{color:#cf222e}.meta{color:#57606a}
@media print{body{padding:0}.pair{break-inside:avoid}}
"""


def _esc(value: Any) -> str:
    return html.escape("" if value is None else str(value))


def _risk_badge(risk: Any) -> str:
    normalized = str(risk or "").lower()
    color = _RISK_COLORS.get(normalized, "#57606a")
    return '<span class="badge" style="background:%s">%s</span>' % (
        color, _esc(risk or "n/a"))


def _summary_table(pairs: List[dict]) -> str:
    columns = [
        "Rank", "Name", "Marker", "Risk", "Product", "Specificity",
        "Search", "Off-targets",
    ]
    rows = []
    for index, pair in enumerate(pairs, 1):
        orderable = _orderable_pair(pair)
        rows.append([
            index,
            orderable["name"],
            "%s %s" % (orderable["marker_type"], orderable["enzyme"] or ""),
            _risk_badge(pair.get("risk")),
            orderable.get("product_size"),
            pair.get("specificity_status"),
            pair.get("search_completeness"),
            pair.get("n_off_target"),
        ])
    header = "".join("<th>%s</th>" % _esc(column) for column in columns)
    body = "".join(
        "<tr>%s</tr>" % "".join("<td>%s</td>" % cell for cell in [
            _esc(row[0]), _esc(row[1]), _esc(row[2]), row[3],
            _esc(row[4]), _esc(row[5]), _esc(row[6]), _esc(row[7]),
        ])
        for row in rows
    )
    return "<table><thead><tr>%s</tr></thead><tbody>%s</tbody></table>" % (
        header, body)


def _products_table(pair: dict) -> str:
    products = list(pair.get("products", []) or [])
    if not products:
        return '<p class="meta">No predicted products.</p>'
    rows = []
    for product in products:
        rows.append("<tr>%s</tr>" % "".join(
            "<td>%s</td>" % _esc(value) for value in [
                product.get("subject"), product.get("start"), product.get("end"),
                product.get("size"), product.get("orientation"),
                "intended" if product.get("on_target") else "off-target",
            ]))
    return ("<table><thead><tr><th>Subject</th><th>Start</th><th>End</th>"
            "<th>Size</th><th>Orientation</th><th>Target</th></tr></thead>"
            "<tbody>%s</tbody></table>" % "".join(rows))


def _dcaps_html(pair: dict) -> str:
    best = _best_dcaps(pair)
    if not best:
        return ""
    digest = best.get("digest") or {}
    changes = best.get("engineered_changes_template") or []
    return """
<div class="warning"><b>Validated dCAPS order candidate</b><br>
Enzyme: {enzyme}; modified primer: {role}; engineered mismatch(es): {count}<br>
Forward: <span class="seq">{forward}</span><br>
Reverse: <span class="seq">{reverse}</span><br>
Predicted digest: ref {ref_fragments} vs alt {alt_fragments}; status {status}<br>
Engineered template changes: <code>{changes}</code></div>
""".format(
        enzyme=_esc(best.get("enzyme")),
        role=_esc(best.get("modified_primer_role")),
        count=_esc(best.get("engineered_mismatches")),
        forward=_esc(best.get("forward")),
        reverse=_esc(best.get("reverse")),
        ref_fragments=_esc(digest.get("allele_a_fragments")),
        alt_fragments=_esc(digest.get("allele_b_fragments")),
        status=_esc(best.get("recommendation_status")),
        changes=_esc(json.dumps(changes, ensure_ascii=False)),
    )


def _pair_detail(pair: dict, rank: int) -> str:
    orderable = _orderable_pair(pair)
    conserved = pair.get("conserved_refs", [])
    reasons = "; ".join(pair.get("risk_reasons", []) or [])
    return """
<div class="pair"><h3>#{rank} · {name} {risk}</h3>
<p><b>Recommended marker:</b> {marker} {enzyme}</p>
<p>Forward: <span class="seq">{forward}</span><br>
Reverse: <span class="seq">{reverse}</span></p>
<p>Product {product} bp; Tm {tmf}/{tmr} °C; GC {gcf}/{gcr}%</p>
<p>Specificity: {specificity}; search: {search}; off-targets: {off}</p>
{dcaps}
<h4>Predicted parent-pair products</h4>{products}
<p><b>Conserved in:</b> {conserved}<br><b>Risk reasons:</b> {reasons}</p>
</div>
""".format(
        rank=rank,
        name=_esc(orderable["name"]),
        risk=_risk_badge(pair.get("risk")),
        marker=_esc(orderable["marker_type"]),
        enzyme=_esc(orderable["enzyme"]),
        forward=_esc(orderable["forward"]),
        reverse=_esc(orderable["reverse"]),
        product=_esc(orderable["product_size"]),
        tmf=_esc(_num(orderable["tm_f"])),
        tmr=_esc(_num(orderable["tm_r"])),
        gcf=_esc(_num(orderable["gc_f"])),
        gcr=_esc(_num(orderable["gc_r"])),
        specificity=_esc(pair.get("specificity_status")),
        search=_esc(pair.get("search_completeness")),
        off=_esc(pair.get("n_off_target")),
        dcaps=_dcaps_html(pair),
        products=_products_table(pair),
        conserved=_esc(", ".join(str(value) for value in conserved)),
        reasons=_esc(reasons),
    )


def html_report(context: dict) -> str:
    title = _esc(_get(context, "title", "primerblast_oss report"))
    template = _esc(_get(context, "template", ""))
    generated = _esc(_get(context, "generated", ""))
    databases = context.get("databases", []) or []
    pairs = list(context.get("pairs", []) or [])
    provenance = {
        "params": context.get("params", {}),
        "provenance": context.get("provenance", {}),
    }
    details = "\n".join(
        _pair_detail(pair, index) for index, pair in enumerate(pairs, 1))
    return """<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{title}</title>
<style>{css}</style></head><body><h1>{title}</h1>
<p class="meta"><b>Template:</b> {template}<br><b>Databases:</b> {databases}<br>
<b>Generated:</b> {generated}</p><h2>Summary</h2>{summary}
<h2>Primer pair details</h2>{details}<h2>Provenance</h2><pre>{provenance}</pre>
</body></html>""".format(
        title=title,
        css=_HTML_CSS,
        template=template or "n/a",
        databases=_esc(", ".join(str(value) for value in databases)) or "n/a",
        generated=generated or "n/a",
        summary=_summary_table(pairs),
        details=details,
        provenance=_esc(json.dumps(
            provenance, indent=2, ensure_ascii=False, default=str)),
    )


if __name__ == "__main__":
    sample = {
        "name": "P1", "forward": "ACGTACGT", "reverse": "TGCATGCA",
        "product_size": 100, "tm_f": 55.0, "tm_r": 55.0,
        "gc_f": 50.0, "gc_r": 50.0, "risk": "low",
        "specificity_status": "specific", "search_completeness": "complete",
        "products": [{"subject": "chr1", "start": 1, "end": 100,
                      "size": 100, "orientation": "F/R", "on_target": True}],
    }
    assert "P1_F" in order_table([sample])
    assert "chr1\t0\t100" in products_to_bed([sample])
    assert "primerblast_oss" in html_report({"pairs": [sample]})
    assert "forward" in pairs_to_csv([sample])
    print("All self-tests passed.")
