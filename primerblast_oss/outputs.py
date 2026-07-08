"""Self-contained output-formatting for primerblast_oss.

This module renders primer-pair results into shareable, static artefacts:
CSV tables, a plain-text oligo order sheet, BED tracks, ASCII off-target maps,
and a printable self-contained HTML report.

Design constraints (intentional):
    * Standard library ONLY (csv, io, html, json).
    * No imports from any other primerblast_oss module -> fully decoupled.
    * Every function consumes PLAIN dicts/lists following the documented schema
      below, so callers can build inputs without importing package types.
    * Nothing here is interactive; outputs are strings (or files) a user can
      save, share, or print to PDF.

Input schema -- a "primer pair record" dict::

    {
      "name": "P1", "forward": "ACGT...", "reverse": "TTGC...",
      "product_size": 368, "tm_f": 60.0, "tm_r": 60.1,
      "gc_f": 55.0, "gc_r": 45.0,
      "left_pos": [start, end], "right_pos": [start, end],   # 1-based
      "risk": "low",                                         # low|medium|high
      "n_off_target": 0, "n_ff": 0, "n_rr": 0, "n_fr_offtarget": 0,
      "tp5_mismatch_min": 0,        # min 3'-5bp mismatches among off-targets
      "snp_in_primer": false, "conserved_refs": ["cameor_v2", "ZW6"],
      "caps_enzyme": "EcoRI", "gel_distinguishable": true,
      "products": [
         {"subject": "chr1", "start": 85338374, "end": 85338741,
          "size": 368, "orientation": "F/R", "on_target": true},
         ...
      ]
    }

All getters tolerate missing keys; absent values render blank / neutral.
"""

from __future__ import annotations

import csv
import html
import io
import json
from typing import Any, Dict, List, Optional, Sequence


# --------------------------------------------------------------------------- #
# Small internal helpers
# --------------------------------------------------------------------------- #

def _get(d: Dict[str, Any], key: str, default: Any = "") -> Any:
    """Return ``d[key]`` or *default* when the key is missing or ``None``."""
    val = d.get(key, default)
    return default if val is None else val


def _num(val: Any) -> str:
    """Render a numeric value compactly, blank when missing/None."""
    if val is None or val == "":
        return ""
    if isinstance(val, float):
        # Trim trailing zeros but keep at least one decimal for Tm-like values.
        return f"{val:.1f}"
    return str(val)


def _tm_diff(pair: Dict[str, Any]) -> str:
    """Absolute Tm difference between the two primers, blank if unavailable."""
    tm_f = pair.get("tm_f")
    tm_r = pair.get("tm_r")
    if isinstance(tm_f, (int, float)) and isinstance(tm_r, (int, float)):
        return f"{abs(float(tm_f) - float(tm_r)):.1f}"
    return ""


def _bool_str(val: Any) -> str:
    """Render a boolean-ish flag as ``yes``/``no``/blank."""
    if val is None or val == "":
        return ""
    return "yes" if bool(val) else "no"


def _synthesis_scale(length: int) -> str:
    """Suggest a standard oligo synthesis scale from the oligo length.

    A pragmatic heuristic for ordering standard desalted primers; longer
    oligos typically move to a larger scale for adequate yield.
    """
    if length <= 0:
        return "25 nmol (std desalt)"
    if length <= 30:
        return "25 nmol (std desalt)"
    if length <= 45:
        return "100 nmol (std desalt)"
    return "250 nmol (PAGE)"


# --------------------------------------------------------------------------- #
# Generic CSV helper
# --------------------------------------------------------------------------- #

def to_csv_generic(rows: List[dict], columns: List[str]) -> str:
    """Render *rows* as CSV using exactly *columns* (in order).

    Missing keys become empty cells. Values are stringified as-is. Uses the
    ``csv`` module so quoting/escaping is always correct.
    """
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(columns)
    for row in rows:
        writer.writerow(["" if row.get(c) is None else row.get(c, "") for c in columns])
    return buf.getvalue()


# --------------------------------------------------------------------------- #
# Experimenter CSV
# --------------------------------------------------------------------------- #

# Column header covering the fields an experimenter cares about.
_CSV_COLUMNS: List[str] = [
    "name", "forward", "reverse", "product_size",
    "tm_f", "tm_r", "tm_diff",
    "gc_f", "gc_r",
    "n_off_target", "n_ff", "n_rr", "n_fr_offtarget",
    "tp5_mismatch_min", "snp_in_primer",
    "conserved_refs", "caps_enzyme", "gel_distinguishable", "risk",
]


def pairs_to_csv(pairs: List[dict]) -> str:
    """Return a CSV string of all primer pairs (one row per pair).

    Tolerates missing keys (rendered blank). ``conserved_refs`` is joined by
    ``';'`` so the CSV stays single-cell per field.
    """
    buf = io.StringIO()
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(_CSV_COLUMNS)

    for pair in pairs:
        conserved = _get(pair, "conserved_refs", [])
        if isinstance(conserved, (list, tuple)):
            conserved_str = ";".join(str(c) for c in conserved)
        else:
            conserved_str = str(conserved)

        writer.writerow([
            _get(pair, "name"),
            _get(pair, "forward"),
            _get(pair, "reverse"),
            _num(pair.get("product_size")),
            _num(pair.get("tm_f")),
            _num(pair.get("tm_r")),
            _tm_diff(pair),
            _num(pair.get("gc_f")),
            _num(pair.get("gc_r")),
            _num(pair.get("n_off_target")),
            _num(pair.get("n_ff")),
            _num(pair.get("n_rr")),
            _num(pair.get("n_fr_offtarget")),
            _num(pair.get("tp5_mismatch_min")),
            _bool_str(pair.get("snp_in_primer")),
            conserved_str,
            _get(pair, "caps_enzyme"),
            _bool_str(pair.get("gel_distinguishable")),
            _get(pair, "risk"),
        ])
    return buf.getvalue()


# --------------------------------------------------------------------------- #
# Oligo order sheet (plain text)
# --------------------------------------------------------------------------- #

def order_table(pairs: List[dict]) -> str:
    """Return a plain-text oligo ORDER sheet, one row per primer.

    Forward and reverse primers are listed on separate rows (``P1_F``,
    ``P1_R``) with sequence 5'->3', length, Tm, and a suggested synthesis
    scale, in nicely aligned columns.
    """
    header = ("Oligo name", "Sequence (5'->3')", "Len", "Tm", "Scale")

    # Build the row data first so we can size columns to content.
    rows: List[Sequence[str]] = []
    for pair in pairs:
        name = str(_get(pair, "name", "?"))
        for suffix, seq_key, tm_key in (("F", "forward", "tm_f"),
                                        ("R", "reverse", "tm_r")):
            seq = str(_get(pair, seq_key, ""))
            length = len(seq)
            tm = pair.get(tm_key)
            rows.append((
                f"{name}_{suffix}",
                seq,
                str(length) if length else "",
                _num(tm),
                _synthesis_scale(length),
            ))

    # Compute per-column widths (header vs. data).
    ncols = len(header)
    widths = [len(header[i]) for i in range(ncols)]
    for row in rows:
        for i in range(ncols):
            widths[i] = max(widths[i], len(row[i]))

    def fmt(cells: Sequence[str]) -> str:
        return "  ".join(str(cells[i]).ljust(widths[i]) for i in range(ncols)).rstrip()

    lines: List[str] = []
    lines.append("Primer synthesis order sheet")
    lines.append("")
    lines.append(fmt(header))
    lines.append("  ".join("-" * widths[i] for i in range(ncols)))
    for row in rows:
        lines.append(fmt(row))
    lines.append("")
    lines.append("Note: scales are suggestions for standard desalted oligos; "
                 "adjust to your supplier and application.")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# BED track of predicted products
# --------------------------------------------------------------------------- #

def products_to_bed(pairs: List[dict], track_name: str = "primerblast_oss") -> str:
    """Return a BED track (0-based half-open) of all predicted products.

    One BED line per predicted product across all pairs::

        chrom  start-1  end  name=<pair>:<orientation>:<size>bp  .  +

    A ``track name=...`` header line is emitted first.
    """
    lines: List[str] = [f'track name={track_name} description="predicted products"']

    for pair in pairs:
        pair_name = str(_get(pair, "name", "?"))
        for prod in _get(pair, "products", []) or []:
            chrom = str(_get(prod, "subject", ""))
            start = prod.get("start")
            end = prod.get("end")
            if chrom == "" or start is None or end is None:
                continue  # skip incomplete product records
            try:
                start0 = int(start) - 1  # BED is 0-based half-open
                end0 = int(end)
            except (TypeError, ValueError):
                continue
            orientation = str(_get(prod, "orientation", ""))
            size = prod.get("size")
            size_str = f"{int(size)}bp" if isinstance(size, (int, float)) else ""
            # BED name has no spaces; build a compact descriptive label.
            name = f"{pair_name}:{orientation}:{size_str}".replace(" ", "_")
            lines.append(f"{chrom}\t{start0}\t{end0}\t{name}\t.\t+")

    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# ASCII off-target map
# --------------------------------------------------------------------------- #

def ascii_offtarget_map(pair: dict, width: int = 56) -> str:
    """Return an ASCII diagram of predicted products for one primer pair.

    Per subject/chromosome, draws the two primer arrows around a span whose
    length is scaled (relative to that subject's largest product) into at most
    *width* columns. Orientation determines arrowheads:

        F/R  ->  ``F >----< R``   (convergent, a real amplicon)
        F/F  ->  ``F >----> F``   (both forward -> mispriming)
        R/R  ->  ``R <----< R``   (both reverse -> mispriming)

    The intended (on_target) product is marked ``(intended)``; others are
    labelled with their orientation as off-targets.
    """
    products = list(_get(pair, "products", []) or [])
    name = str(_get(pair, "name", "?"))

    if not products:
        return f"{name}: (no predicted products)\n"

    # Group products by subject, preserving first-seen order.
    by_subject: Dict[str, List[dict]] = {}
    for prod in products:
        subj = str(_get(prod, "subject", "?"))
        by_subject.setdefault(subj, []).append(prod)

    lines: List[str] = [f"Off-target map for {name}"]

    for subject, prods in by_subject.items():
        # Scale spans relative to the largest product on THIS subject.
        sizes = [p.get("size") for p in prods if isinstance(p.get("size"), (int, float))]
        max_size = max(sizes) if sizes else 1
        max_size = max_size or 1  # guard against zero

        for prod in prods:
            orientation = str(_get(prod, "orientation", "F/R")).upper()
            size = prod.get("size")
            on_target = bool(prod.get("on_target", False))

            # Determine per-product tag on the subject header line.
            if on_target:
                tag = "(intended)"
            else:
                tag = f"(off-target, {orientation})" if orientation else "(off-target)"
            lines.append(f"{subject}  {tag}")

            # Scale the span width; keep a small minimum so arrows are legible.
            if isinstance(size, (int, float)) and max_size:
                span = int(round((float(size) / float(max_size)) * width))
            else:
                span = width
            span = max(span, 4)

            # Choose arrowheads from orientation halves.
            parts = orientation.split("/")
            left_code = parts[0] if parts and parts[0] else "F"
            right_code = parts[1] if len(parts) > 1 and parts[1] else "R"

            left_head = ">" if left_code == "F" else "<"
            right_head = ">" if right_code == "F" else "<"

            dashes = "-" * span
            size_label = f"{int(size)} bp" if isinstance(size, (int, float)) else ""
            lines.append(
                f"  {left_code} {left_head}{dashes}{right_head} {right_code}"
                f"    {size_label}".rstrip()
            )

    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------- #
# Self-contained HTML report
# --------------------------------------------------------------------------- #

# Risk -> CSS colour mapping used for the summary table badges.
_RISK_COLORS = {
    "low": "#1a7f37",     # green
    "medium": "#bf8700",  # orange
    "high": "#cf222e",    # red
}

_HTML_CSS = """\
:root { color-scheme: light; }
* { box-sizing: border-box; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
    Helvetica, Arial, sans-serif;
  color: #1f2328; background: #ffffff; margin: 0; padding: 24px 32px;
  line-height: 1.45; font-size: 14px;
}
h1 { font-size: 22px; margin: 0 0 4px; }
h2 { font-size: 17px; margin: 28px 0 8px; border-bottom: 2px solid #d0d7de;
     padding-bottom: 4px; }
h3 { font-size: 15px; margin: 18px 0 6px; }
.meta { color: #57606a; font-size: 13px; margin: 2px 0; }
table { border-collapse: collapse; width: 100%; margin: 8px 0 4px;
        font-size: 13px; }
th, td { border: 1px solid #d0d7de; padding: 5px 8px; text-align: left;
         vertical-align: top; }
th { background: #f6f8fa; font-weight: 600; }
tr:nth-child(even) td { background: #fafbfc; }
code, .seq { font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo,
             Consolas, monospace; font-size: 12.5px; word-break: break-all; }
.badge { display: inline-block; padding: 1px 9px; border-radius: 10px;
         color: #fff; font-size: 12px; font-weight: 600; }
.pair { border: 1px solid #d0d7de; border-radius: 8px; padding: 12px 16px;
        margin: 12px 0; background: #fff; }
.pair h3 { margin-top: 0; }
.kv { display: flex; flex-wrap: wrap; gap: 4px 24px; margin: 6px 0; }
.kv div { min-width: 160px; }
.kv .label { color: #57606a; font-size: 12px; }
.flags span { display: inline-block; margin-right: 12px; font-size: 12.5px; }
.off { color: #cf222e; } .on { color: #1a7f37; }
footer { margin-top: 32px; border-top: 1px solid #d0d7de; padding-top: 12px;
         color: #57606a; font-size: 12px; }
footer pre { background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px;
             padding: 10px; overflow-x: auto; font-size: 12px; }
@media print {
  body { padding: 0; font-size: 12px; }
  .pair { break-inside: avoid; }
  h2 { break-after: avoid; }
}
"""


def _esc(val: Any) -> str:
    """HTML-escape any value (stringify first)."""
    return html.escape("" if val is None else str(val))


def _risk_badge(risk: Any) -> str:
    """Return an HTML badge span coloured by risk level."""
    r = str(risk or "").lower()
    color = _RISK_COLORS.get(r, "#57606a")
    label = _esc(risk or "n/a")
    return f'<span class="badge" style="background:{color}">{label}</span>'


def _summary_table(pairs: List[dict]) -> str:
    """Build the top-level summary table of all pairs."""
    cols = ["Rank", "Name", "Risk", "Product (bp)", "Tm F/R", "GC F/R",
            "Off-tgt", "F/F", "R/R", "F/R off", "3' mm min"]
    head = "".join(f"<th>{_esc(c)}</th>" for c in cols)
    rows_html: List[str] = []
    for i, pair in enumerate(pairs, start=1):
        tm = f"{_num(pair.get('tm_f'))} / {_num(pair.get('tm_r'))}"
        gc = f"{_num(pair.get('gc_f'))} / {_num(pair.get('gc_r'))}"
        cells = [
            str(i),
            _esc(_get(pair, "name")),
            _risk_badge(pair.get("risk")),
            _esc(_num(pair.get("product_size"))),
            _esc(tm),
            _esc(gc),
            _esc(_num(pair.get("n_off_target"))),
            _esc(_num(pair.get("n_ff"))),
            _esc(_num(pair.get("n_rr"))),
            _esc(_num(pair.get("n_fr_offtarget"))),
            _esc(_num(pair.get("tp5_mismatch_min"))),
        ]
        rows_html.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
    return (f'<table><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(rows_html)}</tbody></table>')


def _products_table(pair: dict) -> str:
    """Build the predicted-products table for one pair's detail section."""
    products = list(_get(pair, "products", []) or [])
    if not products:
        return '<p class="meta">No predicted products.</p>'
    cols = ["Subject", "Start", "End", "Size (bp)", "Orientation", "Target"]
    head = "".join(f"<th>{_esc(c)}</th>" for c in cols)
    rows_html: List[str] = []
    for prod in products:
        on_target = bool(prod.get("on_target", False))
        target_cls = "on" if on_target else "off"
        target_txt = "intended" if on_target else "off-target"
        cells = [
            _esc(_get(prod, "subject")),
            _esc(_num(prod.get("start")) or prod.get("start")),
            _esc(_num(prod.get("end")) or prod.get("end")),
            _esc(_num(prod.get("size"))),
            _esc(_get(prod, "orientation")),
            f'<span class="{target_cls}">{target_txt}</span>',
        ]
        rows_html.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
    return (f'<table><thead><tr>{head}</tr></thead>'
            f'<tbody>{"".join(rows_html)}</tbody></table>')


def _pair_detail(pair: dict, rank: int) -> str:
    """Build one per-pair detail section."""
    name = _esc(_get(pair, "name"))
    left = _get(pair, "left_pos", [])
    right = _get(pair, "right_pos", [])
    left_s = "-".join(str(x) for x in left) if isinstance(left, (list, tuple)) else _esc(left)
    right_s = "-".join(str(x) for x in right) if isinstance(right, (list, tuple)) else _esc(right)

    conserved = _get(pair, "conserved_refs", [])
    if isinstance(conserved, (list, tuple)):
        conserved_s = ", ".join(_esc(c) for c in conserved)
    else:
        conserved_s = _esc(conserved)

    # 3' off-target mismatch minimum is a specificity hint.
    tp5 = _num(pair.get("tp5_mismatch_min"))
    snp = _bool_str(pair.get("snp_in_primer"))
    caps = _esc(_get(pair, "caps_enzyme")) or "n/a"
    gel = _bool_str(pair.get("gel_distinguishable"))

    return f"""\
<div class="pair">
  <h3>#{rank} &middot; {name} &nbsp; {_risk_badge(pair.get("risk"))}</h3>
  <div class="kv">
    <div><span class="label">Forward (5'-&gt;3')</span><br>
         <span class="seq">{_esc(_get(pair, "forward"))}</span></div>
    <div><span class="label">Reverse (5'-&gt;3')</span><br>
         <span class="seq">{_esc(_get(pair, "reverse"))}</span></div>
  </div>
  <div class="kv">
    <div><span class="label">Product size</span><br>{_esc(_num(pair.get("product_size")))} bp</div>
    <div><span class="label">Tm F / R</span><br>{_esc(_num(pair.get("tm_f")))} / {_esc(_num(pair.get("tm_r")))} &deg;C</div>
    <div><span class="label">GC F / R</span><br>{_esc(_num(pair.get("gc_f")))} / {_esc(_num(pair.get("gc_r")))} %</div>
    <div><span class="label">Forward pos</span><br>{_esc(left_s)}</div>
    <div><span class="label">Reverse pos</span><br>{_esc(right_s)}</div>
  </div>
  <h4 style="margin:12px 0 4px;font-size:13px;">Predicted products</h4>
  {_products_table(pair)}
  <div class="flags" style="margin-top:8px;">
    <span><b>CAPS enzyme:</b> {caps}</span>
    <span><b>Gel-distinguishable:</b> {gel or "n/a"}</span>
    <span><b>SNP in primer:</b> {snp or "n/a"}</span>
    <span><b>3' mismatch min (off-tgt):</b> {tp5 or "n/a"}</span>
    <span><b>Conserved in:</b> {conserved_s or "n/a"}</span>
  </div>
</div>"""


def html_report(context: dict) -> str:
    """Return a self-contained, printable HTML report string.

    ``context`` schema::

        {
          "title": str, "template": str, "databases": [...],
          "generated": str, "params": dict, "provenance": dict,
          "pairs": [pair dicts]
        }

    The output has inline CSS, no external resources, a light theme, and is
    designed to print cleanly to PDF.
    """
    title = _esc(_get(context, "title", "primerblast_oss report"))
    template = _esc(_get(context, "template", ""))
    generated = _esc(_get(context, "generated", ""))
    databases = _get(context, "databases", []) or []
    if isinstance(databases, (list, tuple)):
        db_str = ", ".join(_esc(d) for d in databases)
    else:
        db_str = _esc(databases)

    pairs = list(_get(context, "pairs", []) or [])

    # Provenance + params serialised for reproducibility (pretty JSON).
    provenance = _get(context, "provenance", {}) or {}
    params = _get(context, "params", {}) or {}
    prov_blob = {"params": params, "provenance": provenance}
    try:
        prov_json = json.dumps(prov_blob, indent=2, ensure_ascii=False, default=str)
    except (TypeError, ValueError):
        prov_json = str(prov_blob)

    details = "\n".join(_pair_detail(p, i) for i, p in enumerate(pairs, start=1))

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
{_HTML_CSS}</style>
</head>
<body>
<header>
  <h1>{title}</h1>
  <p class="meta"><b>Template:</b> {template or "n/a"}</p>
  <p class="meta"><b>Databases:</b> {db_str or "n/a"}</p>
  <p class="meta"><b>Generated:</b> {generated or "n/a"}</p>
</header>

<h2>Summary ({len(pairs)} primer pair{"s" if len(pairs) != 1 else ""})</h2>
{_summary_table(pairs)}

<h2>Primer pair details</h2>
{details}

<footer>
  <p><b>Provenance &amp; parameters</b> (for reproducibility)</p>
  <pre>{_esc(prov_json)}</pre>
</footer>
</body>
</html>
"""


# --------------------------------------------------------------------------- #
# Self-test / demo
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    # Two sample pair dicts exercising the full schema.
    clean_pair: Dict[str, Any] = {
        "name": "P1",
        "forward": "ACGTGGTCAACGGATTTGCAC",
        "reverse": "TTGCACCAGTTGAGCTTCGAT",
        "product_size": 368,
        "tm_f": 60.0, "tm_r": 60.1, "gc_f": 55.0, "gc_r": 45.0,
        "left_pos": [85338374, 85338394],
        "right_pos": [85338721, 85338741],
        "risk": "low",
        "n_off_target": 0, "n_ff": 0, "n_rr": 0, "n_fr_offtarget": 0,
        "tp5_mismatch_min": 0,
        "snp_in_primer": False,
        "conserved_refs": ["cameor_v2", "ZW6"],
        "caps_enzyme": "EcoRI",
        "gel_distinguishable": True,
        "products": [
            {"subject": "chr1", "start": 85338374, "end": 85338741,
             "size": 368, "orientation": "F/R", "on_target": True},
        ],
    }

    messy_pair: Dict[str, Any] = {
        "name": "P2",
        "forward": "GGATCCAATGCGTTAGCCTGA",
        "reverse": "CTGCAGTTACCGGATTACGGT",
        "product_size": 402,
        "tm_f": 59.4, "tm_r": 61.2, "gc_f": 52.4, "gc_r": 52.4,
        "left_pos": [12045, 12065],
        "right_pos": [12427, 12447],
        "risk": "medium",
        "n_off_target": 1, "n_ff": 1, "n_rr": 0, "n_fr_offtarget": 1,
        "tp5_mismatch_min": 2,
        "snp_in_primer": True,
        "conserved_refs": ["cameor_v2"],
        "caps_enzyme": "",
        "gel_distinguishable": False,
        "products": [
            {"subject": "chr2", "start": 12046, "end": 12447,
             "size": 402, "orientation": "F/R", "on_target": True},
            {"subject": "chr5", "start": 285451473, "end": 285454732,
             "size": 3260, "orientation": "F/F", "on_target": False},
            {"subject": "chr7", "start": 9001000, "end": 9001900,
             "size": 900, "orientation": "F/R", "on_target": False},
        ],
    }

    pairs = [clean_pair, messy_pair]

    print("=" * 70)
    print("CSV")
    print("=" * 70)
    print(pairs_to_csv(pairs))

    print("=" * 70)
    print("ORDER TABLE")
    print("=" * 70)
    print(order_table(pairs))

    print("=" * 70)
    print("BED")
    print("=" * 70)
    print(products_to_bed(pairs))

    print("=" * 70)
    print("ASCII OFF-TARGET MAPS")
    print("=" * 70)
    for p in pairs:
        print(ascii_offtarget_map(p))

    print("=" * 70)
    print("GENERIC CSV helper")
    print("=" * 70)
    print(to_csv_generic(
        [{"a": 1, "b": 2}, {"a": 3}],
        ["a", "b", "c"],
    ))

    # Build and write the HTML report.
    context = {
        "title": "PrimerBLAST-OSS report: PsCIK2 locus",
        "template": "cameor_v2:chr1:85,330,000-85,345,000",
        "databases": ["cameor_v2", "ZW6", "nr/nt"],
        "generated": "2026-07-08 12:00 JST",
        "params": {
            "primer3": {"opt_tm": 60, "opt_size": 20, "product_size_range": "300-500"},
            "specificity": {"min_3p_mismatch": 2, "max_product": 4000},
        },
        "provenance": {
            "tool": "primerblast_oss",
            "version": "0.1.0",
            "primer3_version": "2.6.1",
            "blast_version": "2.15.0+",
            "command": "primerblast-oss run --template PsCIK2.fa",
        },
        "pairs": pairs,
    }
    html_out = html_report(context)
    out_path = "/tmp/pbo_sample_report.html"
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html_out)
    print("=" * 70)
    print(f"HTML written to {out_path} ({len(html_out)} chars)")
    print("=" * 70)
