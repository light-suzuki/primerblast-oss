"""Primer design via the primer3_core binary (Boulder-IO).

We call the primer3_core executable directly rather than primer3-py so the
package has no compiled dependency; the binary ships with the standard
`primer3` package on Debian/Ubuntu.
"""
from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass
class DesignParams:
    """Primer3 design parameters (Primer-BLAST-like defaults)."""

    product_size_ranges: Sequence[Tuple[int, int]] = ((70, 1000),)
    opt_size: int = 20
    min_size: int = 18
    max_size: int = 25
    opt_tm: float = 60.0
    min_tm: float = 57.0
    max_tm: float = 63.0
    max_tm_diff: float = 3.0
    min_gc: float = 20.0
    max_gc: float = 80.0
    max_poly_x: int = 5
    max_self_any_th: float = 45.0
    max_self_end_th: float = 35.0
    max_hairpin_th: float = 24.0
    salt_monovalent: float = 50.0
    salt_divalent: float = 1.5
    dntp_conc: float = 0.6
    dna_conc: float = 50.0
    num_return: int = 10
    target: Optional[Tuple[int, int]] = None
    included_region: Optional[Tuple[int, int]] = None
    excluded_regions: Sequence[Tuple[int, int]] = ()


@dataclass
class PrimerPair:
    index: int
    template_id: str
    forward: str
    reverse: str
    left_start: int
    left_len: int
    right_start: int
    right_len: int
    product_size: int
    tm_f: float
    tm_r: float
    gc_f: float
    gc_r: float
    self_any_th: float = 0.0
    self_end_th: float = 0.0
    penalty: float = 0.0
    specificity: Dict = field(default_factory=dict)

    @property
    def left_3p(self) -> int:
        return self.left_start + self.left_len - 1

    @property
    def right_3p(self) -> int:
        return self.right_start - self.right_len + 1


def clean_sequence(seq: str) -> str:
    """Normalize a template without changing its biological coordinates.

    Whitespace is FASTA formatting and is removed. Canonical bases are retained;
    every ambiguity or unsupported non-whitespace symbol is represented by ``N``
    instead of being deleted. This prevents artificial sequence junctions and
    keeps Primer3 coordinates aligned with the supplied template.
    """
    compact = re.sub(r"\s+", "", seq.upper())
    return "".join(base if base in "ACGT" else "N" for base in compact)


def _ambiguous_regions(seq: str) -> List[Tuple[int, int]]:
    """Return 0-based ``(start, length)`` runs containing non-ACGT bases."""
    regions: List[Tuple[int, int]] = []
    start: Optional[int] = None
    for i, base in enumerate(seq):
        if base not in "ACGT":
            if start is None:
                start = i
        elif start is not None:
            regions.append((start, i - start))
            start = None
    if start is not None:
        regions.append((start, len(seq) - start))
    return regions


def _merge_regions(regions: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping/adjacent 0-based regions for deterministic Boulder-IO."""
    spans = sorted((start, start + length) for start, length in regions
                   if start >= 0 and length > 0)
    merged: List[List[int]] = []
    for start, end in spans:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end - start) for start, end in merged]


def read_fasta(path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    buf: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(buf)))
                name = line[1:].split()[0] if len(line) > 1 else "seq"
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        records.append((name, "".join(buf)))
    return records


def _detect_primer3(explicit: Optional[str]) -> str:
    for cand in (explicit, "primer3_core"):
        if cand and shutil.which(cand):
            return shutil.which(cand)  # type: ignore[return-value]
    raise RuntimeError(
        "primer3_core not found. Install the 'primer3' package or pass primer3_bin."
    )


def _build_boulder(template_id: str, seq: str, p: DesignParams) -> str:
    ranges = " ".join(f"{lo}-{hi}" for lo, hi in p.product_size_ranges)
    lines = [
        f"SEQUENCE_ID={template_id}",
        f"SEQUENCE_TEMPLATE={seq}",
        "PRIMER_TASK=generic",
        "PRIMER_PICK_LEFT_PRIMER=1",
        "PRIMER_PICK_INTERNAL_OLIGO=0",
        "PRIMER_PICK_RIGHT_PRIMER=1",
        f"PRIMER_NUM_RETURN={p.num_return}",
        f"PRIMER_OPT_SIZE={p.opt_size}",
        f"PRIMER_MIN_SIZE={p.min_size}",
        f"PRIMER_MAX_SIZE={p.max_size}",
        f"PRIMER_OPT_TM={p.opt_tm}",
        f"PRIMER_MIN_TM={p.min_tm}",
        f"PRIMER_MAX_TM={p.max_tm}",
        f"PRIMER_PAIR_MAX_DIFF_TM={p.max_tm_diff}",
        f"PRIMER_MIN_GC={p.min_gc}",
        f"PRIMER_MAX_GC={p.max_gc}",
        f"PRIMER_MAX_POLY_X={p.max_poly_x}",
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1",
        f"PRIMER_MAX_SELF_ANY_TH={p.max_self_any_th}",
        f"PRIMER_MAX_SELF_END_TH={p.max_self_end_th}",
        f"PRIMER_MAX_HAIRPIN_TH={p.max_hairpin_th}",
        f"PRIMER_SALT_MONOVALENT={p.salt_monovalent}",
        f"PRIMER_SALT_DIVALENT={p.salt_divalent}",
        f"PRIMER_DNTP_CONC={p.dntp_conc}",
        f"PRIMER_DNA_CONC={p.dna_conc}",
        f"PRIMER_PRODUCT_SIZE_RANGE={ranges}",
        "PRIMER_EXPLAIN_FLAG=1",
    ]
    if p.included_region:
        lines.append(f"SEQUENCE_INCLUDED_REGION={p.included_region[0]},{p.included_region[1]}")
    if p.target:
        lines.append(f"SEQUENCE_TARGET={p.target[0]},{p.target[1]}")

    # Explicitly prevent primers from landing on assembly gaps/masked sequence.
    # These regions are added without changing the coordinate system.
    excluded = _merge_regions([*p.excluded_regions, *_ambiguous_regions(seq)])
    if excluded:
        excl = " ".join(f"{start},{length}" for start, length in excluded)
        lines.append(f"SEQUENCE_EXCLUDED_REGION={excl}")
    lines.append("=")
    return "\n".join(lines) + "\n"


def _parse_boulder(out: str, template_id: str) -> Tuple[List[PrimerPair], str]:
    kv: Dict[str, str] = {}
    for line in out.splitlines():
        if "=" in line and not line.startswith("="):
            k, v = line.split("=", 1)
            kv[k] = v
    n = int(kv.get("PRIMER_PAIR_NUM_RETURNED", "0"))
    pairs: List[PrimerPair] = []
    for i in range(n):
        try:
            lpos = kv[f"PRIMER_LEFT_{i}"].split(",")
            rpos = kv[f"PRIMER_RIGHT_{i}"].split(",")
            pair = PrimerPair(
                index=i,
                template_id=template_id,
                forward=kv[f"PRIMER_LEFT_{i}_SEQUENCE"],
                reverse=kv[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                left_start=int(lpos[0]),
                left_len=int(lpos[1]),
                right_start=int(rpos[0]),
                right_len=int(rpos[1]),
                product_size=int(kv[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]),
                tm_f=float(kv[f"PRIMER_LEFT_{i}_TM"]),
                tm_r=float(kv[f"PRIMER_RIGHT_{i}_TM"]),
                gc_f=float(kv[f"PRIMER_LEFT_{i}_GC_PERCENT"]),
                gc_r=float(kv[f"PRIMER_RIGHT_{i}_GC_PERCENT"]),
                self_any_th=float(kv.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH", 0.0)),
                self_end_th=float(kv.get(f"PRIMER_PAIR_{i}_COMPL_END_TH", 0.0)),
                penalty=float(kv.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0)),
            )
            if set(pair.forward) <= set("ACGT") and set(pair.reverse) <= set("ACGT"):
                pairs.append(pair)
        except KeyError:
            continue
    explain = kv.get("PRIMER_PAIR_EXPLAIN", "")
    return pairs, explain


def design_primers(
    template_id: str,
    sequence: str,
    params: Optional[DesignParams] = None,
    primer3_bin: Optional[str] = None,
) -> Tuple[List[PrimerPair], str]:
    """Design primer pairs for a template. Returns ``(pairs, Primer3 explain)``."""
    params = params or DesignParams()
    seq = clean_sequence(sequence)
    if len(seq) < params.min_size * 2:
        raise ValueError(f"Template '{template_id}' too short ({len(seq)} bp) for design.")
    exe = _detect_primer3(primer3_bin)
    boulder = _build_boulder(template_id, seq, params)
    proc = subprocess.run(
        [exe], input=boulder.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    if proc.returncode != 0:
        raise RuntimeError(f"primer3_core failed: {proc.stderr.decode(errors='ignore')}")
    return _parse_boulder(proc.stdout.decode(errors="ignore"), template_id)
