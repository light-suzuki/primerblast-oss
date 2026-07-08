"""Minimal, dependency-free VCF (v4.x) record parser.

VCF positions are 1-based; a record's span is [POS, POS+len(REF)-1] on the
reference (this module keeps that 1-based inclusive convention). Only the
standard library is used, and both plain and gzip-compressed (`.vcf.gz`)
files are supported (detected by extension).
"""
from __future__ import annotations

import gzip
import io
from dataclasses import dataclass, field
from typing import List, Optional, Sequence


def _open_text(path: str) -> io.TextIOBase:
    """Open a possibly gzip-compressed file as UTF-8 text (by extension)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


@dataclass
class Variant:
    """A single VCF record. POS is 1-based."""

    chrom: str
    pos: int
    id: str
    ref: str
    alt: List[str] = field(default_factory=list)
    qual: Optional[float] = None
    info: str = ""

    @property
    def is_snp(self) -> bool:
        """True when REF and every ALT allele are a single base."""
        return len(self.ref) == 1 and all(len(a) == 1 for a in self.alt)

    @property
    def end(self) -> int:
        """Last reference position spanned, 1-based inclusive."""
        return self.pos + len(self.ref) - 1


def parse_vcf(path: str) -> List[Variant]:
    """Parse a VCF file (plain or `.gz`) into a list of :class:`Variant`.

    `##` meta lines and the `#CHROM` header line are skipped. Each data
    record's first 8 columns (CHROM POS ID REF ALT QUAL FILTER INFO) are
    read; ALT is split on ',' and '.' QUAL becomes None.
    """
    variants: List[Variant] = []
    with _open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue  # covers both '##...' meta and '#CHROM' header
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            try:
                pos = int(cols[1])
            except ValueError:
                continue
            qual_raw = cols[5]
            qual: Optional[float]
            if qual_raw in (".", ""):
                qual = None
            else:
                try:
                    qual = float(qual_raw)
                except ValueError:
                    qual = None
            alt = [a for a in cols[4].split(",") if a] if cols[4] != "." else []
            variants.append(
                Variant(
                    chrom=cols[0],
                    pos=pos,
                    id=cols[2],
                    ref=cols[3],
                    alt=alt,
                    qual=qual,
                    info=cols[7],
                )
            )
    return variants


def variants_in(
    variants: Sequence[Variant], chrom: str, start: int, end: int
) -> List[Variant]:
    """Return variants on `chrom` whose span [pos, end] overlaps the
    1-based inclusive window [start, end]."""
    hits: List[Variant] = []
    for v in variants:
        if v.chrom != chrom:
            continue
        if v.end < start or v.pos > end:
            continue
        hits.append(v)
    return hits


if __name__ == "__main__":
    import os
    import tempfile

    # The on-disk VCF is empty, so exercise the parser on a synthetic file.
    sample = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "##source=synthetic",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "chr1\t100\trs1\tA\tG\t60\tPASS\tDP=30",
            "chr1\t150\trs2\tC\tT,A\t45.5\tPASS\tDP=25",
            "chr1\t200\trs3\tGATC\tG\t99\tPASS\tDP=40",  # deletion (indel)
            "chr1\t500\trs4\tT\tC\t.\tPASS\tDP=10",
            "chr2\t100\trs5\tA\tT\t50\tPASS\tDP=20",
        ]
    )
    tmp = os.path.join(tempfile.gettempdir(), "primerblast_oss_demo.vcf")
    with open(tmp, "w") as fh:
        fh.write(sample + "\n")

    variants = parse_vcf(tmp)
    print(f"parsed {len(variants)} variants")
    for v in variants:
        print(f"  {v.chrom}:{v.pos} {v.ref}>{','.join(v.alt)} "
              f"snp={v.is_snp} end={v.end} qual={v.qual}")

    window = variants_in(variants, "chr1", 150, 205)
    print("variants_in chr1:150-205:")
    for v in window:
        print(f"  {v.chrom}:{v.pos} {v.id}")
    os.remove(tmp)
