"""Minimal, dependency-free BED interval parser/writer.

BED coordinates are 0-based, half-open: an interval [start, end) covers the
bases start .. end-1, so its length is simply end - start. This differs from
GFF3/VCF (1-based inclusive); keep the conventions straight when converting.

Only the standard library is used, and both plain and gzip-compressed
(`.bed.gz`) files are supported (detected by extension).
"""
from __future__ import annotations

import gzip
import io
from dataclasses import dataclass
from typing import List, Sequence


def _open_text(path: str) -> io.TextIOBase:
    """Open a possibly gzip-compressed file as UTF-8 text (by extension)."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


@dataclass
class BedInterval:
    """A BED interval. `start` is 0-based; `end` is exclusive (half-open)."""

    chrom: str
    start: int
    end: int
    name: str = ""
    score: str = "."
    strand: str = "."

    @property
    def length(self) -> int:
        """Interval length in bp (end - start)."""
        return self.end - self.start


def parse_bed(path: str) -> List[BedInterval]:
    """Parse a BED file (plain or `.gz`) into a list of :class:`BedInterval`.

    Rows are tab-split and may carry 3-6 columns (extra columns are ignored);
    `track`, `browser`, and comment (`#`) lines are skipped.
    """
    intervals: List[BedInterval] = []
    with _open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if line.startswith(("#", "track", "browser")):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            try:
                start = int(cols[1])
                end = int(cols[2])
            except ValueError:
                continue
            intervals.append(
                BedInterval(
                    chrom=cols[0],
                    start=start,
                    end=end,
                    name=cols[3] if len(cols) > 3 else "",
                    score=cols[4] if len(cols) > 4 else ".",
                    strand=cols[5] if len(cols) > 5 else ".",
                )
            )
    return intervals


def to_bed_line(
    chrom: str,
    start0: int,
    end: int,
    name: str = "",
    score: str = ".",
    strand: str = ".",
) -> str:
    """Format one 6-column BED line (no trailing newline). `start0` is
    0-based, `end` is exclusive."""
    return "\t".join([chrom, str(start0), str(end), name, score, strand])


def write_bed(intervals: Sequence[BedInterval]) -> str:
    """Serialize intervals to a 6-column, tab-separated BED string
    (newline-terminated per row)."""
    lines = [
        to_bed_line(iv.chrom, iv.start, iv.end, iv.name, iv.score, iv.strand)
        for iv in intervals
    ]
    return "".join(line + "\n" for line in lines)


if __name__ == "__main__":
    import io as _io

    intervals = [
        BedInterval("chr1", 5031, 6598, "Psat.cameor.v2.1g00050", ".", "+"),
        BedInterval("chr1", 25902, 27160, "Psat.cameor.v2.1g00100", ".", "+"),
    ]
    text = write_bed(intervals)
    print("write_bed output:")
    print(text, end="")

    # Round-trip: parse the text back through the same tab-split logic.
    parsed: List[BedInterval] = []
    for line in _io.StringIO(text):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        cols = line.split("\t")
        parsed.append(
            BedInterval(
                cols[0], int(cols[1]), int(cols[2]),
                cols[3] if len(cols) > 3 else "",
                cols[4] if len(cols) > 4 else ".",
                cols[5] if len(cols) > 5 else ".",
            )
        )
    print("round-trip parsed:")
    for iv in parsed:
        print(f"  {iv.chrom} {iv.start}-{iv.end} {iv.name} "
              f"strand={iv.strand} length={iv.length}")
    assert parsed == intervals, "round-trip mismatch"
    print("round-trip OK")
