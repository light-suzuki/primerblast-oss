"""Random-access genome sequence via a `.fai` index (samtools faidx format).

Lets the pipeline pull a template region straight out of a multi-hundred-Mbp
chromosome by seeking, and convert between genomic coordinates and a local
template -- the anchor that makes every primer's strand/coordinates explicit.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple

_COMP = str.maketrans("ACGTNacgtnRYSWKMBDHVryswkmbdhv",
                      "TGCANtgcanYRSWMKVHDByrswmkvhdb")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass
class _FaiEntry:
    length: int
    offset: int
    linebases: int
    linewidth: int


class Genome:
    """A FASTA + .fai index. Coordinates are 1-based inclusive."""

    def __init__(self, fasta: str):
        self.fasta = fasta
        self.fai_path = fasta + ".fai"
        if not Path(self.fai_path).exists():
            raise RuntimeError(
                f"FASTA index not found: {self.fai_path} (run `samtools faidx {fasta}`)")
        self.index: Dict[str, _FaiEntry] = {}
        with open(self.fai_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                name, length, offset, lb, lw = parts[:5]
                self.index[name] = _FaiEntry(int(length), int(offset), int(lb), int(lw))

    def __contains__(self, name: str) -> bool:
        return name in self.index

    def length(self, name: str) -> int:
        return self.index[name].length

    def chroms(self):
        return list(self.index.keys())

    def fetch(self, name: str, start: int, end: int, strand: str = "+") -> str:
        """Fetch bases [start, end] (1-based inclusive). strand '-' returns the
        reverse complement, so the result always reads 5'->3' on that strand."""
        if name not in self.index:
            raise KeyError(f"sequence '{name}' not in {self.fai_path}")
        e = self.index[name]
        start = max(1, start)
        end = min(e.length, end)
        if end < start:
            return ""
        want = end - start + 1
        start0 = start - 1
        byte_start = e.offset + (start0 // e.linebases) * e.linewidth + (start0 % e.linebases)
        n_lines = want // e.linebases + 2
        with open(self.fasta, "rb") as fh:
            fh.seek(byte_start)
            raw = fh.read(want + n_lines * (e.linewidth - e.linebases) + 4)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:want].decode().upper()
        return revcomp(seq) if strand == "-" else seq

    def local_to_genomic(self, name: str, region_start: int, strand: str,
                         local_index0: int) -> int:
        """Map a 0-based index on a template extracted at `region_start`
        (1-based genomic start of the extraction) back to a genomic coordinate.
        For '+' the template runs 5'->3' with genome; for '-' it is revcomp."""
        if strand == "-":
            # region_start here is the genomic END (higher coord) of extraction
            return region_start - local_index0
        return region_start + local_index0
