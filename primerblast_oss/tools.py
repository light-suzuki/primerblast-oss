"""Helpers for building BLAST databases from FASTA."""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple


class FaidxError(RuntimeError):
    """Base class for .fai-indexed FASTA access failures.

    Subclasses distinguish unknown contigs from invalid coordinates so
    callers (CLI, GUI, public API) can point at the exact input to fix.
    """


class FaidxContigError(FaidxError):
    """The requested contig is absent from the .fai index."""


class FaidxCoordError(FaidxError):
    """The requested coordinates are invalid (start < 1 or end < start)."""


def _load_fai(fasta: str) -> Dict[str, Tuple[int, int, int, int]]:
    fai = fasta + ".fai"
    if not Path(fai).exists():
        raise RuntimeError(f"FASTA index not found: {fai} (run `samtools faidx {fasta}`)")
    index: Dict[str, Tuple[int, int, int, int]] = {}
    with open(fai) as fh:
        for line in fh:
            name, length, offset, linebases, linewidth = line.split("\t")[:5]
            index[name] = (int(length), int(offset), int(linebases), int(linewidth))
    return index


def faidx_fetch(fasta: str, name: str, start: int, end: int) -> str:
    """Fetch bases [start, end] (1-based inclusive) from an .fai-indexed FASTA
    by seeking, without loading the whole (multi-hundred-Mbp) sequence.

    Inputs are validated with explicit errors instead of a raw KeyError: a
    missing FASTA or .fai names the offending path, an unknown contig raises
    ``FaidxContigError``, and invalid coordinates (``start < 1`` or
    ``end < start``) raise ``FaidxCoordError``. An end beyond the contig
    length is clamped to the contig end (mirroring the clamp of start to 1);
    a range clamped to empty returns "".
    """
    fasta_p = Path(fasta)
    if not fasta_p.exists():
        raise FaidxError(f"FASTA not found: {fasta}")
    fai = fasta + ".fai"
    if not Path(fai).exists():
        raise FaidxError(
            f"FASTA index not found: {fai} (run `samtools faidx {fasta}`)")
    index = _load_fai(fasta)
    if name not in index:
        raise FaidxContigError(f"contig '{name}' not found in {fai}")
    if start < 1:
        raise FaidxCoordError(
            f"invalid start {start} for contig '{name}' (start must be >= 1)")
    if end < start:
        raise FaidxCoordError(
            f"invalid coordinates {start}-{end} for contig '{name}' "
            "(end must be >= start)")
    length, offset, linebases, linewidth = index[name]
    start = max(1, start)
    end = min(length, end)
    want = end - start + 1
    if want <= 0:
        return ""
    start0 = start - 1
    byte_start = offset + (start0 // linebases) * linewidth + (start0 % linebases)
    # bytes to read = wanted bases + the newline bytes interleaved among them
    n_newlines = want // linebases + 2
    with open(fasta, "rb") as fh:
        fh.seek(byte_start)
        raw = fh.read(want + n_newlines * (linewidth - linebases) + 4)
    seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:want]
    return seq.decode().upper()


def make_blastdb(
    fasta: str,
    out: Optional[str] = None,
    title: Optional[str] = None,
    parse_seqids: bool = True,
    makeblastdb_bin: Optional[str] = None,
) -> str:
    """Build a nucleotide BLAST database. Returns the db path prefix.

    parse_seqids is on by default so downstream tools can extract subject
    regions by accession; the existing local pea DBs were built without it.
    """
    exe = makeblastdb_bin or shutil.which("makeblastdb")
    if not exe:
        raise RuntimeError("makeblastdb not found. Install BLAST+.")
    fasta_p = Path(fasta)
    out = out or str(fasta_p.with_suffix(""))
    title = title or fasta_p.stem
    cmd = [exe, "-in", fasta, "-dbtype", "nucl", "-out", out, "-title", title]
    if parse_seqids:
        cmd.append("-parse_seqids")
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(f"makeblastdb failed: {proc.stderr.decode(errors='ignore')}")
    return out
