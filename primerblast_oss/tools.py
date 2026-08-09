"""Helpers for building BLAST databases from FASTA."""
from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_INDEX_SUFFIXES = (".nhr", ".nin", ".nsq")


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


def _cleanup_temp_files(temp_out: Path) -> List[Path]:
    """Remove every temporary index file; return the ones that could not be."""
    leftovers: List[Path] = []
    for path in temp_out.parent.glob(temp_out.name + ".*"):
        try:
            path.unlink()
        except OSError:
            leftovers.append(path)
    return leftovers


def _raise_cleanup_failed(kind: str, details: str, temp_out: Path,
                          leftovers: List[Path]) -> None:
    message = "makeblastdb %s: %s" % (kind, details)
    if leftovers:
        message += (
            " (could not remove temporary files; leftover temp prefix: %s)"
            % " ".join(str(path) for path in leftovers))
    raise RuntimeError(message)


def make_blastdb(
    fasta: str,
    out: Optional[str] = None,
    title: Optional[str] = None,
    parse_seqids: bool = True,
    makeblastdb_bin: Optional[str] = None,
) -> str:
    """Build a nucleotide BLAST database. Returns the db path prefix.

    The database is built under a temporary prefix in the same directory as
    the final prefix and moved there only after the run exits 0 and the core
    index files (``.nhr``/``.nin``/``.nsq``) are present. A failed or
    interrupted build therefore never leaves a partial database at the final
    prefix, and an existing database at that prefix is only ever replaced by
    a fully built one. If cleanup of the temporary files fails, the leftover
    temp prefix is named in the error message.

    parse_seqids is on by default so downstream tools can extract subject
    regions by accession; the existing local pea DBs were built without it.
    """
    exe = makeblastdb_bin or shutil.which("makeblastdb")
    if not exe:
        raise RuntimeError("makeblastdb not found. Install BLAST+.")
    fasta_p = Path(fasta)
    out = out or str(fasta_p.with_suffix(""))
    title = title or fasta_p.stem
    out_p = Path(out)
    temp_out = out_p.parent / (
        out_p.name + ".tmp.%d.%d" % (os.getpid(), int(time.time() * 1e6)))
    cmd = [exe, "-in", fasta, "-dbtype", "nucl", "-out", str(temp_out),
           "-title", title]
    if parse_seqids:
        cmd.append("-parse_seqids")
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        leftovers = _cleanup_temp_files(temp_out)
        _raise_cleanup_failed(
            "failed", proc.stderr.decode(errors="ignore").strip(), temp_out,
            leftovers)

    missing = [
        suffix for suffix in _INDEX_SUFFIXES
        if not (Path(str(temp_out) + suffix).is_file()
                and Path(str(temp_out) + suffix).stat().st_size > 0)
    ]
    if missing:
        leftovers = _cleanup_temp_files(temp_out)
        _raise_cleanup_failed(
            "finished but did not produce a complete index (missing %s)"
            % ", ".join(missing), "output may be incomplete", temp_out,
            leftovers)

    moved = 0
    for path in sorted(temp_out.parent.glob(temp_out.name + ".*")):
        final = path.parent / (out_p.name + path.suffix)
        try:
            os.replace(str(path), str(final))
        except OSError as error:
            leftovers = _cleanup_temp_files(temp_out)
            _raise_cleanup_failed(
                "failed to move %s to %s: %s" % (path.name, final, error),
                "", temp_out, leftovers)
        moved += 1
    if moved == 0:
        leftovers = _cleanup_temp_files(temp_out)
        _raise_cleanup_failed(
            "produced no index files", "no files matched %s.*" % temp_out.name,
            temp_out, leftovers)
    return out
