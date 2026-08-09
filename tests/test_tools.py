"""Unit tests for .fai-indexed FASTA access (no external tools needed)."""
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from primerblast_oss.tools import (  # noqa: E402
    FaidxContigError, FaidxCoordError, FaidxError, faidx_fetch,
)


_TEMP_DIRS = []  # keep TemporaryDirectory objects alive for the test session


def _fixture():
    """Build a 80 bp two-line FASTA (40 bp/line) with its .fai index.

    Returns ``(fasta_path, fai_path)`` inside a temporary directory.
    """
    tmp = tempfile.TemporaryDirectory()
    _TEMP_DIRS.append(tmp)
    fasta = Path(tmp.name) / "genome.fa"
    # write in binary so line endings stay "\n" on every platform (the .fai
    # linewidth below assumes single-byte newlines)
    fasta.write_bytes((">chr1\n" + "A" * 40 + "\n" + "C" * 40 + "\n").encode())
    fai = fasta.with_suffix(".fa.fai")
    fai.write_text("chr1\t80\t6\t40\t41\n")
    return str(fasta), str(fai)


def test_fetch_full_and_boundary_coordinates():
    fasta, _ = _fixture()
    assert faidx_fetch(fasta, "chr1", 1, 80) == "A" * 40 + "C" * 40
    assert faidx_fetch(fasta, "chr1", 1, 1) == "A"
    assert faidx_fetch(fasta, "chr1", 80, 80) == "C"
    assert faidx_fetch(fasta, "chr1", 40, 41) == "AC"


def test_fetch_across_line_boundary():
    fasta, _ = _fixture()
    assert faidx_fetch(fasta, "chr1", 38, 43) == "AAACCC"


def test_fetch_unknown_contig_raises_specific_error():
    fasta, _ = _fixture()
    try:
        faidx_fetch(fasta, "chrX", 1, 10)
    except FaidxContigError as error:
        assert "chrX" in str(error)
        assert error is not None  # it is a FaidxError, not a raw KeyError
    else:
        raise AssertionError("unknown contig must raise FaidxContigError")


def test_fetch_missing_fasta_names_path():
    try:
        faidx_fetch("/no/such/genome.fa", "chr1", 1, 10)
    except FaidxError as error:
        assert "/no/such/genome.fa" in str(error)
    else:
        raise AssertionError("missing FASTA must raise FaidxError")


def test_fetch_missing_fai_names_path():
    tmp = tempfile.TemporaryDirectory()
    _TEMP_DIRS.append(tmp)
    fasta = Path(tmp.name) / "genome.fa"
    fasta.write_text(">chr1\nAAAA\n")
    try:
        faidx_fetch(str(fasta), "chr1", 1, 4)
    except FaidxError as error:
        assert str(fasta) + ".fai" in str(error)
    else:
        raise AssertionError("missing .fai must raise FaidxError")


def test_fetch_start_greater_than_end_is_rejected():
    fasta, _ = _fixture()
    try:
        faidx_fetch(fasta, "chr1", 60, 50)
    except FaidxCoordError as error:
        assert "60" in str(error) and "50" in str(error)
    else:
        raise AssertionError("start > end must raise FaidxCoordError")


def test_fetch_start_below_one_is_rejected():
    fasta, _ = _fixture()
    try:
        faidx_fetch(fasta, "chr1", 0, 10)
    except FaidxCoordError as error:
        assert "start" in str(error)
    else:
        raise AssertionError("start < 1 must raise FaidxCoordError")


def test_fetch_end_beyond_contig_is_clamped():
    fasta, _ = _fixture()
    assert faidx_fetch(fasta, "chr1", 1, 100) == "A" * 40 + "C" * 40


def test_fetch_range_fully_beyond_contig_returns_empty():
    fasta, _ = _fixture()
    assert faidx_fetch(fasta, "chr1", 90, 100) == ""


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception:
            failed += 1
            print(f"FAIL {fn.__name__}")
            traceback.print_exc()
    print(f"\n{len(fns)-failed}/{len(fns)} passed")
    raise SystemExit(1 if failed else 0)
