"""Unit tests for atomic makeblastdb builds (no external tools needed)."""
import os
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from primerblast_oss import tools  # noqa: E402


def _write_fasta(tmp: Path) -> str:
    fasta = tmp / "genome.fa"
    fasta.write_text(">chr1\nAAAA\n")
    return str(fasta)


def _fake_makeblastdb(monkeypatch, suffixes, returncode=0, stderr=b""):
    """Intercept subprocess.run: write index files next to the -out prefix."""
    def fake_run(cmd, **kwargs):
        temp_out = cmd[cmd.index("-out") + 1]
        for suffix in suffixes:
            Path(temp_out + suffix).write_bytes(b"x" * 4)
        return SimpleNamespace(
            returncode=returncode, stdout=b"", stderr=stderr)

    monkeypatch.setattr(tools.subprocess, "run", fake_run)


def _prefix_files(tmp: Path, prefix: str):
    return sorted(path.suffix for path in tmp.glob(prefix + ".*"))


def test_failed_build_leaves_no_partial_db_at_final_prefix(tmp_path, monkeypatch):
    fasta = _write_fasta(tmp_path)
    out = str(tmp_path / "db")
    _fake_makeblastdb(monkeypatch, [".nhr", ".nin"], returncode=1,
                      stderr=b"makeblastdb: error")
    try:
        tools.make_blastdb(fasta, out=out, makeblastdb_bin="makeblastdb")
    except RuntimeError as error:
        assert "makeblastdb: error" in str(error)
    else:
        raise AssertionError("non-zero exit must raise")
    assert _prefix_files(tmp_path, "db") == []
    assert _prefix_files(tmp_path, "db.tmp.") == []


def test_successful_build_moves_files_and_returns_same_prefix(tmp_path, monkeypatch):
    fasta = _write_fasta(tmp_path)
    out = str(tmp_path / "db")
    _fake_makeblastdb(monkeypatch, [".nhr", ".nin", ".nsq"], returncode=0)
    assert tools.make_blastdb(fasta, out=out, makeblastdb_bin="makeblastdb") == out
    assert sorted(_prefix_files(tmp_path, "db")) == [".nhr", ".nin", ".nsq"]
    assert _prefix_files(tmp_path, "db.tmp.") == []


def test_success_exit_without_full_index_is_treated_as_failure(tmp_path, monkeypatch):
    fasta = _write_fasta(tmp_path)
    out = str(tmp_path / "db")
    _fake_makeblastdb(monkeypatch, [".nhr", ".nsq"], returncode=0)
    try:
        tools.make_blastdb(fasta, out=out, makeblastdb_bin="makeblastdb")
    except RuntimeError as error:
        assert "complete index" in str(error)
        assert ".nin" in str(error)
    else:
        raise AssertionError("missing index files must raise")
    assert _prefix_files(tmp_path, "db") == []
    assert _prefix_files(tmp_path, "db.tmp.") == []


def test_existing_db_is_untouched_by_failed_build(tmp_path, monkeypatch):
    fasta = _write_fasta(tmp_path)
    out = tmp_path / "db"
    (out).with_suffix(".nhr").write_bytes(b"original")
    (out).with_suffix(".nin").write_bytes(b"original")
    _fake_makeblastdb(monkeypatch, [".nhr"], returncode=1, stderr=b"boom")
    try:
        tools.make_blastdb(fasta, out=str(out), makeblastdb_bin="makeblastdb")
    except RuntimeError:
        pass
    else:
        raise AssertionError("non-zero exit must raise")
    assert (out).with_suffix(".nhr").read_bytes() == b"original"
    assert (out).with_suffix(".nin").read_bytes() == b"original"


def test_cleanup_failure_reports_leftover_temp_prefix(tmp_path, monkeypatch):
    fasta = _write_fasta(tmp_path)
    out = str(tmp_path / "db")
    _fake_makeblastdb(monkeypatch, [".nhr"], returncode=1, stderr=b"boom")
    original_unlink = Path.unlink

    def failing_unlink(self):
        raise OSError("simulated cleanup failure")

    monkeypatch.setattr(Path, "unlink", failing_unlink)
    try:
        tools.make_blastdb(fasta, out=out, makeblastdb_bin="makeblastdb")
    except RuntimeError as error:
        assert "leftover temp prefix" in str(error)
        assert "db.tmp." in str(error)
    else:
        raise AssertionError("non-zero exit must raise")
    monkeypatch.setattr(Path, "unlink", original_unlink)


def test_missing_makeblastdb_is_reported():
    try:
        tools.make_blastdb("x.fa", out="db", makeblastdb_bin="")
    except RuntimeError as error:
        assert "makeblastdb not found" in str(error)
    else:
        raise AssertionError("missing binary must raise")
