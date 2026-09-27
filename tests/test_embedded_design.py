import subprocess
import pytest
from primerblast_oss import design, workbench_design
from primerblast_oss.errors import Primer3Error


def test_workbench_adapter_uses_shared_runner_without_changing_settings(monkeypatch):
    seen = []
    def runner(text, **kwargs):
        seen.append((text, kwargs))
        return "PRIMER_PAIR_NUM_RETURNED=0\n=\n"
    monkeypatch.setattr(workbench_design, "run_boulder", runner)
    monkeypatch.setenv("PRIMER3_TIMEOUT_SEC", "12")
    assert workbench_design.design_primers("acgt " * 100, target_start_1based=20,
        target_length=30, primer_max_size=27, primer_salt_monovalent=60) == []
    text, kwargs = seen[0]
    assert "SEQUENCE_TARGET=19,30" in text and "PRIMER_MAX_SIZE=27" in text
    assert "PRIMER_SALT_MONOVALENT=60" in text and kwargs["timeout"] == 12


def test_primer3_success_exit_with_error_is_not_an_empty_design(monkeypatch):
    monkeypatch.setattr(design, "_detect_primer3", lambda _: "synthetic-primer3")
    monkeypatch.setattr(design.subprocess, "run", lambda *a, **k:
        subprocess.CompletedProcess(a, 0, b"PRIMER_ERROR=Invalid target\n=\n", b""))
    with pytest.raises(Primer3Error, match="Invalid target"):
        design.run_boulder("synthetic")
