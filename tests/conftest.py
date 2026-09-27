import pytest


@pytest.fixture(autouse=True)
def isolated_reference_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("SEQWB_REFERENCE_CACHE_DIR", str(tmp_path / "reference-cache"))
