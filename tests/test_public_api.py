"""Contract tests for the stable public library API (no external tools).

These tests pin the versioned surface described in docs/PUBLIC_API.md:
JSON-safe versioned results, structured errors, cancellation, and the
compatibility aliases SnapyGene relies on.
"""
import json

import primerblast_oss
from primerblast_oss import api
from primerblast_oss.design import PrimerPair
from primerblast_oss.errors import (
    BlastError,
    CancelledError,
    InvalidDatabaseError,
    MalformedInputError,
    Primer3Error,
    PrimerblastError,
    SearchIncompleteError,
    ToolMissingError,
)
from primerblast_oss.pipeline import PipelineResult
from primerblast_oss.specificity import (
    Amplicon,
    PrimerHitStats,
    PrimingSite,
    SpecParams,
)


def test_version_and_capabilities_contract():
    assert primerblast_oss.__version__
    caps = primerblast_oss.capabilities()
    assert caps["api_version"] == primerblast_oss.API_VERSION
    assert caps["package_version"] == primerblast_oss.__version__
    assert set(caps["capabilities"]) == {
        "design", "pair_specificity", "pool_in_silico_pcr", "tiling",
        "multiplex", "qpcr_probe_design", "blast_database_creation",
        "thermodynamic_filtering",
    }
    assert caps["capabilities"]["design"] is True
    assert caps["capabilities"]["qpcr_probe_design"] is False


def test_discover_tools_shape():
    found = primerblast_oss.discover_tools()
    assert set(found["tools"]) == {"primer3_core", "blastn", "makeblastdb"}
    for info in found["tools"].values():
        assert set(info) == {"path", "version", "available"}
        assert info["available"] == (info["path"] is not None)
    assert found["complete"] == (found["missing"] == [])
    for name in found["missing"]:
        assert "Install" in found["messages"][name]


def test_error_classes_are_runtime_compatible():
    for cls in (PrimerblastError, ToolMissingError, InvalidDatabaseError,
                Primer3Error, BlastError, MalformedInputError,
                SearchIncompleteError, CancelledError):
        assert issubclass(cls, RuntimeError)
        assert issubclass(cls, PrimerblastError)


def test_legacy_adapter_imports_still_work():
    from primerblast_oss.design import DesignParams
    from primerblast_oss.pipeline import run_pipeline
    from primerblast_oss.specificity import (
        in_silico_pcr,
        pair_specificity,
        spec_params_for_profile,
    )
    assert DesignParams and run_pipeline
    assert spec_params_for_profile("local-strict")
    assert callable(pair_specificity) and callable(in_silico_pcr)


def test_public_symbols_resolve():
    for name in primerblast_oss.__all__:
        assert hasattr(primerblast_oss, name), name


def test_json_safe_roundtrips_through_json():
    amp = Amplicon("chr1", 100, 467, 368, "F", "R", 0, 0, on_target=True)
    amp.__dict__["nearest_gap"] = 42
    pair = PrimerPair(
        index=0, template_id="t", forward="A" * 20, reverse="T" * 20,
        left_start=0, left_len=20, right_start=99, right_len=20,
        product_size=100, tm_f=60.0, tm_r=60.0, gc_f=50.0, gc_r=50.0,
        specificity={"per_db": [{"on_target": [amp], "off_target": []}]},
    )
    payload = api.json_safe({"pairs": [pair]})
    reloaded = json.loads(json.dumps(payload))
    assert reloaded == payload
    target = reloaded["pairs"][0]["specificity"]["per_db"][0]["on_target"][0]
    assert target == {"subject": "chr1", "start": 100, "end": 467, "size": 368,
                      "fwd_primer": "F", "rev_primer": "R",
                      "fwd_mismatch": 0, "rev_mismatch": 0, "on_target": True,
                      "fwd_tp5": 0, "rev_tp5": 0, "fwd_tm": None,
                      "rev_tm": None, "fwd_end3_dg": None, "rev_end3_dg": None,
                      "nearest_gap": 42}


def test_design_and_screen_propagates_missing_tool(monkeypatch):
    import primerblast_oss.design as design

    def missing(*_args, **_kwargs):
        raise ToolMissingError(
            "primer3_core not found. Install the 'primer3' package or pass primer3_bin.")

    monkeypatch.setattr(design, "_detect_primer3", missing)
    try:
        api.design_and_screen("t", "ACGT" * 25, ["db"])
    except PrimerblastError as error:
        assert isinstance(error, RuntimeError)
        assert "primer3_core" in str(error)
    else:
        raise AssertionError("missing primer3_core must raise")


def test_pair_specificity_missing_blastn_raises(monkeypatch):
    import primerblast_oss.specificity as specificity

    def missing(*_args, **_kwargs):
        raise ToolMissingError("blastn not found. Install BLAST+ or pass blastn_bin.")

    monkeypatch.setattr(specificity, "_detect_blastn", missing)
    try:
        api.pair_specificity_result("A" * 20, "T" * 20, "db")
    except ToolMissingError as error:
        assert "blastn" in str(error)
    else:
        raise AssertionError("missing blastn must raise")


def _mocked_screen(monkeypatch):
    """Patch BLAST away so pair specificity runs on fixed priming sites."""
    import primerblast_oss.specificity as specificity

    sites = [
        PrimingSite("F", "chr1", "+", 119, 0, 0, plen=20),
        PrimingSite("R", "chr1", "-", 280, 0, 0, plen=20),
    ]
    stats = {
        "F": PrimerHitStats("F", 1, 1, 1, False, False),
        "R": PrimerHitStats("R", 1, 1, 1, False, False),
    }
    monkeypatch.setattr(specificity, "_detect_blastn", lambda value: "blastn")
    monkeypatch.setattr(
        specificity, "screen_primers_with_stats", lambda *args: (sites, stats))
    return specificity


def test_pair_specificity_result_is_json_safe(monkeypatch):
    _mocked_screen(monkeypatch)
    result = api.pair_specificity_result(
        "A" * 20, "T" * 20, "db", designed_size=200,
        sp=SpecParams(min_product=40, max_product=1000),
    )
    assert result["api_version"] == api.API_VERSION
    assert isinstance(result["on_target"][0], dict)
    assert json.loads(json.dumps(result)) == result


def test_pool_in_silico_pcr_is_json_safe(monkeypatch):
    _mocked_screen(monkeypatch)
    result = api.pool_in_silico_pcr(
        {"F": "A" * 20, "R": "T" * 20}, "db",
        sp=SpecParams(min_product=40, max_product=1000),
    )
    assert result["api_version"] == api.API_VERSION
    assert isinstance(result["products"][0], dict)
    assert json.loads(json.dumps(result)) == result


def test_cancel_check_aborts_long_running_calls():
    for call in (
        lambda: api.design_and_screen("t", "ACGT" * 25, ["db"], cancel_check=lambda: True),
        lambda: api.pair_specificity_result("A" * 20, "T" * 20, "db", cancel_check=lambda: True),
        lambda: api.pool_in_silico_pcr({"F": "A" * 20}, "db", cancel_check=lambda: True),
    ):
        try:
            call()
        except CancelledError:
            continue
        raise AssertionError("cancel_check must raise CancelledError")


def test_multiple_databases_reach_the_pipeline(monkeypatch):
    captured = {}

    def fake_run_pipeline(template_id, sequence, databases, **kwargs):
        captured["databases"] = list(databases)
        return PipelineResult(
            template_id=template_id, template_len=100, pairs=[],
            primer3_explain="", databases=list(databases))

    monkeypatch.setattr(api, "run_pipeline", fake_run_pipeline)
    out = api.design_and_screen("t", "ACGT" * 25, ["dbA", "dbB"])
    assert captured["databases"] == ["dbA", "dbB"]
    assert out["databases"] == ["dbA", "dbB"]


def test_strict_search_raises_on_incomplete_evidence(monkeypatch):
    pair = PrimerPair(
        index=0, template_id="t", forward="A" * 20, reverse="T" * 20,
        left_start=0, left_len=20, right_start=99, right_len=20,
        product_size=100, tm_f=60.0, tm_r=60.0, gc_f=50.0, gc_r=50.0,
        specificity={"search_completeness": "possibly_truncated"},
    )
    result = PipelineResult(
        template_id="t", template_len=100, pairs=[pair],
        primer3_explain="", databases=["db"],
    )
    monkeypatch.setattr(api, "run_pipeline", lambda *args, **kwargs: result)
    try:
        api.design_and_screen("t", "ACGT" * 25, ["db"], strict_search=True)
    except SearchIncompleteError as error:
        assert "possibly_truncated" in str(error)
    else:
        raise AssertionError("strict search must raise on incomplete evidence")
    out = api.design_and_screen("t", "ACGT" * 25, ["db"])
    assert out["api_version"] == api.API_VERSION
    assert (out["pairs"][0]["specificity"]["search_completeness"]
            == "possibly_truncated")


def test_create_database_contract(monkeypatch):
    monkeypatch.setattr(api, "make_blastdb", lambda *_a, **_k: "/tmp/db")
    result = api.create_database("g.fa", out="/tmp/db")
    assert result["api_version"] == api.API_VERSION
    assert result["db"] == "/tmp/db"
    assert result["parse_seqids"] is True
