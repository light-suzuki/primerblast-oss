"""Unit tests for the wet-lab validation panel manifest and analysis."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from benchmarks.validation_analysis import (  # noqa: E402
    compute_metrics,
    format_summary,
    load_manifest,
    main,
    reconcile,
    verify_immutability,
)


def _assay(assay_id, category, status, rank="A", risk_level="low",
           risk_score=90.0, expected_bands=(200,), expected_off=(),
           band_count=None, extra=False, exclusion_reason=None,
           caps_concordant=None):
    prediction = {
        "expected_band_sizes": list(expected_bands),
        "expected_off_target_sizes": list(expected_off),
        "rank": rank,
        "risk_score": risk_score,
        "risk_level": risk_level,
        "search_completeness": "complete",
    }
    if category == "caps_marker":
        prediction["caps"] = {"enzyme": "TaqI"}
    return {
        "assay_id": assay_id,
        "category": category,
        "primers": {"forward": "A" * 20, "reverse": "T" * 20},
        "prediction": prediction,
        "observed": {
            "status": status,
            "band_count": band_count if band_count is not None else (
                2 if extra else 1),
            "band_sizes": [],
            "exclusion_reason": exclusion_reason,
            "caps_concordant": caps_concordant,
        },
    }


def _manifest(assays):
    return {"schema_version": "1.0", "panel_version": "1.0", "assays": assays}


def test_offtarget_sensitivity_and_fpr():
    records = reconcile(_manifest([
        _assay("A1", "recent_paralog", "extra_bands", expected_off=(1400,)),
        _assay("A2", "repeat", "confirmed", expected_off=(200, 200)),
        _assay("A3", "unique_locus", "confirmed"),
    ]))
    metrics = compute_metrics(records)
    assert metrics["n_assessed"] == 3
    assert metrics["offtarget_sensitivity"] == 1.0
    assert metrics["offtarget_false_positive_rate"] == 0.5


def test_missed_extra_band_lowers_sensitivity():
    records = reconcile(_manifest([
        _assay("A1", "organellar", "extra_bands", expected_off=()),
        _assay("A2", "recent_paralog", "extra_bands", expected_off=(900,)),
    ]))
    metrics = compute_metrics(records)
    assert metrics["offtarget_sensitivity"] == 0.5


def test_amplification_failure_and_calibration():
    records = reconcile(_manifest([
        _assay("A1", "unique_locus", "confirmed", rank="A", risk_level="low"),
        _assay("A2", "primer_site_variant", "failed_amplification",
               rank="C", risk_level="medium"),
        _assay("A3", "unique_locus", "confirmed", rank="A", risk_level="low"),
    ]))
    metrics = compute_metrics(records)
    assert metrics["amplification_failure_rate"] == 1 / 3
    assert metrics["amplification_success_rate"] == 2 / 3
    assert metrics["calibration_by_rank"]["A"]["success_rate"] == 1.0
    assert metrics["calibration_by_rank"]["C"]["success_rate"] == 0.0
    assert metrics["calibration_by_risk"]["low"]["success_rate"] == 1.0
    assert metrics["calibration_by_risk"]["medium"]["success_rate"] == 0.0


def test_exclusions_are_listed_not_counted():
    records = reconcile(_manifest([
        _assay("A1", "unique_locus", "confirmed"),
        _assay("A2", "negative_control", "excluded",
               exclusion_reason="contamination"),
    ]))
    metrics = compute_metrics(records)
    assert metrics["n_excluded"] == 1
    assert metrics["n_assessed"] == 1
    assert metrics["exclusions"] == [{
        "assay_id": "A2", "reason": "contamination"}]
    assert metrics["amplification_success_rate"] == 1.0


def test_caps_concordance():
    records = reconcile(_manifest([
        _assay("A1", "caps_marker", "confirmed", caps_concordant=True),
        _assay("A2", "caps_marker", "confirmed", caps_concordant=False),
        _assay("A3", "caps_marker", "pending", caps_concordant=None),
    ]))
    metrics = compute_metrics(records)
    assert metrics["n_caps_assessed"] == 2
    assert metrics["caps_concordance"] == 0.5


def test_repetitive_stratification():
    records = reconcile(_manifest([
        _assay("A1", "repeat", "confirmed"),
        _assay("A2", "tandem_duplication", "failed_amplification"),
        _assay("A3", "unique_locus", "confirmed"),
    ]))
    metrics = compute_metrics(records)
    assert metrics["success_repetitive"] == 0.5
    assert metrics["success_non_repetitive"] == 1.0


def test_hash_mismatch_is_flagged():
    import hashlib

    assay = _assay("A1", "unique_locus", "pending")
    joined = "%s|%s" % (assay["primers"]["forward"], assay["primers"]["reverse"])
    assay["primers"]["sequence_sha256"] = hashlib.sha256(joined.encode()).hexdigest()
    manifest = _manifest([assay])
    assert verify_immutability(manifest) == []
    assay["primers"]["forward"] = "C" * 20
    issues = verify_immutability(manifest)
    assert len(issues) == 1 and "A1" in issues[0]


def test_load_manifest_rejects_unknown_schema():
    import tempfile

    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
        json.dump({"schema_version": "9.9", "assays": []}, fh)
        path = fh.name
    try:
        try:
            load_manifest(path)
        except ValueError as error:
            assert "schema_version" in str(error)
        else:
            raise AssertionError("unknown schema must be rejected")
    finally:
        os.unlink(path)


def test_main_runs_on_example_manifest(capsys):
    example = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "validation_panel", "manifest.example.json")
    assert main(["--manifest", example]) == 0
    output = capsys.readouterr().out
    assert "Validation panel summary" in output
    assert "assays: 8 total" in output
    assert "EXCLUDED: VP-008" in output


if __name__ == "__main__":
    import inspect
    import traceback

    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and not inspect.signature(v).parameters]
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
