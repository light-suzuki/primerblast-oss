"""Analyze a wet-lab validation panel against in-silico predictions.

Computes the panel metrics described in validation_panel/README.md:
off-target sensitivity / false-positive rate, amplification failure rate,
and calibration of rank / risk against observed outcomes, stratified by
locus category and primer-site variant class.

Usage:
    python benchmarks/validation_analysis.py --manifest panel.json \
        [--csv summary.csv]

The manifest schema is versioned (see validation_panel/manifest.example.json).
Pure standard library; no plotting dependency.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

ASSESSED_STATUSES = (
    "confirmed",
    "weak_amplification",
    "failed_amplification",
    "extra_bands",
    "unexpected_sizes",
)
PENDING_STATUSES = ("pending",)
EXCLUDED_STATUSES = ("excluded",)

REPETITIVE_CATEGORIES = {
    "recent_paralog",
    "tandem_duplication",
    "repeat",
    "organellar",
    "high_copy",
}
RANK_ORDER = ["A", "B", "C", "D", "I"]
RISK_ORDER = ["low", "medium", "high"]


def load_manifest(path: str) -> Dict:
    """Load and schema-check a validation panel manifest."""
    with open(path, encoding="utf-8") as fh:
        manifest = json.load(fh)
    if not isinstance(manifest, dict) or "assays" not in manifest:
        raise ValueError("%s: not a validation panel manifest (missing 'assays')" % path)
    schema = str(manifest.get("schema_version", "unknown"))
    if not schema.startswith("1."):
        raise ValueError(
            "%s: unsupported schema_version %s (this tool reads schema 1.x)"
            % (path, schema))
    return manifest


def verify_immutability(manifest: Dict) -> List[str]:
    """Return integrity issues: primer hash mismatches or changed inputs."""
    issues: List[str] = []
    for assay in manifest.get("assays", []):
        assay_id = assay.get("assay_id", "?")
        primers = assay.get("primers", {})
        expected = primers.get("sequence_sha256")
        if not expected:
            continue
        joined = "%s|%s" % (primers.get("forward", ""), primers.get("reverse", ""))
        actual = hashlib.sha256(joined.encode()).hexdigest()
        if actual != expected:
            issues.append(
                "%s: primer sequence_sha256 mismatch (inputs were modified)" % assay_id)
    return issues


def _record(assay: Dict) -> Dict:
    observed = assay.get("observed", {}) or {}
    prediction = assay.get("prediction", {}) or {}
    expected_bands = prediction.get("expected_band_sizes", []) or []
    expected_off = prediction.get("expected_off_target_sizes", []) or []
    expected_band_count = len(expected_bands) or (1 if expected_bands else 0)
    status = observed.get("status", "pending")
    band_count = observed.get("band_count")
    observed_extra = (
        status in ("extra_bands", "unexpected_sizes")
        or (band_count is not None and band_count > max(1, expected_band_count))
    )
    return {
        "assay_id": assay.get("assay_id"),
        "category": assay.get("category", "unknown"),
        "status": status,
        "rank": prediction.get("rank"),
        "risk_score": prediction.get("risk_score"),
        "risk_level": prediction.get("risk_level"),
        "expected_band_count": expected_band_count,
        "expected_off_target_count": len(expected_off),
        "observed_band_count": band_count,
        "observed_band_sizes": observed.get("band_sizes", []) or [],
        "observed_extra": observed_extra,
        "search_completeness": prediction.get("search_completeness"),
        "caps_expected": bool(prediction.get("caps")),
        "caps_concordant": observed.get("caps_concordant"),
        "exclusion_reason": observed.get("exclusion_reason"),
        "replicate_agreement": observed.get("replicate_agreement"),
    }


def reconcile(manifest: Dict) -> List[Dict]:
    """Flatten the manifest into per-assay analysis records."""
    return [_record(assay) for assay in manifest.get("assays", [])]


def _success_rate(subset: List[Dict]) -> Optional[float]:
    n = len(subset)
    if n == 0:
        return None
    confirmed = sum(1 for record in subset if record["status"] == "confirmed")
    return confirmed / n


def compute_metrics(records: List[Dict]) -> Dict:
    """Panel metrics with denominators; None means "no data for this cell"."""
    assessed = [r for r in records if r["status"] in ASSESSED_STATUSES]
    excluded = [r for r in records if r["status"] in EXCLUDED_STATUSES]
    pending = [r for r in records if r["status"] in PENDING_STATUSES]

    with_extra = [r for r in assessed if r["observed_extra"]]
    predicted_off = [r for r in assessed if r["expected_off_target_count"] > 0]
    sensitivity = (
        sum(1 for r in with_extra if r["expected_off_target_count"] > 0)
        / len(with_extra) if with_extra else None
    )
    fpr = (
        sum(1 for r in predicted_off if not r["observed_extra"])
        / len(predicted_off) if predicted_off else None
    )

    intended = [r for r in assessed if r["expected_band_count"] >= 1]
    failed = [r for r in intended if r["status"] == "failed_amplification"]

    calibration_rank = {
        rank: {"n": n, "success_rate": _success_rate(
            [r for r in intended if r["rank"] == rank])}
        for rank in RANK_ORDER
        if (n := sum(1 for r in intended if r["rank"] == rank)) > 0
    }
    calibration_risk = {
        level: {"n": n, "success_rate": _success_rate(
            [r for r in intended if r["risk_level"] == level])}
        for level in RISK_ORDER
        if (n := sum(1 for r in intended if r["risk_level"] == level)) > 0
    }

    repetitive = [r for r in assessed if r["category"] in REPETITIVE_CATEGORIES]
    non_repetitive = [r for r in assessed if r["category"] not in REPETITIVE_CATEGORIES]

    caps = [r for r in assessed if r["caps_expected"]
            and r["caps_concordant"] is not None]
    caps_concordance = (
        sum(1 for r in caps if r["caps_concordant"]) / len(caps) if caps else None
    )

    return {
        "n_assays": len(records),
        "n_assessed": len(assessed),
        "n_pending": len(pending),
        "n_excluded": len(excluded),
        "exclusions": [{
            "assay_id": r["assay_id"],
            "reason": r["exclusion_reason"] or "not stated",
        } for r in records if r["status"] in EXCLUDED_STATUSES],
        "integrity_issues": [],
        "offtarget_sensitivity": sensitivity,
        "offtarget_false_positive_rate": fpr,
        "amplification_failure_rate": (
            len(failed) / len(intended) if intended else None),
        "amplification_success_rate": _success_rate(intended),
        "calibration_by_rank": calibration_rank,
        "calibration_by_risk": calibration_risk,
        "success_repetitive": _success_rate(repetitive),
        "success_non_repetitive": _success_rate(non_repetitive),
        "caps_concordance": caps_concordance,
        "n_caps_assessed": len(caps),
    }


def format_summary(metrics: Dict, integrity_issues: List[str]) -> str:
    pct = lambda value: "n/a" if value is None else "%.0f%%" % (100.0 * value)
    lines = [
        "Validation panel summary",
        "  assays: %d total, %d assessed, %d pending, %d excluded"
        % (metrics["n_assays"], metrics["n_assessed"],
           metrics["n_pending"], metrics["n_excluded"]),
        "Off-target detection (experimentally observed extra bands)",
        "  sensitivity = %s" % pct(metrics["offtarget_sensitivity"]),
        "  false-positive rate of predicted off-targets = %s"
        % pct(metrics["offtarget_false_positive_rate"]),
        "Intended amplification",
        "  success rate = %s" % pct(metrics["amplification_success_rate"]),
        "  failure rate = %s" % pct(metrics["amplification_failure_rate"]),
        "  success by locus class: repetitive %s / non-repetitive %s"
        % (pct(metrics["success_repetitive"]),
           pct(metrics["success_non_repetitive"])),
        "Calibration",
    ]
    for rank in RANK_ORDER:
        if rank in metrics["calibration_by_rank"]:
            cell = metrics["calibration_by_rank"][rank]
            lines.append("  rank %s: n=%d success=%s"
                         % (rank, cell["n"], pct(cell["success_rate"])))
    for level in RISK_ORDER:
        if level in metrics["calibration_by_risk"]:
            cell = metrics["calibration_by_risk"][level]
            lines.append("  risk %s: n=%d success=%s"
                         % (level, cell["n"], pct(cell["success_rate"])))
    if metrics["n_caps_assessed"]:
        lines.append("CAPS predicted-vs-observed concordance: %s (n=%d)"
                     % (pct(metrics["caps_concordance"]),
                        metrics["n_caps_assessed"]))
    for issue in integrity_issues:
        lines.append("INTEGRITY: %s" % issue)
    for exclusion in metrics["exclusions"]:
        lines.append("EXCLUDED: %s (%s)"
                     % (exclusion["assay_id"], exclusion["reason"]))
    return "\n".join(lines)


def write_csv(records: List[Dict], path: str) -> None:
    columns = [
        "assay_id", "category", "status", "rank", "risk_level", "risk_score",
        "expected_bands", "expected_off_targets", "observed_bands",
        "observed_extra", "caps_concordant", "replicate_agreement",
        "search_completeness", "exclusion_reason",
    ]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Analyze a wet-lab validation panel against predictions.")
    parser.add_argument("--manifest", required=True,
                        help="validation panel manifest JSON (schema 1.x)")
    parser.add_argument("--csv", help="optional per-assay CSV output path")
    arguments = parser.parse_args(argv)

    manifest = load_manifest(arguments.manifest)
    integrity_issues = verify_immutability(manifest)
    records = reconcile(manifest)
    metrics = compute_metrics(records)
    metrics["integrity_issues"] = integrity_issues
    print(format_summary(metrics, integrity_issues))
    if arguments.csv:
        write_csv(records, arguments.csv)
        print("wrote %s" % arguments.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
