#!/usr/bin/env python3
"""Fixed release-acceptance checks for the scientific PCR engine.

This is an implementation/release gate, not wet-lab calibration. It combines
the real-tool continuous synthetic benchmark with a small set of pinned
scientific invariants that are especially costly to regress.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from benchmarks.continuous_benchmark import (  # noqa: E402
    FWD,
    REV,
    TARGET_AMPLICON,
    require_tool,
    run_benchmark,
    write_fai,
)
from primerblast_oss import __version__  # noqa: E402
from primerblast_oss.caps import ENZYME_METADATA, caps_scan  # noqa: E402
from primerblast_oss.gel import genotype_pattern_discrimination  # noqa: E402
from primerblast_oss.genome import revcomp  # noqa: E402
from primerblast_oss.provenance import (  # noqa: E402
    PROVENANCE_SCHEMA_VERSION,
    model_versions,
)
from primerblast_oss.specificity import (  # noqa: E402
    SEARCH_COMPLETE,
    SpecParams,
    _priming_sites_from_output,
    _realign_hit_to_site,
    _specificity_verdict,
    pair_specificity,
)


class MemoryGenome:
    def __init__(self, sequence):
        self.sequence = sequence

    def length(self, _subject):
        return len(self.sequence)

    def fetch(self, _subject, start, end, strand="+"):
        sequence = self.sequence[start - 1:end]
        return revcomp(sequence) if strand == "-" else sequence


def _hsp(primer, start, end):
    return list(map(str, [
        "query", "synthetic", 100, len(primer), 0, 0, 1, len(primer),
        start, end, "1e-8", 40, "plus", primer, primer, len(primer),
    ]))


def _build_clean_db(folder: Path):
    fasta = folder / "clean.fa"
    fasta.write_text(">chr_target\n%s\n" % TARGET_AMPLICON, encoding="utf-8")
    write_fai(fasta)
    db = folder / "clean_db"
    subprocess.run([
        require_tool("makeblastdb"), "-in", str(fasta), "-dbtype", "nucl",
        "-out", str(db), "-parse_seqids",
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return fasta, db


def _fixed_checks() -> dict:
    checks = {}

    with tempfile.TemporaryDirectory(prefix="primerblast-release-") as directory:
        root = Path(directory)
        _fasta, db = _build_clean_db(root)
        sp = SpecParams(min_product=40, max_product=220, max_target_seqs=5000)
        clean = pair_specificity(
            FWD, REV, str(db), designed_size=len(TARGET_AMPLICON), sp=sp)
        checks["clean_pair_is_specific"] = (
            clean["specific"] is True
            and clean["n_on_target"] == 1
            and clean["n_off_target"] == 0
        )

    primer = "GCTAGCTACGATCGTACGTA"
    insertion_target = primer[:10] + "A" + primer[10:]
    genome = MemoryGenome("N" * 100 + insertion_target + "N" * 100)
    site, resolved = _realign_hit_to_site(
        _hsp(primer, 101, 121), "F", primer, SpecParams(), genome)
    checks["full_primer_indel_realign"] = bool(
        resolved and site is not None
        and site.total_mismatch == 1
        and site.alignment_source == "full_length_genome"
        and "-" in (site.aligned_query or "")
    )

    malformed = "primer\tchr1"
    _sites, stats = _priming_sites_from_output(
        "ACGTACGTACGTACGTACGT", "F", malformed, SpecParams())
    checks["malformed_search_is_indeterminate"] = (
        stats.completeness != SEARCH_COMPLETE
        and stats.malformed_rows == 1
        and _specificity_verdict(True, stats.completeness) is None
    )

    ref = "A" * 100 + "GAATTC" + "A" * 100
    alt = "A" * 100 + "GACTTC" + "A" * 100
    caps = caps_scan(
        ref, alt, enzymes={"EcoRI": ENZYME_METADATA["EcoRI"]},
        gel_min_gap=25)
    checks["caps_digest_geometry"] = bool(
        len(caps) == 1
        and caps[0].distinguishable
        and sum(caps[0].allele_a_fragments) == len(ref)
        and sum(caps[0].allele_b_fragments) == len(alt)
        and caps[0].allele_a_fragments != caps[0].allele_b_fragments
    )

    separated = genotype_pattern_discrimination(
        [430], [293, 137], background_sizes=[800])
    ambiguous = genotype_pattern_discrimination(
        [430], [293, 137], background_sizes=[430])
    checks["gel_separated_offtarget_is_scoreable"] = (
        separated["distinguishable"] is True)
    checks["comigrating_offtarget_is_ambiguous"] = (
        ambiguous["distinguishable"] is False
        and ambiguous["comparisons"]["AB_vs_BB"]["distinguishable"] is False
    )
    return checks


def run_acceptance(max_seconds: float) -> dict:
    require_tool("blastn")
    require_tool("makeblastdb")
    require_tool("primer3_core")

    started = time.perf_counter()
    continuous = run_benchmark(max_seconds)
    checks = _fixed_checks()
    checks["continuous_real_tool_suite"] = continuous["passed"]
    checks["duplicate_offtarget_detected"] = bool(
        continuous["checks"].get("detects_duplicate_amplicon")
        and continuous["checks"].get("duplicate_pair_not_specific")
        and continuous["checks"].get("duplicate_pair_comigrating"))
    checks["same_primer_product_detected"] = bool(
        continuous["checks"].get("detects_ff_amplicon"))
    checks["thermo_path_runs"] = continuous["checks"].get("thermo_path_runs")

    required = {
        key: value for key, value in checks.items()
        if key != "thermo_path_runs"
    }
    thermo_ok = checks["thermo_path_runs"] is not False
    passed = all(value is True for value in required.values()) and thermo_ok
    return {
        "benchmark": "release_acceptance",
        "passed": passed,
        "package_version": __version__,
        "provenance_schema_version": PROVENANCE_SCHEMA_VERSION,
        "model_versions": model_versions(),
        "source_revision": (
            os.environ.get("GITHUB_SHA")
            or os.environ.get("PRIMERBLAST_GIT_REVISION")
        ),
        "elapsed_seconds": round(time.perf_counter() - started, 4),
        "checks": checks,
        "continuous": continuous,
        "scope_note": (
            "Implementation acceptance on synthetic fixtures; not wet-lab "
            "validation or calibrated experimental success probability."
        ),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-seconds", type=float, default=60.0)
    parser.add_argument("--json-out")
    args = parser.parse_args(argv)

    result = run_acceptance(args.max_seconds)
    payload = json.dumps(result, indent=2, sort_keys=True)
    print(payload)
    if args.json_out:
        Path(args.json_out).write_text(payload + "\n", encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
