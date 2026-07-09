#!/usr/bin/env python3
"""Small continuous benchmark for CI and release checks.

This benchmark is deliberately self-contained: it builds a tiny synthetic FASTA
and BLAST database, then exercises the parts that tend to regress in a
Primer-BLAST-like workflow:

* primer hit pairing into F/R and F/F amplicons
* conservative duplicate/off-target classification
* Primer3 design + BLAST specificity
* optional primer3-py thermodynamic and multiplex checks

It is not a biological validation benchmark. The large, source-backed
PrimerServer2/NCBI comparisons belong in benchmarks/RESULTS.md.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from primerblast_oss.design import DesignParams, design_primers  # noqa: E402
from primerblast_oss.genome import Genome, revcomp               # noqa: E402
from primerblast_oss.pipeline import run_pipeline                # noqa: E402
from primerblast_oss.specificity import (                        # noqa: E402
    SpecParams, in_silico_pcr, pair_specificity,
)


FWD = "ACGTTGCAAGTCCGATCGTA"
REV = "TGACCGTATGCTAGCTTACG"
SPACER = "GATTACA" * 9
TARGET_AMPLICON = FWD + SPACER + revcomp(REV)


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SystemExit(f"required tool not found on PATH: {name}")
    return path


def write_fai(fasta: Path) -> None:
    offset = 0
    rows = []
    with fasta.open() as fh:
        name = None
        seq_start = None
        seq_len = 0
        line_bases = None
        line_width = None
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    rows.append((name, seq_len, seq_start, line_bases, line_width))
                name = line[1:].split()[0]
                seq_start = offset + len(raw)
                seq_len = 0
                line_bases = None
                line_width = None
            else:
                bases = len(line)
                if bases and line_bases is None:
                    line_bases = bases
                    line_width = len(raw)
                seq_len += bases
            offset += len(raw)
        if name is not None:
            rows.append((name, seq_len, seq_start, line_bases or seq_len, line_width or seq_len + 1))
    with (fasta.with_suffix(fasta.suffix + ".fai")).open("w") as out:
        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")


def build_fixture(tmp: Path) -> tuple[Path, Path]:
    fasta = tmp / "synthetic.fa"
    duplicate = ("N" * 50) + TARGET_AMPLICON + ("N" * 50)
    ff_product = FWD + ("C" * 80) + revcomp(FWD)
    design_template = ("ACGTGCAATGCTAGCTAGGCTAATCGGATCGTACGATCGATGCTAGCATCGATGCA" * 8)
    with fasta.open("w") as fh:
        fh.write(f">chr_target\n{TARGET_AMPLICON}\n")
        fh.write(f">chr_duplicate\n{duplicate}\n")
        fh.write(f">chr_ff\n{ff_product}\n")
        fh.write(f">chr_design\n{design_template}\n")
    write_fai(fasta)
    db = tmp / "synthetic_db"
    subprocess.run([
        require_tool("makeblastdb"), "-in", str(fasta), "-dbtype", "nucl",
        "-out", str(db), "-parse_seqids",
    ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return fasta, db


def timed(name: str, fn, timings: dict):
    t0 = time.perf_counter()
    result = fn()
    timings[name] = round(time.perf_counter() - t0, 4)
    return result


def run_benchmark(max_seconds: float) -> dict:
    require_tool("blastn")
    require_tool("makeblastdb")
    require_tool("primer3_core")

    timings = {}
    checks = {}
    with tempfile.TemporaryDirectory(prefix="primerblast-oss-bench-") as d:
        tmp = Path(d)
        fasta, db = build_fixture(tmp)
        genome = Genome(str(fasta))
        sp = SpecParams(min_product=40, max_product=220, max_target_seqs=5000)

        check_res = timed(
            "in_silico_pcr",
            lambda: in_silico_pcr({"F": FWD, "R": REV}, str(db), sp=sp),
            timings,
        )
        products = {(p.subject, p.size, p.orientation) for p in check_res["products"]}
        checks["detects_target_amplicon"] = ("chr_target", 103, "F/R") in products
        checks["detects_duplicate_amplicon"] = ("chr_duplicate", 103, "F/R") in products
        checks["detects_ff_amplicon"] = ("chr_ff", 120, "F/F") in products

        pair_res = timed(
            "pair_specificity",
            lambda: pair_specificity(FWD, REV, str(db), designed_size=103, sp=sp),
            timings,
        )
        checks["duplicate_pair_not_specific"] = pair_res["specific"] is False
        checks["duplicate_pair_comigrating"] = pair_res["n_comigrating"] >= 1

        thermo_res = timed(
            "pair_specificity_thermo",
            lambda: pair_specificity(FWD, REV, str(db), designed_size=103, sp=sp, genome=genome),
            timings,
        )
        checks["thermo_path_runs"] = "thermo_evaluated" in thermo_res

        design_seq = TARGET_AMPLICON + ("ACGT" * 40)
        pairs = timed(
            "primer3_design",
            lambda: design_primers("bench_template", design_seq,
                                   DesignParams(product_size_ranges=[(80, 180)], num_return=3))[0],
            timings,
        )
        checks["primer3_returns_pairs"] = len(pairs) > 0

        pipeline_res = timed(
            "pipeline",
            lambda: run_pipeline("bench_template", design_seq, [str(db)],
                                 design_params=DesignParams(product_size_ranges=[(80, 180)], num_return=2),
                                 spec_params=sp),
            timings,
        )
        checks["pipeline_returns_pairs"] = len(pipeline_res.pairs) > 0

        try:
            from primerblast_oss import dimers
            if dimers.available():
                mux = timed(
                    "multiplex",
                    lambda: dimers.analyze_multiplex([("F", FWD), ("R", REV)]),
                    timings,
                )
                checks["multiplex_runs"] = mux["n_pairs_checked"] == 1
            else:
                checks["multiplex_runs"] = None
        except Exception as exc:  # noqa: BLE001
            checks["multiplex_runs"] = False
            checks["multiplex_error"] = str(exc)

    total = round(sum(timings.values()), 4)
    passed = all(v is not False for v in checks.values()) and total <= max_seconds
    return {
        "benchmark": "continuous_synthetic",
        "passed": passed,
        "max_seconds": max_seconds,
        "total_timed_seconds": total,
        "timings": timings,
        "checks": checks,
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-seconds", type=float, default=30.0)
    ap.add_argument("--json-out")
    args = ap.parse_args(argv)

    result = run_benchmark(args.max_seconds)
    text = json.dumps(result, indent=2, sort_keys=True)
    if args.json_out:
        Path(args.json_out).write_text(text + "\n")
    print(text)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
