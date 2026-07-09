#!/usr/bin/env python3
"""Benchmark primerblast-oss on a local pea genome.

Extracts a real genomic region as the PCR template, designs primers, and
screens their specificity against one or more cultivar databases -- then
reports timing and how many predicted products each pair yields.
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from primerblast_oss.design import DesignParams
from primerblast_oss.specificity import SpecParams
from primerblast_oss.pipeline import run_pipeline
from primerblast_oss import report as report_mod
from primerblast_oss.tools import faidx_fetch

# Path to a directory holding local BLAST databases + the matching FASTA(.fai).
# Override for your machine: `export PBO_DBDIR=/path/to/blastdb`.
DBDIR = os.environ.get("PBO_DBDIR", "/home/user/.codex/blast_databases")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", default=f"{DBDIR}/pisum_v2.fa")
    ap.add_argument("--chrom", default="chr1")
    ap.add_argument("--start", type=int, default=85337500)
    ap.add_argument("--end", type=int, default=85339200)
    ap.add_argument("--db", action="append",
                    default=None, help="repeatable; default = pisum_v2 self-search")
    ap.add_argument("--product-min", type=int, default=300)
    ap.add_argument("--product-max", type=int, default=800)
    ap.add_argument("--num-return", type=int, default=8)
    ap.add_argument("--json-out")
    args = ap.parse_args()

    dbs = args.db or [f"{DBDIR}/pisum_v2"]

    t0 = time.time()
    template = faidx_fetch(args.fasta, args.chrom, args.start, args.end)
    t_extract = time.time() - t0
    tid = f"{args.chrom}:{args.start}-{args.end}"
    print(f"# template {tid}  {len(template)} bp  (extracted in {t_extract*1000:.0f} ms)")
    print(f"# databases: {', '.join(dbs)}")

    dp = DesignParams(
        product_size_ranges=[(args.product_min, args.product_max)],
        num_return=args.num_return,
    )
    sp = SpecParams()

    t1 = time.time()
    result = run_pipeline(tid, template, dbs, design_params=dp, spec_params=sp)
    elapsed = time.time() - t1
    print(f"# design + specificity for {len(result.pairs)} pairs in {elapsed:.1f} s "
          f"({elapsed/max(1,len(result.pairs)):.2f} s/pair)\n")

    print(report_mod.to_text(result))

    if args.json_out:
        with open(args.json_out, "w") as fh:
            fh.write(report_mod.to_json(result))
        print(f"\n# wrote {args.json_out}")


if __name__ == "__main__":
    main()
