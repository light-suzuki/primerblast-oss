#!/usr/bin/env python3
"""Compare legacy per-primer BLAST launches with the batched screening path."""
from __future__ import annotations

import argparse
import json
import shutil
import time
import tracemalloc

from primerblast_oss.specificity import (
    SpecParams,
    _detect_blastn,
    priming_sites_with_stats,
    screen_primers_with_stats,
)


def measure(function):
    tracemalloc.start()
    started = time.perf_counter()
    function()
    elapsed = time.perf_counter() - started
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return elapsed, peak / (1024 * 1024)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--reference-label", default="local-reference")
    parser.add_argument("--primer", required=True)
    parser.add_argument("--sizes", default="10,100,500")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--json-out")
    args = parser.parse_args()

    blastn = _detect_blastn(shutil.which("blastn"))
    sp = SpecParams(num_threads=args.threads)
    rows = []
    for count in [int(x) for x in args.sizes.split(",") if x.strip()]:
        primers = {"P%05d" % i: args.primer for i in range(count)}

        legacy_s, legacy_peak = measure(lambda: [
            priming_sites_with_stats(seq, name, args.db, sp, blastn)
            for name, seq in primers.items()
        ])
        batch_s, batch_peak = measure(lambda: screen_primers_with_stats(
            primers, args.db, sp, blastn))

        rows.append({
            "n_primers": count,
            "legacy_seconds": round(legacy_s, 4),
            "batch_seconds": round(batch_s, 4),
            "speedup": round(legacy_s / batch_s, 3) if batch_s else None,
            "legacy_python_peak_mib": round(legacy_peak, 3),
            "batch_python_peak_mib": round(batch_peak, 3),
        })

    payload = {"database": args.reference_label, "results": rows}
    print(json.dumps(payload, indent=2))
    if args.json_out:
        with open(args.json_out, "w") as handle:
            json.dump(payload, handle, indent=2)


if __name__ == "__main__":
    main()
