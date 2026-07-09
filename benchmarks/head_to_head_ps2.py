#!/usr/bin/env python3
"""Automated multi-locus head-to-head: primerblast-oss vs PrimerServer2.

For each of N sequence windows spread across a genome, this:

  1. designs a primer pair with Primer3 (via primerblast-oss),
  2. predicts that pair's amplicons with **primerblast-oss** in-silico PCR
     (BLAST-hit pairing + optional primer3-py thermodynamic gate), and
  3. predicts the same pair's amplicons with **PrimerServer2** `primertool check`,

then compares the two predictions per locus on amplicon **count**, **product
sizes**, and **genomic coordinates**, and prints a concordance summary.

Both tools run against the *same* local BLAST database, so this measures whether
primerblast-oss's specificity engine agrees with a mature, published local tool
(PrimerServer2) at scale — not just on a hand-picked locus.

Usage:
    python benchmarks/head_to_head_ps2.py \
        --genome /path/tair10.fa --db /path/tair10.fa \
        --primertool /path/to/primertool \
        --n-loci 30 --window 1500 --product-size 100-500

`--genome` must be `.fai`-indexed and its subject names must match the BLAST db.
PrimerServer2 needs `samtools` and BLAST+ on PATH.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from primerblast_oss.genome import Genome                                  # noqa: E402
from primerblast_oss.design import DesignParams, design_primers            # noqa: E402
from primerblast_oss.specificity import SpecParams, in_silico_pcr          # noqa: E402


def spaced_loci(genome: Genome, n: int, window: int, max_n_frac: float = 0.02):
    """Pick up to n low-N windows spread evenly across all chromosomes.

    Deterministic (no RNG): walks each chromosome at a fixed stride and keeps
    windows whose N-fraction is below max_n_frac, until n are collected.
    """
    chroms = list(genome.chroms())
    lengths = {c: genome.length(c) for c in chroms}
    total = sum(lengths.values())
    loci = []
    # budget windows per chromosome proportional to its length
    for c in chroms:
        share = max(1, round(n * lengths[c] / total))
        usable = lengths[c] - 2 * window
        if usable <= 0:
            continue
        stride = max(window, usable // (share + 1))
        pos = window
        while pos < lengths[c] - window and len(loci) < n:
            seq = genome.fetch(c, pos, pos + window, "+").upper()
            if seq and seq.count("N") / len(seq) <= max_n_frac:
                loci.append((f"{c}_{pos}", c, pos, seq))
            pos += stride
        if len(loci) >= n:
            break
    return loci[:n]


def design_top_pair(locus_id: str, seq: str, size_ranges, primer3_bin=None):
    dp = DesignParams(product_size_ranges=size_ranges, num_return=5)
    pairs, _ = design_primers(locus_id, seq, dp, primer3_bin)
    return pairs[0] if pairs else None


def oss_amplicons(fwd: str, rev: str, db: str, genome, thermo_params, gate,
                  size_min: int, size_max: int, profile_sp=None, blastn_bin=None):
    sp = profile_sp or SpecParams()
    sp.min_product, sp.max_product = size_min, size_max
    res = in_silico_pcr({"F": fwd, "R": rev}, db, sp=sp,
                        blastn_bin=blastn_bin, genome=genome,
                        thermo_params=thermo_params, thermo_gate=gate)
    amps = []
    for a in res["products"]:
        amps.append((int(a.size), f"{a.subject}:{a.start}-{a.end}"))
    return amps


def ps2_amplicons(site_id: str, fwd: str, rev: str, db_fasta: str, primertool: str,
                  tm_diff: float, cpu: int, size_min: int, size_max: int):
    """Run PrimerServer2 `check` for one pair and parse its predicted amplicons."""
    with tempfile.TemporaryDirectory() as td:
        q = os.path.join(td, "q.txt")
        out = os.path.join(td, "q.json")
        with open(q, "w") as fh:
            fh.write(f"{site_id} {fwd} {rev}\n")
        cmd = [primertool, "check", q, db_fasta, "-o", out, "-t",
               os.path.join(td, "q.tsv"), "-p", str(cpu), "--Tm-diff", str(int(tm_diff)),
               "--checking-size-min", str(size_min),
               "--checking-size-max", str(size_max),
               "--amplicon-num-max", "1000"]
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        d = json.load(open(out))
    site = d["primers"][site_id]
    dbkey = next((k for k in site if isinstance(site[k], dict)
                  and any(kk.startswith("PRIMER_PAIR_0_AMPLICON") for kk in site[k])), None)
    amps = []
    if dbkey:
        for amp in site[dbkey].get("PRIMER_PAIR_0_AMPLICONS", []):
            amps.append((int(amp["product_size"]), amp["region"]))
    return amps


def _key(amps):
    return sorted(amps)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help=".fai-indexed genome FASTA")
    ap.add_argument("--db", help="BLAST db prefix (default: --genome value)")
    ap.add_argument("--primertool", default=os.environ.get("PRIMERTOOL", "primertool"),
                    help="PrimerServer2 primertool executable")
    ap.add_argument("--n-loci", type=int, default=30)
    ap.add_argument("--window", type=int, default=1500)
    ap.add_argument("--product-size", default="100-500")
    ap.add_argument("--tm-diff", type=float, default=20.0,
                    help="PrimerServer2 --Tm-diff (matches our thermo floor)")
    ap.add_argument("--check-size-min", type=int, default=50,
                    help="off-target amplicon size window, shared by both tools")
    ap.add_argument("--check-size-max", type=int, default=2000)
    ap.add_argument("--specificity-profile", default="local-strict",
                    choices=["local-strict", "ncbi"],
                    help="primerblast-oss mismatch profile for the check")
    ap.add_argument("--cpu", type=int, default=4)
    ap.add_argument("--no-thermo", action="store_true",
                    help="disable primerblast-oss thermodynamic gate (mismatch model only)")
    ap.add_argument("--out", help="write full per-locus JSON here")
    ap.add_argument("--primer3-bin")
    ap.add_argument("--blastn-bin")
    a = ap.parse_args(argv)

    db = a.db or a.genome
    lo, hi = a.product_size.split("-")
    size_ranges = [(int(lo), int(hi))]

    from primerblast_oss.specificity import spec_params_for_profile
    profile_sp = spec_params_for_profile(a.specificity_profile)

    genome = Genome(a.genome)
    thermo_params, gate = None, True
    if not a.no_thermo:
        from primerblast_oss import thermo as _thermo
        if _thermo.available():
            thermo_params = _thermo.ThermoParams()
        else:
            print("warning: primer3-py not installed; running mismatch-model only",
                  file=sys.stderr)

    loci = spaced_loci(genome, a.n_loci, a.window)
    print(f"selected {len(loci)} low-N windows of {a.window} bp across "
          f"{len(list(genome.chroms()))} sequences\n", file=sys.stderr)

    rows = []
    n_designed = count_match = size_match = coord_match = 0
    for i, (lid, chrom, pos, seq) in enumerate(loci, 1):
        pair = design_top_pair(lid, seq, size_ranges, a.primer3_bin)
        if pair is None:
            rows.append({"locus": lid, "designed": False})
            print(f"[{i:2}/{len(loci)}] {lid:16} design: no pair", file=sys.stderr)
            continue
        n_designed += 1
        sid = f"L{i}"
        oss = oss_amplicons(pair.forward, pair.reverse, db, genome,
                            thermo_params, gate, a.check_size_min, a.check_size_max,
                            profile_sp=profile_sp, blastn_bin=a.blastn_bin)
        ps2 = ps2_amplicons(sid, pair.forward, pair.reverse, a.genome,
                            a.primertool, a.tm_diff, a.cpu,
                            a.check_size_min, a.check_size_max)
        cnt = len(oss) == len(ps2)
        sz = sorted(s for s, _ in oss) == sorted(s for s, _ in ps2)
        crd = _key(oss) == _key(ps2)
        count_match += cnt
        size_match += sz
        coord_match += crd
        rows.append({"locus": lid, "designed": True,
                     "forward": pair.forward, "reverse": pair.reverse,
                     "oss": oss, "ps2": ps2,
                     "count_match": cnt, "size_match": sz, "coord_match": crd})
        flag = "OK " if crd else ("~sz" if sz else ("~ct" if cnt else "DIFF"))
        print(f"[{i:2}/{len(loci)}] {lid:16} oss={len(oss)} ps2={len(ps2)}  "
              f"{flag}  sizes oss={sorted(s for s,_ in oss)} ps2={sorted(s for s,_ in ps2)}",
              file=sys.stderr)

    print("\n" + "=" * 64)
    print(f"loci evaluated (pair designed): {n_designed}/{len(loci)}")
    if n_designed:
        print(f"amplicon-count concordance : {count_match}/{n_designed} "
              f"({100*count_match/n_designed:.0f}%)")
        print(f"product-size concordance   : {size_match}/{n_designed} "
              f"({100*size_match/n_designed:.0f}%)")
        print(f"exact-coordinate concordance: {coord_match}/{n_designed} "
              f"({100*coord_match/n_designed:.0f}%)")
    print("=" * 64)

    if a.out:
        with open(a.out, "w") as fh:
            json.dump({"n_loci": len(loci), "n_designed": n_designed,
                       "count_match": count_match, "size_match": size_match,
                       "coord_match": coord_match, "tm_diff": a.tm_diff,
                       "check_size_min": a.check_size_min,
                       "check_size_max": a.check_size_max,
                       "specificity_profile": a.specificity_profile,
                       "thermo": thermo_params is not None, "rows": rows},
                      fh, indent=2, default=str)
        print(f"wrote {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
