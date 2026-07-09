"""Command-line interface for primerblast-oss.

Subcommands:
  design            region + product size  ->  primer pairs, checked for specificity
  check             primer sequences       ->  all predicted PCR products (in-silico PCR)
  multiplex         primer pool            ->  primer-dimer compatibility across the pool
  multiplex-design  multi-target FASTA     ->  design + pick a mutually compatible set
  tile              region + amplicon size ->  overlapping amplicons covering a region
  assay      gene / interval / SNP  ->  design + specificity + variants + CAPS + risk
  markers    QTL interval           ->  evenly spaced markers across the interval
  makedb     FASTA                  ->  BLAST nucleotide database
"""
from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List, Tuple

from . import __version__
from .design import DesignParams, read_fasta
from .specificity import (
    SPECIFICITY_PROFILES, SpecParams, in_silico_pcr, spec_params_for_profile,
)
from .pipeline import run_pipeline
from .tiling import design_tiling
from .tools import make_blastdb
from . import report as R
from . import outputs as OUT


# --------------------------------------------------------------------------- #
# shared helpers
# --------------------------------------------------------------------------- #
def _parse_size_ranges(text: str) -> List[Tuple[int, int]]:
    ranges = []
    for chunk in text.replace(",", " ").split():
        try:
            lo, hi = chunk.split("-")
            lo_i, hi_i = int(lo), int(hi)
        except ValueError as exc:
            raise ValueError(
                f"invalid size range '{chunk}'; use LOW-HIGH, e.g. 150-500"
            ) from exc
        if lo_i <= 0 or hi_i < lo_i:
            raise ValueError(
                f"invalid size range '{chunk}'; require 0 < LOW <= HIGH"
            )
        ranges.append((lo_i, hi_i))
    return ranges


def _spec_from_args(a) -> SpecParams:
    max_target_seqs = a.max_target_seqs
    if a.exhaustive and max_target_seqs is None:
        max_target_seqs = 50000
    return spec_params_for_profile(
        a.specificity_profile,
        max_total_mismatch=a.max_total_mismatch,
        max_3prime_mismatch=a.max_3prime_mismatch,
        three_prime_window=a.three_prime_window,
        require_3prime_terminal_match=(
            False if a.no_3prime_terminal
            else True if a.require_3prime_terminal
            else None
        ),
        min_product=a.min_product, max_product=a.max_product,
        gel_min_gap_bp=a.gel_min_gap,
        word_size=a.word_size,
        evalue=a.evalue,
        max_target_seqs=max_target_seqs,
        num_threads=a.num_threads,
        high_copy_hit_threshold=a.high_copy_hit_threshold,
        high_copy_site_threshold=a.high_copy_site_threshold,
    )


def _add_spec_args(ap: argparse.ArgumentParser) -> None:
    s = ap.add_argument_group("specificity (blast)")
    s.add_argument("--db", action="append", required=True,
                   help="BLAST nucleotide db path (repeatable for multi-db screening)")
    s.add_argument("--specificity-profile", choices=sorted(SPECIFICITY_PROFILES),
                   default="local-strict",
                   help="preset mismatch model; individual options below override it")
    s.add_argument("--max-total-mismatch", type=int)
    s.add_argument("--max-3prime-mismatch", type=int)
    s.add_argument("--three-prime-window", type=int)
    s.add_argument("--min-product", type=int, default=40)
    s.add_argument("--max-product", type=int, default=4000)
    s.add_argument("--gel-min-gap", type=int, default=50,
                   help="size gap (bp) needed to resolve two products on a gel")
    s.add_argument("--word-size", type=int, default=7)
    s.add_argument("--evalue", type=float, default=30000.0)
    s.add_argument("--max-target-seqs", type=int,
                   help="BLAST -max_target_seqs; raise this for repetitive genomes")
    s.add_argument("--num-threads", type=int, default=4,
                   help="blastn worker threads")
    s.add_argument("--high-copy-hit-threshold", type=int, default=10000,
                   help="raw BLAST HSP count reported per primer (informational)")
    s.add_argument("--high-copy-site-threshold", type=int, default=500,
                   help="priming-site count that flags a repeat-prone primer")
    s.add_argument("--exhaustive", action="store_true",
                   help="use a higher BLAST hit cap (50000 unless --max-target-seqs is set)")
    s.add_argument("--no-3prime-terminal", action="store_true",
                   help="do not require the 3'-terminal base to match")
    s.add_argument("--require-3prime-terminal", action="store_true",
                   help="require the 3'-terminal base to match, overriding the profile")
    s.add_argument("--blastn-bin")
    t = ap.add_argument_group("thermodynamics (optional, needs primer3-py)")
    t.add_argument("--genome-fasta",
                   help=".fai-indexed genome FASTA enabling thermodynamic site scoring "
                        "(subject names must match the BLAST db)")
    t.add_argument("--no-thermo", action="store_true",
                   help="disable thermodynamic scoring even if a genome is available")
    t.add_argument("--no-thermo-gate", action="store_true",
                   help="annotate Tm/dG but do not drop thermodynamically non-viable sites")
    t.add_argument("--min-anneal-tm", type=float, default=40.0,
                   help="duplex Tm floor (degC) for a site to count as priming")
    t.add_argument("--max-3p-dg", type=float, default=-5.0,
                   help="3'-end stability ΔG ceiling (kcal/mol) for a viable site")


def _thermo_setup(a, genome=None):
    """Return (genome, thermo_params, thermo_gate) from CLI args.

    `genome` may be pre-supplied (assay/markers already open one); otherwise a
    --genome-fasta path is opened. Thermo is off when disabled or unavailable."""
    from .thermo import ThermoParams, available
    if getattr(a, "no_thermo", False) or not available():
        return genome, None, True
    if genome is None and getattr(a, "genome_fasta", None):
        from .genome import Genome
        genome = Genome(a.genome_fasta)
    tp = ThermoParams(min_anneal_tm=a.min_anneal_tm, max_3p_dg=a.max_3p_dg)
    return genome, tp, (not getattr(a, "no_thermo_gate", False))


def _add_dimer_args(ap: argparse.ArgumentParser) -> None:
    d = ap.add_argument_group("primer-dimer / hairpin (optional, needs primer3-py)")
    d.add_argument("--dimer-dg-warn", type=float, default=-8.0,
                   help="flag structures with ΔG at or below this kcal/mol")
    d.add_argument("--dimer-tm-warn", type=float, default=45.0,
                   help="flag weaker structures at or above this Tm")
    d.add_argument("--dimer-dg-at-tm-warn", type=float, default=-4.0,
                   help="ΔG floor paired with --dimer-tm-warn")


def _dimer_params_from_args(a):
    from .dimers import DimerParams
    return DimerParams(
        dg_warn=a.dimer_dg_warn,
        tm_warn=a.dimer_tm_warn,
        dg_at_tm_warn=a.dimer_dg_at_tm_warn,
    )


def _add_out_args(ap: argparse.ArgumentParser, formats=("text", "json", "tsv")) -> None:
    o = ap.add_argument_group("output")
    o.add_argument("--format", choices=list(formats), default="text")
    o.add_argument("--out", help="write to file instead of stdout")


def _emit(text: str, out: str) -> None:
    if out:
        with open(out, "w") as fh:
            fh.write(text + "\n")
        print(f"wrote {out}", file=sys.stderr)
    else:
        print(text)


def _templates(a) -> List[Tuple[str, str]]:
    if getattr(a, "template_fasta", None):
        return read_fasta(a.template_fasta)
    return [(a.template_id, a.template)]


# --------------------------------------------------------------------------- #
# design
# --------------------------------------------------------------------------- #
def _cmd_design(a) -> int:
    target = None
    if a.target:
        x, y = a.target.split(",")
        target = (int(x), int(y))
    dp = DesignParams(
        product_size_ranges=_parse_size_ranges(a.product_size),
        opt_size=a.opt_size, min_size=a.min_size, max_size=a.max_size,
        opt_tm=a.opt_tm, min_tm=a.min_tm, max_tm=a.max_tm,
        min_gc=a.min_gc, max_gc=a.max_gc, num_return=a.num_return, target=target,
    )
    sp = _spec_from_args(a)
    genome, tp, gate = _thermo_setup(a)
    dimer_params = _dimer_params_from_args(a)
    outs = []
    for tid, seq in _templates(a):
        res = run_pipeline(tid, seq, a.db, design_params=dp, spec_params=sp,
                           primer3_bin=a.primer3_bin, blastn_bin=a.blastn_bin,
                           size_tolerance=a.size_tolerance,
                           genome=genome, thermo_params=tp, thermo_gate=gate,
                           dimer_params=dimer_params)
        outs.append(R.to_json(res) if a.format == "json"
                    else R.to_tsv(res) if a.format == "tsv"
                    else R.to_text(res))
    _emit(("\n" if a.format == "json" else "\n\n").join(outs), a.out)
    return 0


# --------------------------------------------------------------------------- #
# check (in-silico PCR)
# --------------------------------------------------------------------------- #
def _collect_primers(a) -> Dict[str, str]:
    primers: Dict[str, str] = {}
    if a.primers_fasta:
        for name, seq in read_fasta(a.primers_fasta):
            primers[name] = seq.upper()
    if a.forward:
        primers["F"] = a.forward.upper()
    if a.reverse:
        primers["R"] = a.reverse.upper()
    for i, spec in enumerate(a.primer or [], 1):
        if "=" in spec:
            name, seq = spec.split("=", 1)
        else:
            name, seq = f"P{i}", spec
        primers[name] = seq.upper()
    if not primers:
        raise SystemExit("check: provide primers via --forward/--reverse, "
                         "--primer, or --primers-fasta")
    return primers


def _cmd_check(a) -> int:
    primers = _collect_primers(a)
    sp = _spec_from_args(a)
    genome, tp, gate = _thermo_setup(a)
    results = [in_silico_pcr(primers, db, sp=sp, blastn_bin=a.blastn_bin,
                             genome=genome, thermo_params=tp, thermo_gate=gate)
               for db in a.db]
    if a.format == "json":
        text = json.dumps(R.insilico_to_dict(results, primers), indent=2, default=str)
    else:
        text = R.insilico_to_text(results, primers)
    _emit(text, a.out)
    return 0


# --------------------------------------------------------------------------- #
# multiplex (primer-dimer compatibility, no BLAST)
# --------------------------------------------------------------------------- #
def _cmd_multiplex(a) -> int:
    from . import dimers
    if not dimers.available():
        raise SystemExit("multiplex: needs primer3-py (pip install .[thermo])")
    primers = _collect_primers(a)
    dp = _dimer_params_from_args(a)
    res = dimers.analyze_multiplex(list(primers.items()), dp)
    if a.format == "json":
        out = {
            "n_primers": res["n_primers"], "compatible": res["compatible"],
            "n_concerning": res["n_concerning"], "worst_dg": res["worst_dg"],
            "concerning": [{"a": s.a, "b": s.b, "tm": s.tm, "dg": s.dg} for s in res["concerning"]],
        }
        _emit(json.dumps(out, indent=2), a.out)
        return 0
    lines = [f"multiplex check: {res['n_primers']} primers, "
             f"{res['n_pairs_checked']} cross-dimers evaluated",
             f"verdict: {'COMPATIBLE' if res['compatible'] else 'INCOMPATIBLE'} "
             f"({res['n_concerning']} concerning; worst ΔG {res['worst_dg']} kcal/mol)"]
    for s in res["concerning"]:
        lines.append(f"  {s.a} x {s.b}: cross-dimer Tm {s.tm} ΔG {s.dg} kcal/mol")
    _emit("\n".join(lines), a.out)
    return 0


# --------------------------------------------------------------------------- #
# multiplex-design (design + pick a mutually compatible set, one pair per target)
# --------------------------------------------------------------------------- #
def _cmd_multiplex_design(a) -> int:
    from . import dimers
    if not dimers.available():
        raise SystemExit("multiplex-design: needs primer3-py (pip install .[thermo])")

    dp = DesignParams(
        product_size_ranges=_parse_size_ranges(a.product_size),
        opt_size=a.opt_size, min_size=a.min_size, max_size=a.max_size,
        opt_tm=a.opt_tm, min_tm=a.min_tm, max_tm=a.max_tm,
        min_gc=a.min_gc, max_gc=a.max_gc, num_return=a.num_return,
    )
    sp = _spec_from_args(a)
    genome, tp, gate = _thermo_setup(a)
    dimer_params = _dimer_params_from_args(a)

    templates = _templates(a)
    candidates: List[Tuple[str, List[Tuple[str, str]]]] = []
    per_target: Dict[str, List] = {}
    for tid, seq in templates:
        res = run_pipeline(tid, seq, a.db, design_params=dp, spec_params=sp,
                           primer3_bin=a.primer3_bin, blastn_bin=a.blastn_bin,
                           size_tolerance=a.size_tolerance,
                           genome=genome, thermo_params=tp, thermo_gate=gate,
                           dimer_params=dimer_params)
        pairs = res.pairs
        if a.require_specific:
            good = [p for p in pairs if p.specificity.get("rank") in ("A", "B")]
            pairs = good
        pairs = pairs[:a.candidates_per_target]
        per_target[tid] = pairs
        candidates.append((tid, [(p.forward, p.reverse) for p in pairs]))

    selection = dimers.select_multiplex_set(candidates, dimer_params)
    if selection is None:
        raise SystemExit("multiplex-design: dimer analysis unavailable")

    # describe the chosen set and re-check it every-vs-every for a final ΔG
    chosen_primers: List[Tuple[str, str]] = []
    enriched = []
    for sel in selection["selection"]:
        tid = sel["target"]
        idx = sel.get("candidate_index")
        if idx is None:
            enriched.append({"target": tid, "assigned": False})
            continue
        pair = per_target[tid][idx]
        chosen_primers += [(f"{tid}_F", pair.forward), (f"{tid}_R", pair.reverse)]
        enriched.append({
            "target": tid, "assigned": True, "candidate_index": idx,
            "forward": pair.forward, "reverse": pair.reverse,
            "product_size": pair.product_size,
            "rank": pair.specificity.get("rank"),
            "score": pair.specificity.get("score"),
        })
    pool = dimers.analyze_multiplex(chosen_primers, dimer_params) if chosen_primers else None

    if a.format == "json":
        out = {
            "n_targets": selection["n_targets"],
            "n_assigned": selection["n_assigned"],
            "complete": selection["complete"],
            "unassigned": selection["unassigned"],
            "candidates_per_target": a.candidates_per_target,
            "selection": enriched,
            "pool_check": (None if pool is None else {
                "compatible": pool["compatible"], "n_concerning": pool["n_concerning"],
                "worst_dg": pool["worst_dg"],
                "concerning": [{"a": s.a, "b": s.b, "tm": s.tm, "dg": s.dg}
                               for s in pool["concerning"]],
            }),
        }
        _emit(json.dumps(out, indent=2, default=str), a.out)
        return 0

    verdict = ("COMPLETE" if selection["complete"]
               else f"PARTIAL ({selection['n_assigned']}/{selection['n_targets']} assigned)")
    lines = [f"multiplex design: {selection['n_targets']} targets, "
             f"up to {a.candidates_per_target} candidates each",
             f"verdict: {verdict}"]
    for e in enriched:
        if not e["assigned"]:
            lines.append(f"  {e['target']}: UNASSIGNED "
                         f"(no candidate avoids cross-dimers with the rest)")
        else:
            lines.append(f"  {e['target']}: cand #{e['candidate_index']}  "
                         f"{e['product_size']}bp  rank {e['rank']}  "
                         f"F {e['forward']} / R {e['reverse']}")
    if pool is not None:
        lines.append(f"pool cross-dimer check: "
                     f"{'COMPATIBLE' if pool['compatible'] else 'INCOMPATIBLE'} "
                     f"({pool['n_concerning']} concerning; worst ΔG {pool['worst_dg']} kcal/mol)")
        for s in pool["concerning"]:
            lines.append(f"    {s.a} x {s.b}: Tm {s.tm} ΔG {s.dg} kcal/mol")
    if selection["unassigned"]:
        lines.append("note: raise --candidates-per-target or relax --product-size "
                     "to give the selector more room.")
    _emit("\n".join(lines), a.out)
    return 0


# --------------------------------------------------------------------------- #
# tile
# --------------------------------------------------------------------------- #
def _cmd_tile(a) -> int:
    region = None
    if a.region:
        x, y = a.region.split(",")
        region = (int(x), int(y))
    dp = DesignParams(
        opt_size=a.opt_size, min_size=a.min_size, max_size=a.max_size,
        opt_tm=a.opt_tm, min_tm=a.min_tm, max_tm=a.max_tm,
        min_gc=a.min_gc, max_gc=a.max_gc,
    )
    sp = _spec_from_args(a)
    genome, tp, gate = _thermo_setup(a)
    dimer_params = _dimer_params_from_args(a)
    outs = []
    for tid, seq in _templates(a):
        tiles = design_tiling(
            tid, seq, a.db, region=region,
            amplicon_min=a.amplicon_min, amplicon_max=a.amplicon_max,
            overlap=a.overlap, design_params=dp, spec_params=sp,
            primer3_bin=a.primer3_bin, blastn_bin=a.blastn_bin,
            size_tolerance=a.size_tolerance,
            candidates_per_tile=a.candidates_per_tile,
            genome=genome, thermo_params=tp, thermo_gate=gate,
            dimer_params=dimer_params,
        )
        from .design import clean_sequence
        reg = region or (0, len(clean_sequence(seq)) - 1)
        if a.format == "json":
            outs.append(json.dumps(R.tiling_to_dict(tiles, tid, reg, a.db), indent=2, default=str))
        else:
            outs.append(R.tiling_to_text(tiles, tid, reg))
    _emit(("\n" if a.format == "json" else "\n\n").join(outs), a.out)
    return 0


# --------------------------------------------------------------------------- #
# assay / markers (breeding pipeline)
# --------------------------------------------------------------------------- #
def _load_variants(vcf_path):
    if not vcf_path:
        return []
    from .vcf import parse_vcf
    return parse_vcf(vcf_path)


def _render_assay(result, fmt: str, provenance=None) -> str:
    pairs = result.get("pairs", [])
    if fmt == "json":
        import json as _json
        if provenance is not None:
            result = {**result, "provenance": provenance}
        return _json.dumps(result, indent=2, default=str)
    if fmt == "csv":
        return OUT.pairs_to_csv(pairs)
    if fmt == "bed":
        return OUT.products_to_bed(pairs)
    if fmt == "order":
        return OUT.order_table(pairs)
    if fmt == "html":
        ctx = {
            "title": f"Assay: {result['target']['name']}",
            "template": f"{result['target']['chrom']}:{result['target']['start']}-{result['target']['end']}",
            "databases": result["databases"],
            "generated": (provenance or {}).get("generated", ""),
            "params": {}, "provenance": provenance or {}, "pairs": pairs,
        }
        return OUT.html_report(ctx)
    return R.assay_to_text(result)


def _cmd_assay(a) -> int:
    from .genome import Genome
    from .regions import resolve_gene, resolve_interval, resolve_snp
    from .assay import run_assay
    from .provenance import make_manifest

    genome = Genome(a.genome)
    caps_snp = None
    if a.gene:
        region = resolve_gene(a.gff3, a.gene, feature=a.gene_feature, flank=0)
    elif a.interval:
        chrom, span = a.interval.split(":")
        s, e = span.split("-")
        region = resolve_interval(chrom, int(s), int(e), name=a.name or a.interval)
    elif a.snp:
        chrom, pos = a.snp.split(":")
        region = resolve_snp(chrom, int(pos), flank=a.flank or 250, name=a.name)
        if a.alt:
            caps_snp = {"genomic_pos": int(pos), "alt": a.alt}
    else:
        raise SystemExit("assay: give --gene (with --gff3), --interval, or --snp")

    dp = DesignParams(product_size_ranges=_parse_size_ranges(a.product_size),
                      opt_tm=a.opt_tm, min_tm=a.min_tm, max_tm=a.max_tm,
                      min_gc=a.min_gc, max_gc=a.max_gc, num_return=a.num_return)
    sp = _spec_from_args(a)
    variants = _load_variants(a.vcf)
    _g, tp, gate = _thermo_setup(a, genome=genome)
    dimer_params = _dimer_params_from_args(a)

    result = run_assay(region, genome, a.db, flank=a.flank, design_params=dp,
                       spec_params=sp, variants=variants, caps_snp=caps_snp,
                       primer3_bin=a.primer3_bin, blastn_bin=a.blastn_bin,
                       thermo_params=tp, thermo_gate=gate,
                       dimer_params=dimer_params)
    prov = make_manifest({"design": dp.__dict__, "spec": sp.__dict__, "flank": a.flank},
                         a.db, template_info=result["target"])
    _emit(_render_assay(result, a.format, prov), a.out)
    return 0


def _cmd_markers(a) -> int:
    from .genome import Genome
    from .regions import resolve_interval
    from .assay import design_qtl_markers
    import json as _json

    genome = Genome(a.genome)
    chrom, span = a.interval.split(":")
    s, e = span.split("-")
    qtl = resolve_interval(chrom, int(s), int(e), name=a.name or "QTL")
    dp = DesignParams(product_size_ranges=_parse_size_ranges(a.product_size),
                      opt_tm=a.opt_tm, min_tm=a.min_tm, max_tm=a.max_tm,
                      num_return=a.num_return)
    sp = _spec_from_args(a)
    markers = design_qtl_markers(qtl, genome, a.db, n_markers=a.n_markers,
                                 spacing=a.spacing, marker_flank=a.marker_flank,
                                 design_params=dp, spec_params=sp,
                                 primer3_bin=a.primer3_bin, blastn_bin=a.blastn_bin)
    if a.format == "json":
        _emit(_json.dumps(markers, indent=2, default=str), a.out)
    else:
        lines = [f"QTL {a.interval}: {len(markers)} markers"]
        for m in markers:
            if m.get("pairs"):
                p = m["pairs"][0]
                lines.append(f"  {m['marker']} @ {m.get('anchor')}: risk {p['risk']}  "
                             f"{p['product_size']}bp  F {p['forward']} / R {p['reverse']}")
            else:
                lines.append(f"  {m['marker']} @ {m.get('anchor')}: "
                             f"{m.get('error','no pair found')}")
        _emit("\n".join(lines), a.out)
    return 0


# --------------------------------------------------------------------------- #
# makedb
# --------------------------------------------------------------------------- #
def _cmd_makedb(a) -> int:
    out = make_blastdb(a.infile, out=a.out_db, title=a.title,
                       parse_seqids=not a.no_parse_seqids,
                       makeblastdb_bin=a.makeblastdb_bin)
    print(f"built database: {out}")
    return 0


# --------------------------------------------------------------------------- #
# parser
# --------------------------------------------------------------------------- #
def _add_template_args(sp) -> None:
    g = sp.add_argument_group("template")
    src = g.add_mutually_exclusive_group(required=True)
    src.add_argument("--template", help="template DNA sequence (string)")
    src.add_argument("--template-fasta", help="FASTA file with one or more templates")
    g.add_argument("--template-id", default="template")


def _add_design_knobs(sp) -> None:
    d = sp.add_argument_group("primer3")
    d.add_argument("--opt-size", type=int, default=20)
    d.add_argument("--min-size", type=int, default=18)
    d.add_argument("--max-size", type=int, default=25)
    d.add_argument("--opt-tm", type=float, default=60.0)
    d.add_argument("--min-tm", type=float, default=57.0)
    d.add_argument("--max-tm", type=float, default=63.0)
    d.add_argument("--min-gc", type=float, default=20.0)
    d.add_argument("--max-gc", type=float, default=80.0)
    d.add_argument("--size-tolerance", type=int, default=10)
    d.add_argument("--primer3-bin")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="primerblast-oss",
        description="Local, offline Primer-BLAST for plant breeding: primer design, "
                    "in-silico PCR, multiplex-set design, region tiling, and full "
                    "gene/interval/SNP assays with multi-reference specificity.",
        epilog="Run 'primerblast-oss <subcommand> --help' for a subcommand's options.")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # design
    d = sub.add_parser("design", help="design primer pairs and check specificity")
    _add_template_args(d)
    d.add_argument("--product-size", default="70-1000",
                   help="range(s), e.g. '150-500' or '150-500,500-1000'")
    d.add_argument("--num-return", type=int, default=10)
    d.add_argument("--target", help="focus region on template: start,length (0-based)")
    _add_design_knobs(d)
    _add_spec_args(d)
    _add_dimer_args(d)
    _add_out_args(d, formats=("text", "json", "tsv"))
    d.set_defaults(func=_cmd_design)

    # check
    c = sub.add_parser("check", help="in-silico PCR: predict products for given primers")
    c.add_argument("--forward")
    c.add_argument("--reverse")
    c.add_argument("--primer", action="append", help="NAME=SEQ or SEQ (repeatable)")
    c.add_argument("--primers-fasta")
    _add_spec_args(c)
    _add_out_args(c, formats=("text", "json"))
    c.set_defaults(func=_cmd_check)

    # multiplex (primer-dimer compatibility of a primer pool; no BLAST)
    mx = sub.add_parser("multiplex",
                        help="check primer-dimer compatibility across a pool of primers")
    mx.add_argument("--forward")
    mx.add_argument("--reverse")
    mx.add_argument("--primer", action="append", help="NAME=SEQ or SEQ (repeatable)")
    mx.add_argument("--primers-fasta")
    _add_dimer_args(mx)
    _add_out_args(mx, formats=("text", "json"))
    mx.set_defaults(func=_cmd_multiplex)

    # multiplex-design (design each target, then pick one mutually compatible pair each)
    md = sub.add_parser(
        "multiplex-design",
        help="design primers for several targets and pick a compatible set (one pair each)")
    _add_template_args(md)
    md.add_argument("--product-size", default="80-300",
                    help="range(s), e.g. '80-300' (short amplicons suit multiplex)")
    md.add_argument("--num-return", type=int, default=10,
                    help="pairs primer3 returns per target before specificity ranking")
    md.add_argument("--candidates-per-target", type=int, default=5,
                    help="top-ranked pairs per target offered to the set selector")
    md.add_argument("--require-specific", action="store_true",
                    help="require rank A/B candidates; leave targets without one unassigned")
    _add_design_knobs(md)
    _add_spec_args(md)
    _add_dimer_args(md)
    _add_out_args(md, formats=("text", "json"))
    md.set_defaults(func=_cmd_multiplex_design)

    # tile
    t = sub.add_parser("tile", help="cover a whole region with overlapping amplicons")
    _add_template_args(t)
    t.add_argument("--region", help="sub-region to cover: start,end (0-based inclusive)")
    t.add_argument("--amplicon-min", type=int, default=400)
    t.add_argument("--amplicon-max", type=int, default=800)
    t.add_argument("--overlap", type=int, default=40)
    t.add_argument("--candidates-per-tile", type=int, default=8,
                   help="primer pairs evaluated per window (lower = faster)")
    _add_design_knobs(t)
    _add_spec_args(t)
    _add_dimer_args(t)
    _add_out_args(t, formats=("text", "json"))
    t.set_defaults(func=_cmd_tile)

    # assay (breeding pipeline: gene/interval/SNP -> primers + variants + CAPS + risk)
    y = sub.add_parser("assay", help="full assay: design + specificity + variants + CAPS + risk")
    tg = y.add_argument_group("target (choose one)")
    tg.add_argument("--gene", help="gene id (requires --gff3)")
    tg.add_argument("--gene-feature", default="cds", choices=["gene", "mrna", "exon", "cds"])
    tg.add_argument("--interval", help="chrom:start-end (1-based)")
    tg.add_argument("--snp", help="chrom:pos (design a CAPS amplicon spanning it)")
    tg.add_argument("--alt", help="alt allele base at --snp (enables CAPS scan)")
    tg.add_argument("--name", help="target name")
    y.add_argument("--genome", required=True, help="reference FASTA (.fai indexed)")
    y.add_argument("--gff3", help="GFF3 annotation (for --gene)")
    y.add_argument("--vcf", help="VCF of variants (SNP-under-primer + amplicon SNPs)")
    y.add_argument("--flank", type=int, default=200)
    y.add_argument("--product-size", default="100-800")
    y.add_argument("--num-return", type=int, default=10)
    _add_design_knobs(y)
    _add_spec_args(y)
    _add_dimer_args(y)
    _add_out_args(y, formats=("text", "json", "csv", "bed", "order", "html"))
    y.set_defaults(func=_cmd_assay)

    # markers (QTL interval -> evenly spaced markers)
    k = sub.add_parser("markers", help="design evenly spaced markers across a QTL interval")
    k.add_argument("--interval", required=True, help="chrom:start-end QTL interval")
    k.add_argument("--genome", required=True)
    k.add_argument("--name")
    k.add_argument("--n-markers", type=int, default=0)
    k.add_argument("--spacing", type=int, default=0, help="bp between markers (alt to --n-markers)")
    k.add_argument("--marker-flank", type=int, default=300)
    k.add_argument("--product-size", default="100-500")
    k.add_argument("--num-return", type=int, default=6)
    _add_design_knobs(k)
    _add_spec_args(k)
    _add_out_args(k, formats=("text", "json"))
    k.set_defaults(func=_cmd_markers)

    # makedb
    m = sub.add_parser("makedb", help="build a BLAST nucleotide database from FASTA")
    m.add_argument("infile", help="input FASTA")
    m.add_argument("--out-db", help="output db prefix")
    m.add_argument("--title")
    m.add_argument("--no-parse-seqids", action="store_true")
    m.add_argument("--makeblastdb-bin")
    m.set_defaults(func=_cmd_makedb)

    return ap


def main(argv=None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, KeyError, RuntimeError) as exc:
        ap.error(str(exc))


if __name__ == "__main__":
    raise SystemExit(main())
