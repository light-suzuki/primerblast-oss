"""Command-line interface for primerblast-oss."""
from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List, Mapping, Optional, Tuple

from . import __version__
from .design import DesignParams, read_fasta
from .specificity import (
    SPECIFICITY_PROFILES,
    SpecParams,
    in_silico_pcr,
    spec_params_for_profile,
)
from .pipeline import (
    resolve_genome_for_database,
    run_pipeline,
    thermo_metadata,
)
from .tiling import design_tiling
from .tools import make_blastdb
from . import outputs as OUT
from . import report as R


def _parse_size_ranges(text: str) -> List[Tuple[int, int]]:
    ranges = []
    for chunk in text.replace(",", " ").split():
        try:
            low, high = chunk.split("-")
            low_i, high_i = int(low), int(high)
        except ValueError as error:
            raise ValueError(
                "invalid size range '%s'; use LOW-HIGH, e.g. 150-500" % chunk
            ) from error
        if low_i <= 0 or high_i < low_i:
            raise ValueError(
                "invalid size range '%s'; require 0 < LOW <= HIGH" % chunk)
        ranges.append((low_i, high_i))
    return ranges


def _spec_from_args(arguments) -> SpecParams:
    max_target_seqs = arguments.max_target_seqs
    if arguments.exhaustive and max_target_seqs is None:
        max_target_seqs = 50000
    return spec_params_for_profile(
        arguments.specificity_profile,
        max_total_mismatch=arguments.max_total_mismatch,
        max_3prime_mismatch=arguments.max_3prime_mismatch,
        three_prime_window=arguments.three_prime_window,
        require_3prime_terminal_match=(
            False if arguments.no_3prime_terminal
            else True if arguments.require_3prime_terminal
            else None
        ),
        min_product=arguments.min_product,
        max_product=arguments.max_product,
        gel_min_gap_bp=arguments.gel_min_gap,
        word_size=arguments.word_size,
        evalue=arguments.evalue,
        max_target_seqs=max_target_seqs,
        num_threads=arguments.num_threads,
        high_copy_hit_threshold=arguments.high_copy_hit_threshold,
        high_copy_site_threshold=arguments.high_copy_site_threshold,
    )


def _add_spec_args(parser: argparse.ArgumentParser) -> None:
    specificity = parser.add_argument_group("specificity (blast)")
    specificity.add_argument(
        "--db", action="append", required=True,
        help="BLAST nucleotide DB path (repeatable for multi-DB screening)")
    specificity.add_argument(
        "--specificity-profile", choices=sorted(SPECIFICITY_PROFILES),
        default="local-strict",
        help="preset mismatch model; individual options below override it")
    specificity.add_argument("--max-total-mismatch", type=int)
    specificity.add_argument("--max-3prime-mismatch", type=int)
    specificity.add_argument("--three-prime-window", type=int)
    specificity.add_argument("--min-product", type=int, default=40)
    specificity.add_argument("--max-product", type=int, default=4000)
    specificity.add_argument(
        "--gel-min-gap", type=int, default=50,
        help="size gap (bp) needed to resolve two products on a gel")
    specificity.add_argument("--word-size", type=int, default=7)
    specificity.add_argument("--evalue", type=float, default=30000.0)
    specificity.add_argument(
        "--max-target-seqs", type=int,
        help="BLAST -max_target_seqs; raise this for repetitive genomes")
    specificity.add_argument(
        "--num-threads", type=int, default=4, help="blastn worker threads")
    specificity.add_argument(
        "--high-copy-hit-threshold", type=int, default=10000,
        help="raw BLAST HSP count that makes a search repeat-limited")
    specificity.add_argument(
        "--high-copy-site-threshold", type=int, default=500,
        help="priming-site count that makes a search repeat-limited")
    specificity.add_argument(
        "--exhaustive", action="store_true",
        help="use a higher BLAST hit cap (50000 unless explicitly set)")
    specificity.add_argument(
        "--no-3prime-terminal", action="store_true",
        help="do not require the 3'-terminal base to match")
    specificity.add_argument(
        "--require-3prime-terminal", action="store_true",
        help="require the 3'-terminal base to match, overriding the profile")
    specificity.add_argument("--blastn-bin")

    thermo = parser.add_argument_group(
        "thermodynamics (optional, needs primer3-py)")
    thermo.add_argument(
        "--db-genome", action="append", default=[], metavar="DB=FASTA",
        help="associate one BLAST DB with its .fai-indexed FASTA; repeat per DB")
    thermo.add_argument(
        "--genome-fasta",
        help="legacy FASTA association for the first/design DB only; never reused "
             "for secondary DBs")
    thermo.add_argument(
        "--no-thermo", action="store_true",
        help="disable thermodynamic site scoring")
    thermo.add_argument(
        "--no-thermo-gate", action="store_true",
        help="annotate Tm/dG but do not drop non-viable sites")
    thermo.add_argument(
        "--min-anneal-tm", type=float, default=40.0,
        help="duplex Tm floor (degC) for a site to count as priming")
    thermo.add_argument(
        "--max-3p-dg", type=float, default=-5.0,
        help="3'-end stability ΔG ceiling (kcal/mol) for a viable site")


def _parse_db_genome_specs(specifications: List[str],
                           databases: List[str]) -> Dict[str, str]:
    mappings: Dict[str, str] = {}
    for specification in specifications or []:
        if "=" not in specification:
            raise ValueError(
                "invalid --db-genome '%s'; use exact DB=FASTA" % specification)
        database, fasta = specification.split("=", 1)
        database, fasta = database.strip(), fasta.strip()
        if not database or not fasta:
            raise ValueError(
                "invalid --db-genome '%s'; DB and FASTA are required" % specification)
        if database not in databases:
            raise ValueError(
                "--db-genome DB '%s' does not exactly match a supplied --db" % database)
        if database in mappings and mappings[database] != fasta:
            raise ValueError(
                "multiple different FASTA files supplied for DB '%s'" % database)
        mappings[database] = fasta
    return mappings


def _thermo_setup(arguments, design_genome=None):
    """Return ``(genomes_by_db, ThermoParams|None, gate)``.

    An association is always explicit. ``design_genome`` and legacy
    ``--genome-fasta`` apply only to the first DB. ``--db-genome`` may override
    that association and is the required interface for secondary assemblies.
    """
    from .thermo import ThermoParams, available

    databases = list(arguments.db)
    fasta_by_db = _parse_db_genome_specs(
        getattr(arguments, "db_genome", []), databases)
    genomes_by_db: Dict[str, object] = {}
    if design_genome is not None and databases:
        genomes_by_db[databases[0]] = design_genome

    legacy_fasta = getattr(arguments, "genome_fasta", None)
    if legacy_fasta and databases and databases[0] not in fasta_by_db:
        fasta_by_db[databases[0]] = legacy_fasta

    if fasta_by_db:
        from .genome import Genome
        for database, fasta in fasta_by_db.items():
            genomes_by_db[database] = Genome(fasta)

    if getattr(arguments, "no_thermo", False):
        return genomes_by_db, None, True
    if not available():
        if genomes_by_db:
            print(
                "warning: primer3-py is unavailable; thermodynamic scoring skipped",
                file=sys.stderr,
            )
        return genomes_by_db, None, True

    params = ThermoParams(
        min_anneal_tm=arguments.min_anneal_tm,
        max_3p_dg=arguments.max_3p_dg,
    )
    gate = not getattr(arguments, "no_thermo_gate", False)
    missing = [database for database in databases if database not in genomes_by_db]
    if missing:
        print(
            "warning: no associated FASTA for DB(s): %s; thermodynamic scoring "
            "will be skipped for those DBs. Use --db-genome DB=FASTA."
            % ", ".join(missing),
            file=sys.stderr,
        )
    return genomes_by_db, params, gate


def _genome_paths(genomes_by_db: Mapping[str, object]) -> Dict[str, Optional[str]]:
    return {
        database: getattr(genome, "fasta", None)
        for database, genome in genomes_by_db.items()
    }


def _add_dimer_args(parser: argparse.ArgumentParser) -> None:
    dimers = parser.add_argument_group(
        "primer-dimer / hairpin (optional, needs primer3-py)")
    dimers.add_argument("--dimer-dg-warn", type=float, default=-8.0)
    dimers.add_argument("--dimer-tm-warn", type=float, default=45.0)
    dimers.add_argument("--dimer-dg-at-tm-warn", type=float, default=-4.0)


def _dimer_params_from_args(arguments):
    from .dimers import DimerParams
    return DimerParams(
        dg_warn=arguments.dimer_dg_warn,
        tm_warn=arguments.dimer_tm_warn,
        dg_at_tm_warn=arguments.dimer_dg_at_tm_warn,
    )


def _add_out_args(parser: argparse.ArgumentParser,
                  formats=("text", "json", "tsv")) -> None:
    output = parser.add_argument_group("output")
    output.add_argument("--format", choices=list(formats), default="text")
    output.add_argument("--out", help="write to file instead of stdout")


def _emit(text: str, output_path: Optional[str]) -> None:
    if output_path:
        with open(output_path, "w") as handle:
            handle.write(text + "\n")
        print("wrote %s" % output_path, file=sys.stderr)
    else:
        print(text)


def _templates(arguments) -> List[Tuple[str, str]]:
    if getattr(arguments, "template_fasta", None):
        return read_fasta(arguments.template_fasta)
    return [(arguments.template_id, arguments.template)]


def _cmd_design(arguments) -> int:
    target = None
    if arguments.target:
        start, length = arguments.target.split(",")
        target = (int(start), int(length))
    design = DesignParams(
        product_size_ranges=_parse_size_ranges(arguments.product_size),
        opt_size=arguments.opt_size,
        min_size=arguments.min_size,
        max_size=arguments.max_size,
        opt_tm=arguments.opt_tm,
        min_tm=arguments.min_tm,
        max_tm=arguments.max_tm,
        min_gc=arguments.min_gc,
        max_gc=arguments.max_gc,
        num_return=arguments.num_return,
        target=target,
    )
    specificity = _spec_from_args(arguments)
    genomes_by_db, thermo_params, thermo_gate = _thermo_setup(arguments)
    dimer_params = _dimer_params_from_args(arguments)
    outputs = []
    for template_id, sequence in _templates(arguments):
        result = run_pipeline(
            template_id,
            sequence,
            arguments.db,
            design_params=design,
            spec_params=specificity,
            primer3_bin=arguments.primer3_bin,
            blastn_bin=arguments.blastn_bin,
            size_tolerance=arguments.size_tolerance,
            genomes_by_db=genomes_by_db,
            thermo_params=thermo_params,
            thermo_gate=thermo_gate,
            dimer_params=dimer_params,
        )
        outputs.append(
            R.to_json(result) if arguments.format == "json"
            else R.to_tsv(result) if arguments.format == "tsv"
            else R.to_text(result)
        )
    _emit(("\n" if arguments.format == "json" else "\n\n").join(outputs),
          arguments.out)
    return 0


def _collect_primers(arguments) -> Dict[str, str]:
    primers: Dict[str, str] = {}
    if arguments.primers_fasta:
        for name, sequence in read_fasta(arguments.primers_fasta):
            primers[name] = sequence.upper()
    if arguments.forward:
        primers["F"] = arguments.forward.upper()
    if arguments.reverse:
        primers["R"] = arguments.reverse.upper()
    for index, specification in enumerate(arguments.primer or [], 1):
        if "=" in specification:
            name, sequence = specification.split("=", 1)
        else:
            name, sequence = "P%s" % index, specification
        primers[name] = sequence.upper()
    if not primers:
        raise SystemExit(
            "check: provide primers via --forward/--reverse, --primer, or --primers-fasta")
    return primers


def _cmd_check(arguments) -> int:
    primers = _collect_primers(arguments)
    specificity = _spec_from_args(arguments)
    genomes_by_db, thermo_params, thermo_gate = _thermo_setup(arguments)
    results = []
    for database in arguments.db:
        genome, association = resolve_genome_for_database(
            database, arguments.db, genomes_by_db=genomes_by_db)
        result = in_silico_pcr(
            primers,
            database,
            sp=specificity,
            blastn_bin=arguments.blastn_bin,
            genome=genome,
            thermo_params=thermo_params,
            thermo_gate=thermo_gate,
        )
        result.update(thermo_metadata(
            genome, thermo_params, thermo_gate, association))
        results.append(result)
    if arguments.format == "json":
        text = json.dumps(
            R.insilico_to_dict(results, primers), indent=2, default=str)
    else:
        text = R.insilico_to_text(results, primers)
    _emit(text, arguments.out)
    return 0


def _cmd_multiplex(arguments) -> int:
    from . import dimers
    if not dimers.available():
        raise SystemExit("multiplex: needs primer3-py (pip install .[thermo])")
    primers = _collect_primers(arguments)
    params = _dimer_params_from_args(arguments)
    result = dimers.analyze_multiplex(list(primers.items()), params)
    if arguments.format == "json":
        output = {
            "n_primers": result["n_primers"],
            "compatible": result["compatible"],
            "n_concerning": result["n_concerning"],
            "worst_dg": result["worst_dg"],
            "concerning": [
                {"a": structure.a, "b": structure.b,
                 "tm": structure.tm, "dg": structure.dg}
                for structure in result["concerning"]
            ],
        }
        _emit(json.dumps(output, indent=2), arguments.out)
        return 0
    lines = [
        "multiplex check: %s primers, %s cross-dimers evaluated" % (
            result["n_primers"], result["n_pairs_checked"]),
        "verdict: %s (%s concerning; worst ΔG %s kcal/mol)" % (
            "COMPATIBLE" if result["compatible"] else "INCOMPATIBLE",
            result["n_concerning"], result["worst_dg"]),
    ]
    for structure in result["concerning"]:
        lines.append("  %s x %s: cross-dimer Tm %s ΔG %s kcal/mol" % (
            structure.a, structure.b, structure.tm, structure.dg))
    _emit("\n".join(lines), arguments.out)
    return 0


def _cmd_multiplex_design(arguments) -> int:
    from . import dimers
    if not dimers.available():
        raise SystemExit(
            "multiplex-design: needs primer3-py (pip install .[thermo])")
    design = DesignParams(
        product_size_ranges=_parse_size_ranges(arguments.product_size),
        opt_size=arguments.opt_size,
        min_size=arguments.min_size,
        max_size=arguments.max_size,
        opt_tm=arguments.opt_tm,
        min_tm=arguments.min_tm,
        max_tm=arguments.max_tm,
        min_gc=arguments.min_gc,
        max_gc=arguments.max_gc,
        num_return=arguments.num_return,
    )
    specificity = _spec_from_args(arguments)
    genomes_by_db, thermo_params, thermo_gate = _thermo_setup(arguments)
    dimer_params = _dimer_params_from_args(arguments)

    candidates: List[Tuple[str, List[Tuple[str, str]]]] = []
    per_target: Dict[str, List] = {}
    for template_id, sequence in _templates(arguments):
        result = run_pipeline(
            template_id,
            sequence,
            arguments.db,
            design_params=design,
            spec_params=specificity,
            primer3_bin=arguments.primer3_bin,
            blastn_bin=arguments.blastn_bin,
            size_tolerance=arguments.size_tolerance,
            genomes_by_db=genomes_by_db,
            thermo_params=thermo_params,
            thermo_gate=thermo_gate,
            dimer_params=dimer_params,
        )
        pairs = result.pairs
        if arguments.require_specific:
            pairs = [
                pair for pair in pairs
                if pair.specificity.get("rank") in ("A", "B")
            ]
        pairs = pairs[:arguments.candidates_per_target]
        per_target[template_id] = pairs
        candidates.append((
            template_id, [(pair.forward, pair.reverse) for pair in pairs]))

    selection = dimers.select_multiplex_set(candidates, dimer_params)
    if selection is None:
        raise SystemExit("multiplex-design: dimer analysis unavailable")

    chosen_primers: List[Tuple[str, str]] = []
    enriched = []
    for selected in selection["selection"]:
        template_id = selected["target"]
        candidate_index = selected.get("candidate_index")
        if candidate_index is None:
            enriched.append({"target": template_id, "assigned": False})
            continue
        pair = per_target[template_id][candidate_index]
        chosen_primers.extend([
            ("%s_F" % template_id, pair.forward),
            ("%s_R" % template_id, pair.reverse),
        ])
        enriched.append({
            "target": template_id,
            "assigned": True,
            "candidate_index": candidate_index,
            "forward": pair.forward,
            "reverse": pair.reverse,
            "product_size": pair.product_size,
            "rank": pair.specificity.get("rank"),
            "score": pair.specificity.get("score"),
            "search_completeness": pair.specificity.get("search_completeness"),
        })
    pool = (
        dimers.analyze_multiplex(chosen_primers, dimer_params)
        if chosen_primers else None
    )

    if arguments.format == "json":
        output = {
            "n_targets": selection["n_targets"],
            "n_assigned": selection["n_assigned"],
            "complete": selection["complete"],
            "unassigned": selection["unassigned"],
            "candidates_per_target": arguments.candidates_per_target,
            "thermo_genomes": _genome_paths(genomes_by_db),
            "selection": enriched,
            "pool_check": None if pool is None else {
                "compatible": pool["compatible"],
                "n_concerning": pool["n_concerning"],
                "worst_dg": pool["worst_dg"],
                "concerning": [
                    {"a": structure.a, "b": structure.b,
                     "tm": structure.tm, "dg": structure.dg}
                    for structure in pool["concerning"]
                ],
            },
        }
        _emit(json.dumps(output, indent=2, default=str), arguments.out)
        return 0

    verdict = (
        "COMPLETE" if selection["complete"]
        else "PARTIAL (%s/%s assigned)" % (
            selection["n_assigned"], selection["n_targets"])
    )
    lines = [
        "multiplex design: %s targets, up to %s candidates each" % (
            selection["n_targets"], arguments.candidates_per_target),
        "verdict: %s" % verdict,
    ]
    for entry in enriched:
        if not entry["assigned"]:
            lines.append("  %s: UNASSIGNED" % entry["target"])
        else:
            lines.append(
                "  %s: cand #%s  %sbp  rank %s  F %s / R %s" % (
                    entry["target"], entry["candidate_index"],
                    entry["product_size"], entry["rank"],
                    entry["forward"], entry["reverse"]))
    if pool is not None:
        lines.append("pool cross-dimer check: %s (%s concerning; worst ΔG %s)" % (
            "COMPATIBLE" if pool["compatible"] else "INCOMPATIBLE",
            pool["n_concerning"], pool["worst_dg"]))
    _emit("\n".join(lines), arguments.out)
    return 0


def _cmd_tile(arguments) -> int:
    region = None
    if arguments.region:
        start, end = arguments.region.split(",")
        region = (int(start), int(end))
    design = DesignParams(
        opt_size=arguments.opt_size,
        min_size=arguments.min_size,
        max_size=arguments.max_size,
        opt_tm=arguments.opt_tm,
        min_tm=arguments.min_tm,
        max_tm=arguments.max_tm,
        min_gc=arguments.min_gc,
        max_gc=arguments.max_gc,
    )
    specificity = _spec_from_args(arguments)
    genomes_by_db, thermo_params, thermo_gate = _thermo_setup(arguments)
    dimer_params = _dimer_params_from_args(arguments)
    outputs = []
    for template_id, sequence in _templates(arguments):
        tiles = design_tiling(
            template_id,
            sequence,
            arguments.db,
            region=region,
            amplicon_min=arguments.amplicon_min,
            amplicon_max=arguments.amplicon_max,
            overlap=arguments.overlap,
            design_params=design,
            spec_params=specificity,
            primer3_bin=arguments.primer3_bin,
            blastn_bin=arguments.blastn_bin,
            size_tolerance=arguments.size_tolerance,
            candidates_per_tile=arguments.candidates_per_tile,
            genomes_by_db=genomes_by_db,
            thermo_params=thermo_params,
            thermo_gate=thermo_gate,
            dimer_params=dimer_params,
        )
        from .design import clean_sequence
        requested_region = region or (0, len(clean_sequence(sequence)) - 1)
        if arguments.format == "json":
            output = R.tiling_to_dict(
                tiles, template_id, requested_region, arguments.db)
            output["thermo_genomes"] = _genome_paths(genomes_by_db)
            outputs.append(json.dumps(output, indent=2, default=str))
        else:
            outputs.append(R.tiling_to_text(
                tiles, template_id, requested_region))
    _emit(("\n" if arguments.format == "json" else "\n\n").join(outputs),
          arguments.out)
    return 0


def _load_variants(vcf_path):
    if not vcf_path:
        return []
    from .vcf import parse_vcf
    return parse_vcf(vcf_path)


def _render_assay(result, output_format: str, provenance=None) -> str:
    pairs = result.get("pairs", [])
    if output_format == "json":
        output = dict(result)
        if provenance is not None:
            output["provenance"] = provenance
        return json.dumps(output, indent=2, default=str)
    if output_format == "csv":
        return OUT.pairs_to_csv(pairs)
    if output_format == "bed":
        return OUT.products_to_bed(pairs)
    if output_format == "order":
        return OUT.order_table(pairs)
    if output_format == "html":
        context = {
            "title": "Assay: %s" % result["target"]["name"],
            "template": "%s:%s-%s" % (
                result["target"]["chrom"], result["target"]["start"],
                result["target"]["end"]),
            "databases": result["databases"],
            "generated": (provenance or {}).get("generated", ""),
            "params": {},
            "provenance": provenance or {},
            "pairs": pairs,
        }
        return OUT.html_report(context)
    return R.assay_to_text(result)


def _cmd_assay(arguments) -> int:
    from .genome import Genome
    from .regions import resolve_gene, resolve_interval, resolve_snp
    from .assay import run_assay
    from .provenance import make_manifest

    design_genome = Genome(arguments.genome)
    caps_snp = None
    if arguments.gene:
        region = resolve_gene(
            arguments.gff3, arguments.gene,
            feature=arguments.gene_feature, flank=0)
    elif arguments.interval:
        chromosome, span = arguments.interval.split(":")
        start, end = span.split("-")
        region = resolve_interval(
            chromosome, int(start), int(end),
            name=arguments.name or arguments.interval)
    elif arguments.snp:
        chromosome, position = arguments.snp.split(":")
        region = resolve_snp(
            chromosome, int(position), flank=arguments.flank or 250,
            name=arguments.name)
        if arguments.alt:
            caps_snp = {
                "genomic_pos": int(position), "alt": arguments.alt}
    else:
        raise SystemExit(
            "assay: give --gene (with --gff3), --interval, or --snp")

    design = DesignParams(
        product_size_ranges=_parse_size_ranges(arguments.product_size),
        opt_tm=arguments.opt_tm,
        min_tm=arguments.min_tm,
        max_tm=arguments.max_tm,
        min_gc=arguments.min_gc,
        max_gc=arguments.max_gc,
        num_return=arguments.num_return,
    )
    specificity = _spec_from_args(arguments)
    variants = _load_variants(arguments.vcf)
    genomes_by_db, thermo_params, thermo_gate = _thermo_setup(
        arguments, design_genome=design_genome)
    dimer_params = _dimer_params_from_args(arguments)

    result = run_assay(
        region,
        design_genome,
        arguments.db,
        flank=arguments.flank,
        design_params=design,
        spec_params=specificity,
        variants=variants,
        caps_snp=caps_snp,
        primer3_bin=arguments.primer3_bin,
        blastn_bin=arguments.blastn_bin,
        genomes_by_db=genomes_by_db,
        thermo_params=thermo_params,
        thermo_gate=thermo_gate,
        dimer_params=dimer_params,
    )
    manifest = make_manifest(
        {
            "design": design.__dict__,
            "spec": specificity.__dict__,
            "flank": arguments.flank,
            "thermo_genomes": _genome_paths(genomes_by_db),
        },
        arguments.db,
        template_info=result["target"],
    )
    _emit(_render_assay(result, arguments.format, manifest), arguments.out)
    return 0


def _cmd_markers(arguments) -> int:
    from .genome import Genome
    from .regions import resolve_interval
    from .assay import design_qtl_markers

    design_genome = Genome(arguments.genome)
    chromosome, span = arguments.interval.split(":")
    start, end = span.split("-")
    qtl = resolve_interval(
        chromosome, int(start), int(end), name=arguments.name or "QTL")
    design = DesignParams(
        product_size_ranges=_parse_size_ranges(arguments.product_size),
        opt_tm=arguments.opt_tm,
        min_tm=arguments.min_tm,
        max_tm=arguments.max_tm,
        num_return=arguments.num_return,
    )
    specificity = _spec_from_args(arguments)
    genomes_by_db, thermo_params, thermo_gate = _thermo_setup(
        arguments, design_genome=design_genome)
    dimer_params = _dimer_params_from_args(arguments)
    markers = design_qtl_markers(
        qtl,
        design_genome,
        arguments.db,
        n_markers=arguments.n_markers,
        spacing=arguments.spacing,
        marker_flank=arguments.marker_flank,
        design_params=design,
        spec_params=specificity,
        primer3_bin=arguments.primer3_bin,
        blastn_bin=arguments.blastn_bin,
        genomes_by_db=genomes_by_db,
        thermo_params=thermo_params,
        thermo_gate=thermo_gate,
        dimer_params=dimer_params,
    )
    if arguments.format == "json":
        _emit(json.dumps({
            "interval": arguments.interval,
            "thermo_genomes": _genome_paths(genomes_by_db),
            "markers": markers,
        }, indent=2, default=str), arguments.out)
    else:
        lines = ["QTL %s: %s markers" % (arguments.interval, len(markers))]
        for marker in markers:
            if marker.get("pairs"):
                pair = marker["pairs"][0]
                lines.append(
                    "  %s @ %s: risk %s  %sbp  F %s / R %s" % (
                        marker["marker"], marker.get("anchor"), pair["risk"],
                        pair["product_size"], pair["forward"], pair["reverse"]))
            else:
                lines.append("  %s @ %s: %s" % (
                    marker["marker"], marker.get("anchor"),
                    marker.get("error", "no pair found")))
        _emit("\n".join(lines), arguments.out)
    return 0


def _cmd_makedb(arguments) -> int:
    output = make_blastdb(
        arguments.infile,
        out=arguments.out_db,
        title=arguments.title,
        parse_seqids=not arguments.no_parse_seqids,
        makeblastdb_bin=arguments.makeblastdb_bin,
    )
    print("built database: %s" % output)
    return 0


def _add_template_args(parser) -> None:
    template = parser.add_argument_group("template")
    source = template.add_mutually_exclusive_group(required=True)
    source.add_argument("--template", help="template DNA sequence")
    source.add_argument("--template-fasta", help="FASTA with one or more templates")
    template.add_argument("--template-id", default="template")


def _add_design_knobs(parser) -> None:
    design = parser.add_argument_group("primer3")
    design.add_argument("--opt-size", type=int, default=20)
    design.add_argument("--min-size", type=int, default=18)
    design.add_argument("--max-size", type=int, default=25)
    design.add_argument("--opt-tm", type=float, default=60.0)
    design.add_argument("--min-tm", type=float, default=57.0)
    design.add_argument("--max-tm", type=float, default=63.0)
    design.add_argument("--min-gc", type=float, default=20.0)
    design.add_argument("--max-gc", type=float, default=80.0)
    design.add_argument("--size-tolerance", type=int, default=10)
    design.add_argument("--primer3-bin")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="primerblast-oss",
        description=(
            "Local, offline Primer-BLAST for plant breeding with multi-reference "
            "specificity, per-assembly thermodynamics, tiling, and assays."),
        epilog="Run 'primerblast-oss <subcommand> --help' for options.",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__)
    subcommands = parser.add_subparsers(dest="cmd", required=True)

    design = subcommands.add_parser(
        "design", help="design primer pairs and check specificity")
    _add_template_args(design)
    design.add_argument("--product-size", default="70-1000")
    design.add_argument("--num-return", type=int, default=10)
    design.add_argument("--target", help="start,length on template (0-based)")
    _add_design_knobs(design)
    _add_spec_args(design)
    _add_dimer_args(design)
    _add_out_args(design, formats=("text", "json", "tsv"))
    design.set_defaults(func=_cmd_design)

    check = subcommands.add_parser(
        "check", help="in-silico PCR for supplied primers")
    check.add_argument("--forward")
    check.add_argument("--reverse")
    check.add_argument("--primer", action="append", help="NAME=SEQ or SEQ")
    check.add_argument("--primers-fasta")
    _add_spec_args(check)
    _add_out_args(check, formats=("text", "json"))
    check.set_defaults(func=_cmd_check)

    multiplex = subcommands.add_parser(
        "multiplex", help="check primer-dimer compatibility across a primer pool")
    multiplex.add_argument("--forward")
    multiplex.add_argument("--reverse")
    multiplex.add_argument("--primer", action="append", help="NAME=SEQ or SEQ")
    multiplex.add_argument("--primers-fasta")
    _add_dimer_args(multiplex)
    _add_out_args(multiplex, formats=("text", "json"))
    multiplex.set_defaults(func=_cmd_multiplex)

    multiplex_design = subcommands.add_parser(
        "multiplex-design", help="design one compatible pair per target")
    _add_template_args(multiplex_design)
    multiplex_design.add_argument("--product-size", default="80-300")
    multiplex_design.add_argument("--num-return", type=int, default=10)
    multiplex_design.add_argument("--candidates-per-target", type=int, default=5)
    multiplex_design.add_argument("--require-specific", action="store_true")
    _add_design_knobs(multiplex_design)
    _add_spec_args(multiplex_design)
    _add_dimer_args(multiplex_design)
    _add_out_args(multiplex_design, formats=("text", "json"))
    multiplex_design.set_defaults(func=_cmd_multiplex_design)

    tile = subcommands.add_parser(
        "tile", help="cover a region with overlapping amplicons")
    _add_template_args(tile)
    tile.add_argument("--region", help="start,end (0-based inclusive)")
    tile.add_argument("--amplicon-min", type=int, default=400)
    tile.add_argument("--amplicon-max", type=int, default=800)
    tile.add_argument("--overlap", type=int, default=40)
    tile.add_argument("--candidates-per-tile", type=int, default=8)
    _add_design_knobs(tile)
    _add_spec_args(tile)
    _add_dimer_args(tile)
    _add_out_args(tile, formats=("text", "json"))
    tile.set_defaults(func=_cmd_tile)

    assay = subcommands.add_parser(
        "assay", help="design + specificity + variants + CAPS + risk")
    target = assay.add_argument_group("target (choose one)")
    target.add_argument("--gene", help="gene id (requires --gff3)")
    target.add_argument(
        "--gene-feature", default="cds",
        choices=["gene", "mrna", "exon", "cds"])
    target.add_argument("--interval", help="chrom:start-end (1-based)")
    target.add_argument("--snp", help="chrom:pos")
    target.add_argument("--alt", help="alternate allele at --snp")
    target.add_argument("--name")
    assay.add_argument(
        "--genome", required=True,
        help="design-reference FASTA (.fai indexed); associated with first --db")
    assay.add_argument("--gff3")
    assay.add_argument("--vcf")
    assay.add_argument("--flank", type=int, default=200)
    assay.add_argument("--product-size", default="100-800")
    assay.add_argument("--num-return", type=int, default=10)
    _add_design_knobs(assay)
    _add_spec_args(assay)
    _add_dimer_args(assay)
    _add_out_args(
        assay, formats=("text", "json", "csv", "bed", "order", "html"))
    assay.set_defaults(func=_cmd_assay)

    markers = subcommands.add_parser(
        "markers", help="design evenly spaced markers across a QTL interval")
    markers.add_argument("--interval", required=True)
    markers.add_argument(
        "--genome", required=True,
        help="design-reference FASTA associated with first --db")
    markers.add_argument("--name")
    markers.add_argument("--n-markers", type=int, default=0)
    markers.add_argument("--spacing", type=int, default=0)
    markers.add_argument("--marker-flank", type=int, default=300)
    markers.add_argument("--product-size", default="100-500")
    markers.add_argument("--num-return", type=int, default=6)
    _add_design_knobs(markers)
    _add_spec_args(markers)
    _add_dimer_args(markers)
    _add_out_args(markers, formats=("text", "json"))
    markers.set_defaults(func=_cmd_markers)

    makedb = subcommands.add_parser(
        "makedb", help="build a BLAST nucleotide database from FASTA")
    makedb.add_argument("infile")
    makedb.add_argument("--out-db")
    makedb.add_argument("--title")
    makedb.add_argument("--no-parse-seqids", action="store_true")
    makedb.add_argument("--makeblastdb-bin")
    makedb.set_defaults(func=_cmd_makedb)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    arguments = parser.parse_args(argv)
    try:
        return arguments.func(arguments)
    except (ValueError, KeyError, RuntimeError) as error:
        parser.error(str(error))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
