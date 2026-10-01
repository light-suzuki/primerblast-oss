"""Shared experiment application layer for CLI, GUI and local agents."""
from __future__ import annotations
import json
import tempfile
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple
from .design import DesignParams, clean_sequence, read_fasta
from .specificity import SpecParams, in_silico_pcr
from .pipeline import run_pipeline
from .tiling import design_tiling
from .tools import make_blastdb
from . import report as R
from . import outputs as OUT

def _f(params: Dict, key: str, default):
    """Fetch a value, treating '' / None as 'use default'."""
    val = params.get(key, default)
    if val is None or val == "":
        return default
    return val


def _parse_size_ranges(text: str) -> List[Tuple[int, int]]:
    ranges = []
    for chunk in str(text).replace(",", " ").split():
        lo, hi = chunk.split("-")
        ranges.append((int(lo), int(hi)))
    return ranges or [(70, 1000)]


def _design_params(p: Dict, *, want_ranges: bool = True) -> DesignParams:
    kwargs = dict(
        opt_size=int(_f(p, "opt_size", 20)),
        min_size=int(_f(p, "min_size", 18)),
        max_size=int(_f(p, "max_size", 25)),
        opt_tm=float(_f(p, "opt_tm", 60.0)),
        min_tm=float(_f(p, "min_tm", 57.0)),
        max_tm=float(_f(p, "max_tm", 63.0)),
        min_gc=float(_f(p, "min_gc", 20.0)),
        max_gc=float(_f(p, "max_gc", 80.0)),
        num_return=int(_f(p, "num_return", 10)),
    )
    if want_ranges:
        kwargs["product_size_ranges"] = _parse_size_ranges(_f(p, "product_size", "70-1000"))
    target = _f(p, "target", None)
    if target:
        x, y = str(target).split(",")
        kwargs["target"] = (int(x), int(y))
    return DesignParams(**kwargs)


def _spec_params(p: Dict) -> SpecParams:
    return SpecParams(
        max_total_mismatch=int(_f(p, "max_total_mismatch", 4)),
        max_3prime_mismatch=int(_f(p, "max_3prime_mismatch", 1)),
        three_prime_window=int(_f(p, "three_prime_window", 5)),
        require_3prime_terminal_match=not bool(p.get("no_3prime_terminal", False)),
        min_product=int(_f(p, "min_product", 40)),
        max_product=int(_f(p, "max_product", 4000)),
        gel_min_gap_bp=int(_f(p, "gel_min_gap", 50)),
        word_size=int(_f(p, "word_size", 7)),
        num_threads=int(_f(p, "num_threads", 2)),
        max_target_seqs=int(_f(p, "max_target_seqs", 5000)),
    )


def _associated_genomes(p: Dict, design_genome=None) -> Dict:
    from .genome import Genome

    databases = _databases(p)
    mappings = p.get("db_genomes") or {}
    if isinstance(mappings, str):
        parsed = {}
        for line in mappings.splitlines():
            if line.strip():
                database, separator, fasta = line.partition("=")
                if not separator or not fasta.strip():
                    raise ValueError("Use DB=FASTA, one association per line.")
                parsed[database.strip()] = fasta.strip()
        mappings = parsed
    if not isinstance(mappings, dict) or any(db not in databases for db in mappings):
        raise ValueError("FASTA associations must name selected databases.")
    genomes = {database: Genome(str(fasta)) for database, fasta in mappings.items()}
    if design_genome is None and _f(p, "genome", None):
        design_genome = Genome(str(p["genome"]))
    if design_genome is not None:
        existing = genomes.get(databases[0])
        if existing is not None and Path(existing.fasta).resolve() != Path(design_genome.fasta).resolve():
            raise ValueError("The first database FASTA must match the design genome.")
        genomes[databases[0]] = design_genome
    return genomes


def _databases(p: Dict) -> List[str]:
    dbs = [d for d in (p.get("db") or []) if str(d).strip()]
    if not dbs:
        raise ValueError("At least one BLAST database is required.")
    return dbs


def _templates(p: Dict) -> List[Tuple[str, str]]:
    """Return [(id, seq), ...] from a raw sequence or pasted FASTA text."""
    text = str(p.get("template", "")).strip()
    if not text:
        raise ValueError("A template sequence (or FASTA) is required.")
    if text.lstrip().startswith(">"):
        records: List[Tuple[str, str]] = []
        name, buf = None, []
        for line in text.splitlines():
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(buf)))
                header = line[1:].strip()
                if not header:
                    raise ValueError("FASTA records require a nonempty identifier.")
                name = header.split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            records.append((name, "".join(buf)))
        return records or [("template", "")]
    return [(str(p.get("template_id", "template") or "template"), text)]


# --------------------------------------------------------------------------- #
# mode handlers -> JSON-serializable dict
# --------------------------------------------------------------------------- #
def _find_gene_seqid(gff3_path: str, gene: str) -> Optional[str]:
    """Bound the parse only when normalized gene matches share one seqid.

    Cross-chromosome aliases must reach the full parser's ambiguity check."""
    from .annotation_index import gene_annotation
    try:
        selected = gene_annotation(gff3_path, gene).gene(gene)
    except (OSError, ValueError):
        return None
    return selected.seqid if selected else None


def _gene_to_template(p: Dict):
    """Resolve a gene name/ID (+GFF3 +genome), retaining its genomic context."""
    from .genome import Genome
    from .regions import resolve_gene, extract_template

    gene = _f(p, "gene", None)
    gff3 = _f(p, "gff3", None)
    genome_path = _f(p, "genome", None)
    if not gene:
        raise ValueError("Enter a gene name or ID.")
    if not gff3:
        raise ValueError("A GFF3 annotation path is required for gene-based design.")
    if not genome_path:
        raise ValueError("A reference genome FASTA (.fai indexed) path is required.")
    genome = Genome(genome_path)
    seqid = _find_gene_seqid(gff3, gene)   # bound the parse to one chromosome
    region = resolve_gene(gff3, gene, feature=_f(p, "gene_feature", "cds"),
                          flank=0, gff3_seqid=seqid)
    tmpl = extract_template(genome, region, flank=int(_f(p, "flank", 0)))
    return tmpl


def _run_design(p: Dict) -> Dict:
    dp = _design_params(p)
    sp = _spec_params(p)
    dbs = _databases(p)
    size_tol = int(_f(p, "size_tolerance", 10))
    if _f(p, "source", "sequence") == "gene":
        genomic_template = _gene_to_template(p)
        templates = [(genomic_template.id, genomic_template.seq)]
    else:
        genomic_template = None
        templates = _templates(p)
    results = []
    for tid, seq in templates:
        res = run_pipeline(tid, seq, dbs, design_params=dp, spec_params=sp,
                           size_tolerance=size_tol, genomes_by_db=_associated_genomes(p))
        d = R.to_dict(res)
        d["template_sequence"] = clean_sequence(seq)
        if genomic_template is not None:
            from .annotations import template_annotations
            template = genomic_template
            d["template"] = {"sequence": template.seq, "chrom": template.region.chrom,
                             "start": template.ext_start, "end": template.ext_end,
                             "anchor": template.anchor_coord, "strand": template.anchor_strand}
            d["template"]["annotations"] = template_annotations(_f(p, "gff3", None), d["template"])
        d["tsv"] = R.to_tsv(res)
        results.append(d)
    return {"mode": "design", "templates": results}


def _run_check(p: Dict) -> Dict:
    from .primer_evidence import oligo_hypotheses, check_evidence_report
    from .sequence_tools import dna_input, product_sequence
    sp = _spec_params(p)
    dbs = _databases(p)
    primers: Dict[str, str] = {}
    if p.get("forward"):
        primers["F"] = dna_input(p["forward"])
    if p.get("reverse"):
        primers["R"] = dna_input(p["reverse"])
    for i, spec in enumerate(p.get("primers") or [], 1):
        spec = str(spec).strip()
        if not spec:
            continue
        if "=" in spec:
            name, seq = spec.split("=", 1)
        else:
            name, seq = f"P{i}", spec
        name = name.strip()
        if not name or name in primers:
            raise ValueError("Primer names must be nonempty and unique.")
        primers[name] = dna_input(seq)
    if not primers:
        raise ValueError("Provide at least a forward/reverse primer or a primer list.")
    genomes = _associated_genomes(p)
    orientation = _f(p, "input_orientation", "as_supplied")
    hypotheses = oligo_hypotheses(primers, orientation)
    annotation_map = p.get("db_gff3") or {}
    if not isinstance(annotation_map, dict) or any(db not in dbs for db in annotation_map):
        raise ValueError("db_gff3 must map selected databases to matching GFF3 paths.")
    annotation_map = dict(annotation_map)
    if _f(p, "gff3", None) or _f(p, "annotation_gff3", None):
        annotation_map.setdefault(dbs[0], _f(p, "gff3", None) or p["annotation_gff3"])
    results = []
    for oligos in hypotheses:
        raw = [in_silico_pcr(oligos, db, sp=sp, genome=genomes.get(db)) for db in dbs]
        reports = R.insilico_to_dict(raw, oligos)["results"]
        for result in reports:
            result["oligos"] = oligos
            result["reverse_complemented_inputs"] = [name for name in primers if primers[name] != oligos[name]]
            for amplicon in result["products"]:
                product_sequence(amplicon, result["db"], genomes.get(result["db"]), annotation_map.get(result["db"]))
            result["fasta"] = "".join(a.get("fasta", "") for a in result["products"])
        results.extend(reports)
    return check_evidence_report(results, primers, orientation)


def _run_blast(p: Dict) -> Dict:
    from .sequence_tools import nucleotide_search
    return nucleotide_search(_templates(p), _databases(p), task=_f(p, "task", "blastn"),
                             evalue=float(_f(p, "evalue", 10)),
                             max_target_seqs=int(_f(p, "max_target_seqs", 100)),
                             num_threads=int(_f(p, "num_threads", 2)))


def _run_primer3(p: Dict) -> Dict:
    from dataclasses import asdict
    from .design import design_primers
    from .sequence_tools import dna_input, fasta_record
    params = _design_params(p)
    templates = []
    for name, sequence in _templates(p):
        sequence = clean_sequence(dna_input(sequence))
        pairs, explain = (design_primers(name, sequence, params, primer3_bin=p["primer3_bin"])
                          if p.get("primer3_bin") else design_primers(name, sequence, params))
        rows = []
        for pair in pairs:
            row = asdict(pair)
            row["sequence"] = sequence[pair.left_start:pair.right_start + 1]
            row["fasta"] = fasta_record("%s_pair%d_reference" % (name, pair.index + 1), row["sequence"])
            rows.append(row)
        templates.append({"template_id": name, "template_sequence": sequence,
                          "template_len": len(sequence), "pairs": rows, "primer3_explain": explain})
    return {"mode": "primer3", "specificity_status": "not_evaluated", "templates": templates}


def _run_tile(p: Dict) -> Dict:
    dp = _design_params(p, want_ranges=False)
    sp = _spec_params(p)
    dbs = _databases(p)
    size_tol = int(_f(p, "size_tolerance", 10))
    region = None
    if _f(p, "region", None):
        x, y = str(p["region"]).split(",")
        region = (int(x), int(y))
    out = []
    for tid, seq in _templates(p):
        tiles = design_tiling(
            tid, seq, dbs, region=region,
            amplicon_min=int(_f(p, "amplicon_min", 400)),
            amplicon_max=int(_f(p, "amplicon_max", 800)),
            overlap=int(_f(p, "overlap", 40)),
            design_params=dp, spec_params=sp, size_tolerance=size_tol,
            candidates_per_tile=int(_f(p, "candidates_per_tile", 8)),
            genomes_by_db=_associated_genomes(p),
        )
        reg = region or (0, len(clean_sequence(seq)) - 1)
        out.append(R.tiling_to_dict(tiles, tid, reg, dbs))
    return {"mode": "tile", "templates": out}


def _run_assay(p: Dict) -> Dict:
    from .genome import Genome
    from .regions import resolve_gene, resolve_interval, resolve_snp
    from .assay import run_assay
    from .provenance import make_manifest
    from .vcf import parse_vcf

    genome_path = _f(p, "genome", None)
    if not genome_path:
        raise ValueError("A reference genome FASTA (.fai-indexed) path is required.")
    genome = Genome(genome_path)
    dbs = _databases(p)
    flank = int(_f(p, "flank", 200))
    caps_snp = None
    if _f(p, "gene", None):
        if not _f(p, "gff3", None):
            raise ValueError("--gene requires a GFF3 annotation path.")
        region = resolve_gene(p["gff3"], p["gene"],
                              feature=_f(p, "gene_feature", "cds"), flank=0,
                              gff3_seqid=_find_gene_seqid(p["gff3"], p["gene"]))
    elif _f(p, "interval", None):
        chrom, span = str(p["interval"]).split(":")
        s, e = span.split("-")
        region = resolve_interval(chrom, int(s), int(e),
                                  name=_f(p, "name", p["interval"]), strand=_f(p, "strand", "+"))
    elif _f(p, "snp", None):
        chrom, pos = str(p["snp"]).split(":")
        region = resolve_snp(chrom, int(pos), flank=flank or 250,
                             name=_f(p, "name", None))
        if _f(p, "alt", None):
            caps_snp = {"genomic_pos": int(pos), "alt": str(p["alt"]).upper()}
    else:
        raise ValueError("Choose a target: gene (+GFF3), interval, or SNP.")

    dp = _design_params(p)
    sp = _spec_params(p)
    variants = parse_vcf(p["vcf"]) if _f(p, "vcf", None) else []
    result = run_assay(region, genome, dbs, flank=flank, design_params=dp,
                       spec_params=sp, variants=variants, caps_snp=caps_snp,
                       genomes_by_db=_associated_genomes(p, genome),
                       gel_ladder=_f(p, "gel_ladder", "auto"),
                       custom_ladder_bands=[int(v) for v in str(_f(p, "ladder_bands", "")).split(",") if v.strip()],
                       gel_percent=float(p["gel_percent"]) if _f(p, "gel_percent", "auto") != "auto" else None,
                       aspcr_pairs_to_screen=int(_f(p, "aspcr_pairs_to_screen", 2)))
    prov = make_manifest({"design": dp.__dict__, "spec": sp.__dict__, "flank": flank},
                         dbs, template_info=result["target"])
    result["provenance"] = prov
    from .annotations import template_annotations
    result["template"]["annotations"] = template_annotations(_f(p, "annotation_gff3", None) or _f(p, "gff3", None), result["template"])
    pairs = result.get("pairs", [])
    result["exports"] = {
        "csv": _safe(OUT.pairs_to_csv, pairs),
        "order": _safe(OUT.order_table, pairs),
        "bed": _safe(OUT.products_to_bed, pairs),
        "ascii": "\n\n".join(_safe(OUT.ascii_offtarget_map, pr) or "" for pr in pairs),
    }
    result["mode"] = "assay"
    from .gel import best_analysis_from_assay, virtual_gel_svg
    analysis = best_analysis_from_assay(result)
    if analysis:
        result["exports"]["svg"] = virtual_gel_svg(analysis)
    return result


def _run_sequence(p: Dict) -> Dict:
    """Reuse the CLI orchestration so GUI and CLI share coverage semantics."""
    from .cli import build_parser, _cmd_sequence

    with tempfile.TemporaryDirectory(prefix="primerblast-sequence-") as folder:
        root = Path(folder)
        output = root / "result.json"
        args = ["sequence", "--format", "json", "--out", str(output)]
        source = _f(p, "source", "sequence")
        if source in ("gene", "interval"):
            args += ["--" + source, str(_f(p, source, "")), "--genome", str(_f(p, "genome", ""))]
        else:
            template = str(_f(p, "template", ""))
            if template.lstrip().startswith(">"):
                fasta = root / "template.fa"
                fasta.write_text(template, encoding="utf-8")
                args += ["--template-fasta", str(fasta)]
            else:
                args += ["--template", template]
        for database in _databases(p):
            args += ["--db", database]
        mapping = p.get("db_genomes") or {}
        if isinstance(mapping, str):
            entries = [line.strip() for line in mapping.splitlines() if line.strip()]
        else:
            entries = [str(db) + "=" + str(fasta) for db, fasta in mapping.items()]
        # Validate mappings with the same contract as the other GUI modes.
        _associated_genomes(p)
        for entry in entries:
            args += ["--db-genome", entry]
        keys = ("template_id", "gff3", "gene_feature", "strand", "flank", "overlap",
                "amplicon_size", "candidates_per_tile", "size_tolerance", "opt_size",
                "min_size", "max_size", "opt_tm", "min_tm", "max_tm", "min_gc", "max_gc",
                "max_total_mismatch", "max_3prime_mismatch", "three_prime_window",
                "min_product", "max_product", "gel_min_gap", "word_size", "num_threads", "max_target_seqs")
        for key in keys:
            if _f(p, key, None) is not None:
                args += ["--" + key.replace("_", "-"), str(p[key])]
        for flag in ("m13_tails", "no_3prime_terminal"):
            if p.get(flag):
                args.append("--" + flag.replace("_", "-"))
        try:
            parsed = build_parser().parse_args(args)
        except SystemExit as exc:
            raise ValueError("Invalid sequencing parameters; check size ranges and numeric values") from exc
        _cmd_sequence(parsed)
        result = json.loads(output.read_text(encoding="utf-8"))
        from .sequencing import order_table
        plans = result.get("plans", [result])
        from .annotations import template_annotations
        for plan in plans:
            if plan.get("template"):
                plan["template"]["annotations"] = template_annotations(_f(p, "annotation_gff3", None) or _f(p, "gff3", None), plan["template"])
        result["exports"] = {"order": "\n".join(order_table(plan) for plan in plans)}
        return result


def _run_markers(p: Dict) -> Dict:
    from .genome import Genome
    from .regions import resolve_interval
    from .assay import design_qtl_markers

    genome_path = _f(p, "genome", None)
    if not genome_path:
        raise ValueError("A reference genome FASTA path is required.")
    genome = Genome(genome_path)
    dbs = _databases(p)
    chrom, span = str(_f(p, "interval", "")).split(":")
    s, e = span.split("-")
    qtl = resolve_interval(chrom, int(s), int(e), name=_f(p, "name", "QTL"))
    dp = _design_params(p)
    sp = _spec_params(p)
    markers = design_qtl_markers(
        qtl, genome, dbs,
        n_markers=int(_f(p, "n_markers", 0)),
        spacing=int(_f(p, "spacing", 0)),
        marker_flank=int(_f(p, "marker_flank", 300)),
        design_params=dp, spec_params=sp,
        genomes_by_db=_associated_genomes(p, genome),
    )
    return {"mode": "markers", "interval": p.get("interval"), "markers": markers}


def _run_makedb(p: Dict) -> Dict:
    infile = _f(p, "infile", None)
    if not infile:
        raise ValueError("An input FASTA path is required.")
    out = make_blastdb(infile, out=_f(p, "out_db", None), title=_f(p, "title", None),
                       parse_seqids=not bool(p.get("no_parse_seqids", False)))
    return {"mode": "makedb", "db": out,
            "message": f"Built BLAST database: {out}"}


def _safe(fn: Callable, *args):
    try:
        return fn(*args)
    except Exception:  # exports are best-effort; never fail the whole job
        return None


HANDLERS: Dict[str, Callable[[Dict], Dict]] = {
    "blast": _run_blast,
    "primer3": _run_primer3,
    "design": _run_design,
    "check": _run_check,
    "tile": _run_tile,
    "sequence": _run_sequence,
    "assay": _run_assay,
    "markers": _run_markers,
    "makedb": _run_makedb,
}


def execute(operation, params, *, allow_db_write=False):
    """Single application boundary; no HTTP, browser or model provider required."""
    if operation not in HANDLERS:
        raise ValueError("Unknown operation: " + str(operation))
    if not isinstance(params, dict):
        raise ValueError("params must be a JSON object")
    if operation == "makedb" and not allow_db_write:
        raise ValueError("Database creation requires explicit --allow-db-write")
    return HANDLERS[operation](params)
