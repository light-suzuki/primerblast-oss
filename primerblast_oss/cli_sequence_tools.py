"""File/stdin inputs and tabular exports for standalone sequence commands."""
import csv
import io
import json
import sys
from pathlib import Path

from .sequence_tools import dna_input, fasta_record, nucleotide_search
from .workflows import _templates, execute


def distinct_outputs(*paths):
    resolved = [Path(path).resolve() for path in paths if path]
    if len(resolved) != len(set(resolved)):
        raise ValueError("Report and FASTA output paths must be different.")


def sequence_params(sequence, fasta, name):
    text = sequence if fasta is None else (
        sys.stdin.read() if fasta == "-" else Path(fasta).read_text(encoding="utf-8-sig"))
    text = text.lstrip("\ufeff").strip()
    if fasta is not None and not text.startswith(">"):
        raise ValueError("FASTA input must contain records beginning with >")
    params = {"template": text, "template_id": name}
    records = _templates(params)
    # Validate the entire batch before running an external tool.
    for _, seq in records:
        dna_input(seq)
    return params, records


def tsv_rows(header, rows):
    stream = io.StringIO()
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return stream.getvalue().rstrip("\n")


def blast_tsv(data):
    names = {q["id"]: q["name"] for q in data["queries"]}
    header, rows = None, []
    for result in data["results"]:
        lines = list(csv.reader(io.StringIO(result["tsv"]), delimiter="\t"))
        header = ["db", "query_name"] + lines[0] + ["strand"]
        rows.extend([result["db"], names[line[0]]] + line + ["+" if int(line[8]) <= int(line[9]) else "-"]
                    for line in lines[1:] if line)
    return tsv_rows(header, rows)


def blast_text(data):
    names = {q["id"]: q["name"] for q in data["queries"]}
    lines = ["BLAST %s (both strands; HSP alignments, without PCR filters)" % data["task"]]
    for result in data["results"]:
        lines.append("\n# %s: %d HSP(s)" % (result["db"], len(result["hits"])))
        for hit in result["hits"]:
            lines.extend([
                "  %s:%s-%s -> %s:%s-%s (%s) identity=%s%% length=%s E=%s" % (
                    names[hit["qseqid"]], hit["qstart"], hit["qend"], hit["sseqid"],
                    hit["sstart"], hit["send"], hit["strand"], hit["pident"], hit["length"], hit["evalue"]),
                "    Query   " + hit["qseq"], "    Subject " + hit["sseq"],
            ])
    return "\n".join(lines)


def cmd_blast(args):
    from .cli import _emit
    _, records = sequence_params(args.query, args.query_fasta, args.query_id)
    data = nucleotide_search(records, args.db, task=args.task, evalue=args.evalue,
                             max_target_seqs=args.max_target_seqs,
                             num_threads=args.num_threads, blastn_bin=args.blastn_bin)
    for result in data["results"]:
        if result["stderr"]:
            print(result["stderr"], file=sys.stderr)
        if result["at_target_limit"]:
            print("warning: BLAST target limit reached in %s for %s; search may be incomplete" % (
                result["db"], ", ".join(result["at_target_limit"])), file=sys.stderr)
    output = (json.dumps(data, indent=2) if args.format == "json" else
              blast_tsv(data) if args.format == "tsv" else blast_text(data))
    _emit(output, args.out)
    return 0


def primer3_tsv(data):
    return tsv_rows(
        ["template_id", "pair", "forward_5to3", "reverse_5to3", "product_size",
         "left_start_0based", "right_end_0based", "tm_f", "tm_r", "gc_f", "gc_r", "penalty",
         "specificity_status"],
        ([template["template_id"], p["index"] + 1, p["forward"], p["reverse"], p["product_size"],
          p["left_start"], p["right_start"], p["tm_f"], p["tm_r"], p["gc_f"], p["gc_r"],
          p["penalty"], "not_evaluated"] for template in data["templates"] for p in template["pairs"]))


def primer3_text(data):
    lines = ["Primer3 only: specificity not evaluated"]
    for template in data["templates"]:
        lines.append("\n# %s: %d pair(s)" % (template["template_id"], len(template["pairs"])))
        for p in template["pairs"]:
            lines.extend(["  pair %s: %s bp  Tm %s / %s" % (
                p["index"] + 1, p["product_size"], p["tm_f"], p["tm_r"]),
                "    F 5'-" + p["forward"] + "-3'", "    R 5'-" + p["reverse"] + "-3'"])
        lines.append("  Primer3 explain: " + template["primer3_explain"])
    return "\n".join(lines)


def cmd_primer3(args):
    from .cli import _emit, _parse_size_ranges
    distinct_outputs(args.out, args.primers_out, args.products_fasta)
    params, _ = sequence_params(args.template, args.template_fasta, args.template_id)
    _parse_size_ranges(args.product_size)
    params.update({key: getattr(args, key) for key in (
        "product_size", "num_return", "target", "opt_size", "min_size", "max_size",
        "opt_tm", "min_tm", "max_tm", "min_gc", "max_gc", "primer3_bin")})
    data = execute("primer3", params)
    for template in data["templates"]:
        if not template["pairs"]:
            print("%s: no primer pairs; %s" % (template["template_id"], template["primer3_explain"]),
                  file=sys.stderr)
    if args.primers_out:
        records = [fasta_record("template%d_pair%d_%s" % (i, p["index"] + 1, side), p[key])
                   for i, template in enumerate(data["templates"], 1) for p in template["pairs"]
                   for side, key in (("F", "forward"), ("R", "reverse"))]
        Path(args.primers_out).write_text("".join(records), encoding="utf-8")
    if args.products_fasta:
        records = [fasta_record("template%d_pair%d_reference" % (i, p["index"] + 1), p["sequence"])
                   for i, template in enumerate(data["templates"], 1) for p in template["pairs"]]
        Path(args.products_fasta).write_text("".join(records), encoding="utf-8")
    output = (json.dumps(data, indent=2) if args.format == "json" else
              primer3_tsv(data) if args.format == "tsv" else primer3_text(data))
    _emit(output, args.out)
    return 0


def add_parsers(subcommands):
    from .cli import _add_design_knobs, _add_template_args, _add_out_args
    blast = subcommands.add_parser("blast", help="ordinary nucleotide BLAST (no PCR filters)")
    source = blast.add_mutually_exclusive_group(required=True)
    source.add_argument("--query", help="DNA or pasted FASTA")
    source.add_argument("--query-fasta", help="FASTA file; - reads stdin")
    blast.add_argument("--query-id", default="query")
    blast.add_argument("--db", action="append", required=True, help="BLAST DB prefix; repeat per reference")
    blast.add_argument("--task", choices=("blastn", "megablast", "blastn-short"), default="blastn")
    blast.add_argument("--evalue", type=float, default=10.0)
    blast.add_argument("--max-target-seqs", type=int, default=100)
    blast.add_argument("--num-threads", type=int, default=2)
    blast.add_argument("--blastn-bin")
    _add_out_args(blast)
    blast.set_defaults(func=cmd_blast, format="tsv")

    primer3 = subcommands.add_parser("primer3", help="design primers without a specificity database")
    _add_template_args(primer3, stdin=True)
    primer3.add_argument("--product-size", default="70-1000")
    primer3.add_argument("--num-return", type=int, default=10)
    primer3.add_argument("--target", help="start,length on template (0-based)")
    _add_design_knobs(primer3, size_tolerance=False)
    primer3.add_argument("--primers-out", help="write candidate oligos as FASTA; all pairs form a pool")
    primer3.add_argument("--products-fasta", help="write designed reference product spans as FASTA")
    _add_out_args(primer3)
    primer3.set_defaults(func=cmd_primer3, format="tsv")


def hypothesis_numbers(data):
    keys = dict.fromkeys(tuple(sorted(r["oligos"].items())) for r in data["results"])
    return {key: i for i, key in enumerate(keys, 1)}


def check_tsv(data):
    numbers = hypothesis_numbers(data)
    return tsv_rows(
        ["db", "hypothesis", "changed_inputs", "search_completeness", "subject", "start_1based",
         "end_1based", "size", "left_primer", "right_primer", "left_oligo_5to3", "right_oligo_5to3",
         "left_direction", "right_direction", "left_mismatches", "right_mismatches", "sequence_status",
         "label_orientation", "original_search_complete", "as_supplied_locus_observed", "annotation_status", "genes"],
        ([r["db"], numbers[tuple(sorted(r["oligos"].items()))], ",".join(r["reverse_complemented_inputs"]), r["search_completeness"],
          p["subject"], p["start"], p["end"], p["size"], p["fwd_primer"], p["rev_primer"],
          r["oligos"][p["fwd_primer"]], r["oligos"][p["rev_primer"]], "+:right", "-:left",
          p["fwd_mismatch"], p["rev_mismatch"], p["input_evidence"]["sequence_status"],
          p["input_evidence"]["label_orientation"], p["input_evidence"]["original_search_complete"],
          p["input_evidence"]["as_supplied_locus_observed"], p.get("annotations", {}).get("status", "not_provided"),
          ",".join(g.get("name") or g.get("id", "") for g in p.get("annotations", {}).get("genes", []))]
         for r in data["results"] for p in r["products"]))


def check_text(data, raw_groups, show_forms=False):
    from . import report
    lines = ["Input-error diagnosis: undetermined; predictions use the selected search conditions."]
    for summary in data["input_assessments"]:
        lines.append("# %s: as-supplied two-primer products=%s; single-primer products=%s; "
                     "reverse-complement candidates=%s; original search complete=%s" % (
                         summary["db"], summary["distinct_primer_products"], summary["same_primer_products"],
                         summary["reverse_complement_candidates"], summary["original_search_complete"]))
    if show_forms:
        for name, forms in data["input_sequence_forms"].items():
            lines.extend("%s %s: %s" % (name, label, forms[key]) for key, label in (
                ("input_5to3", "input 5'-3'"), ("reverse", "reverse character order"),
                ("complement_3to5", "complement 3'-5'"), ("reverse_complement_5to3", "reverse complement 5'-3'")))
    for rows, oligos, changed in raw_groups:
        lines.append("\n## " + ("AS SUPPLIED" if not changed else "CHANGED OLIGOS (reverse complement): " + ", ".join(changed)))
        lines.append(report.insilico_to_text(rows, oligos))
    lines.append("\n## Per-product input evidence (+:right / -:left; 1-based coordinates)")
    for r in data["results"]:
        for p in r["products"]:
            e = p["input_evidence"]
            lines.append("  %s %s:%s-%s  %s + -> / %s - <-  %s; labels=%s; changed=%s" % (
                r["db"], p["subject"], p["start"], p["end"], p["fwd_primer"], p["rev_primer"],
                e["sequence_status"], e["label_orientation"], ",".join(e["changed_primers"]) or "none"))
            if p.get("annotations"):
                lines.append("    annotations=%s; genes=%s" % (p["annotations"]["status"],
                    ",".join(g.get("name") or g.get("id", "") for g in p["annotations"].get("genes", [])) or "-"))
    return "\n".join(lines)


def export_check_fasta(data, path):
    records, missing = [], 0
    db_numbers = {db: i for i, db in enumerate(dict.fromkeys(r["db"] for r in data["results"]), 1)}
    numbers = hypothesis_numbers(data)
    for r in data["results"]:
        hypothesis = numbers[tuple(sorted(r["oligos"].items()))]
        for i, p in enumerate(r["products"], 1):
            if not p.get("sequence"):
                missing += 1
                continue
            status = p["input_evidence"]["sequence_status"]
            name = "db%d_hypothesis%d_product%d_%s_%s:%s-%s_reference_plus" % (
                db_numbers[r["db"]], hypothesis, i, status, p["subject"], p["start"], p["end"])
            records.append(fasta_record(name, p["sequence"]))
    Path(path).write_text("".join(records), encoding="utf-8")
    if missing:
        print("warning: reference FASTA unavailable for %s predicted product(s); "
              "export is partial. Supply matching indexed FASTA or blastdbcmd." % missing, file=sys.stderr)
