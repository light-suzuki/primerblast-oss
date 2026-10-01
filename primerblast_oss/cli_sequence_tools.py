"""File/stdin primer inputs and PCR-check result exports."""
import csv
import io
import sys
from pathlib import Path

from .sequence_tools import dna_input, fasta_record
from .workflows import _templates


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
