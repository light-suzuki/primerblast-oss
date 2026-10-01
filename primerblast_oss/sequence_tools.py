"""Local nucleotide search and sequence exports (1-based inclusive coordinates)."""
import math
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from .annotations import template_annotations


def dna_input(text):
    sequence = re.sub(r"\s+", "", str(text)).upper()
    if not sequence or re.search(r"[^ACGTRYSWKMBDHVN]", sequence):
        raise ValueError("Enter DNA using ACGT or IUPAC bases only.")
    return sequence


def fasta_record(name, sequence):
    name = re.sub(r"\s+", "_", str(name))
    return ">" + name + "\n" + "\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n"


def nucleotide_search(records, databases, *, task="blastn", evalue=10.0,
                      max_target_seqs=100, num_threads=2):
    """Ordinary BLAST HSPs; do not apply PCR priming/3'-end filters."""
    if task not in ("blastn", "megablast", "blastn-short"):
        raise ValueError("Choose blastn, megablast or blastn-short.")
    if not math.isfinite(evalue) or evalue <= 0 or not 1 <= max_target_seqs <= 10000 or not 1 <= num_threads <= 128:
        raise ValueError("Invalid BLAST search limits.")
    exe = shutil.which("blastn")
    if not exe:
        raise RuntimeError("blastn not found. Install BLAST+.")
    columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
    queries = [{"id": "query%d" % i, "name": name, "sequence": dna_input(seq)}
               for i, (name, seq) in enumerate(records, 1)]
    results = []
    with tempfile.TemporaryDirectory(prefix="primerblast-search-") as folder:
        query = Path(folder) / "query.fa"
        query.write_text("".join(fasta_record(q["id"], q["sequence"]) for q in queries), encoding="utf-8")
        for database in databases:
            command = [exe, "-query", str(query), "-db", database, "-task", task,
                       "-strand", "both", "-evalue", str(evalue), "-max_target_seqs", str(max_target_seqs),
                       "-num_threads", str(num_threads), "-outfmt", "6 " + columns]
            process = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if process.returncode:
                raise RuntimeError("BLAST failed: " + process.stderr.decode(errors="replace").strip())
            hits = []
            for line in process.stdout.decode().splitlines():
                row = dict(zip(columns.split(), line.split("\t")))
                for key in ("length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send"):
                    row[key] = int(row[key])
                for key in ("pident", "evalue", "bitscore"):
                    row[key] = float(row[key])
                row["strand"] = "+" if row["sstart"] <= row["send"] else "-"
                hits.append(row)
            limited = [q["id"] for q in queries if len({h["sseqid"] for h in hits if h["qseqid"] == q["id"]}) >= max_target_seqs]
            results.append({"db": database, "hits": hits, "at_target_limit": limited,
                            "search_scope": "reported_hsps", "stderr": process.stderr.decode(errors="replace").strip(),
                            "tsv": "\t".join(columns.split()) + "\n" + process.stdout.decode()})
    return {"mode": "blast", "task": task, "queries": queries, "results": results,
            "parameters": {"strand": "both", "evalue": evalue, "max_target_seqs": max_target_seqs, "num_threads": num_threads}}


def product_sequence(product, database, genome=None, annotation=None):
    """Export reference bases; no invented gene names or substituted oligo bases."""
    sequence = None
    source = "genome_fasta" if genome else "blast_database"
    try:
        if genome:
            sequence = genome.fetch(product["subject"], product["start"], product["end"])
        else:
            exe = shutil.which("blastdbcmd")
            if not exe:
                raise RuntimeError("Supply the matching genome FASTA (.fai) or install blastdbcmd.")
            proc = subprocess.run([exe, "-db", database, "-entry", product["subject"],
                                   "-range", "%s-%s" % (product["start"], product["end"]),
                                   "-strand", "plus", "-outfmt", "%s"],
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if proc.returncode:
                raise RuntimeError(proc.stderr.decode(errors="replace").strip())
            sequence = dna_input(proc.stdout.decode())
        if len(sequence) != product["size"]:
            raise ValueError("Extracted sequence length does not match product coordinates.")
        product.update(sequence=sequence, sequence_status="available", sequence_source=source,
                       fasta=fasta_record("%s:%s-%s_reference_plus" % (product["subject"], product["start"], product["end"]), sequence))
    except (OSError, RuntimeError, ValueError, KeyError) as error:
        product.update(sequence_status="unavailable", sequence_error=str(error))
    context = {"chrom": product["subject"], "start": product["start"], "end": product["end"],
               "anchor": product["start"], "strand": "+", "sequence": sequence or ""}
    try:
        product["annotations"] = template_annotations(annotation, context)
    except (OSError, ValueError) as error:
        product["annotations"] = {"status": "unavailable", "genes": [], "error": str(error)}
