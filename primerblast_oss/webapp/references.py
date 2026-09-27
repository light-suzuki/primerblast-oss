"""Small reference catalog: explicit profiles or unambiguous sibling files."""
import json
import os
import re
import gzip
from pathlib import Path


def infer_gene_id_format(genome, gff3):
    """Inspect a bounded annotation sample; never return gene sequences."""
    if not gff3 or not Path(gff3).is_file():
        return None
    from ..gff3 import _parse_attributes
    formats, chromosomes = set(), set()
    opener = gzip.open if str(gff3).endswith(".gz") else open
    try:
        with opener(gff3, "rt", encoding="utf-8", errors="replace") as handle:
            consumed, samples = 0, 0
            while consumed < 1024 * 1024 and samples < 32:
                line = handle.readline(16384)
                if not line:
                    break
                consumed += len(line)
                if len(line) == 16384 and not line.endswith("\n"):
                    continue
                fields = line.rstrip().split("\t")
                if len(fields) != 9 or fields[2].lower() != "gene":
                    continue
                identifier = _parse_attributes(fields[8]).get("ID", "")
                identifier = re.sub(r"^gene:", "", identifier, flags=re.I)
                match = re.fullmatch(r"(.+?)(\d+)([gG])(\d+)", identifier)
                if not match:
                    return None
                prefix, chrom, separator, number = match.groups()
                formats.add((prefix, separator, len(number), len(chrom)))
                chromosomes.add(chrom); samples += 1
        if not formats or len(formats) != 1:
            return None
        prefix, separator, digits, chrom_width = next(iter(formats))
        with open(str(genome) + ".fai", encoding="utf-8") as handle:
            for index, line in enumerate(handle):
                if index >= 10000:
                    break
                match = re.fullmatch(r"(?:chr|chromosome)?(\d+)(?:LG\d+)?", line.split("\t")[0], re.I)
                if match:
                    chromosomes.add(str(int(match.group(1))).zfill(chrom_width))
                    if len(chromosomes) > 64:
                        return None
        return {"prefix": prefix, "separator": separator, "digits": digits,
                "chromosomes": sorted(chromosomes, key=int), "source": "annotation_sample"}
    except (OSError, ValueError, EOFError):
        return None


def _profile(name, genome, gff3=None, database=None, gene_id_format=None):
    if gene_id_format is not None:
        if (not isinstance(gene_id_format, dict) or
                not isinstance(gene_id_format.get("prefix"), str) or
                not isinstance(gene_id_format.get("separator"), str) or
                type(gene_id_format.get("digits")) is not int or
                not 1 <= gene_id_format["digits"] <= 12 or
                not isinstance(gene_id_format.get("chromosomes"), list) or
                not gene_id_format["chromosomes"] or
                any(not isinstance(chrom, str) or not re.fullmatch(r"\d+", chrom)
                    for chrom in gene_id_format["chromosomes"])):
            raise ValueError("invalid gene_id_format in reference profile")
        gene_id_format = dict(gene_id_format, source="registered")
    genome = str(Path(genome).expanduser())
    gff3 = str(Path(gff3).expanduser()) if gff3 else None
    database = str(Path(database).expanduser()) if database else None
    missing = []
    if not Path(genome).is_file():
        missing.append("FASTA")
    if not Path(genome + ".fai").is_file():
        missing.append("FASTA index (.fai)")
    if gff3 and not Path(gff3).is_file():
        missing.append("GFF3")
    if database and not Path(str(database) + ".nin").is_file():
        missing.append("BLAST database")
    return {"name": str(name), "genome": genome, "gff3": str(gff3) if gff3 else "",
            "database": str(database) if database else "", "available": not missing,
            "missing": missing, "gene_id_format": gene_id_format or
            (infer_gene_id_format(genome, gff3) if not missing else None)}


def reference_catalog(databases, config_path=None):
    path = Path(config_path or os.environ.get("PRIMERBLAST_REFERENCES") or
                Path.home() / ".codex" / "primerblast-oss" / "references.json")
    references, warnings = [], []
    if path.is_file():
        try:
            profiles = json.loads(path.read_text(encoding="utf-8"))
            if not isinstance(profiles, list):
                raise ValueError("expected a list of reference profiles")
            for profile in profiles:
                if not isinstance(profile, dict) or not profile.get("genome"):
                    raise ValueError("each profile requires genome")
                references.append(_profile(profile.get("name") or Path(profile["genome"]).name,
                                           profile["genome"], profile.get("gff3"), profile.get("database"), profile.get("gene_id_format")))
        except (OSError, ValueError, TypeError) as error:
            warnings.append(f"Reference catalog {path}: {error}")
    seen = {profile["genome"] for profile in references}
    for database in databases:
        prefix = Path(database["path"])
        stem = prefix.with_suffix("") if prefix.suffix in (".fa", ".fasta", ".fna") else prefix
        candidates = [prefix] + [Path(str(stem) + ext) for ext in (".fa", ".fasta", ".fna")]
        fasta = {str(candidate) for candidate in candidates
                 if candidate.is_file() and Path(str(candidate) + ".fai").is_file()}
        if len(fasta) != 1:
            continue
        genome = next(iter(fasta))
        if genome in seen:
            continue
        annotations = [Path(str(stem) + ext) for ext in (".gff3", ".gff3.gz", ".gff", ".gff.gz")]
        annotations = [str(candidate) for candidate in annotations if candidate.is_file()]
        references.append(_profile(database["name"], genome,
                                   annotations[0] if len(annotations) == 1 else None, str(prefix)))
        seen.add(genome)
    return {"references": references, "warnings": warnings}
