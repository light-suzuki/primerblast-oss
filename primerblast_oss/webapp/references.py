"""Small reference catalog: explicit profiles or unambiguous sibling files."""
import json
import os
from pathlib import Path


def _profile(name, genome, gff3=None, database=None):
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
            "missing": missing}


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
                                           profile["genome"], profile.get("gff3"), profile.get("database")))
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
