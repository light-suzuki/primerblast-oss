"""Reproducibility manifest for tools, parameters, BLAST DBs and FASTA files."""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from . import __version__


def _run_version(command: List[str]) -> str:
    executable = shutil.which(command[0])
    if not executable:
        return "not found"
    try:
        process = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            timeout=15)
        return (
            process.stdout.decode(errors="ignore").strip().splitlines()[0]
            if process.stdout else ""
        )
    except Exception as error:  # noqa: BLE001
        return "error: %s" % error


def tool_versions() -> Dict[str, str]:
    return {
        "primerblast_oss": __version__,
        "python": sys.version.split()[0],
        "primer3_core": _run_version(["primer3_core", "--version"]),
        "blastn": _run_version(["blastn", "-version"]),
        "makeblastdb": _run_version(["makeblastdb", "-version"]),
    }


def db_fingerprint(db_path: str) -> Dict[str, Optional[str]]:
    info: Dict[str, Optional[str]] = {
        "db": db_path, "index_bytes": None, "mtime": None}
    for extension in (".nin", ".nsq", ".ndb"):
        path = Path(db_path + extension)
        if path.exists():
            stat = path.stat()
            info["index_bytes"] = str(stat.st_size)
            info["mtime"] = datetime.fromtimestamp(
                stat.st_mtime, timezone.utc).isoformat()
            break
    return info


def sha1_file(path: str, limit_bytes: int = 4_000_000) -> Optional[str]:
    file_path = Path(path)
    if not file_path.exists():
        return None
    digest = hashlib.sha1()
    with file_path.open("rb") as handle:
        digest.update(handle.read(limit_bytes))
    return digest.hexdigest()


def fasta_fingerprint(path: Optional[str]) -> Dict[str, Optional[str]]:
    if not path:
        return {"fasta": None, "bytes": None, "mtime": None, "sha1_prefix": None}
    file_path = Path(path)
    if not file_path.exists():
        return {"fasta": path, "bytes": None, "mtime": None, "sha1_prefix": None}
    stat = file_path.stat()
    return {
        "fasta": path,
        "bytes": str(stat.st_size),
        "mtime": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        "sha1_prefix": sha1_file(path),
    }


def make_manifest(params: Dict, databases: Sequence[str],
                  template_info: Optional[Dict] = None,
                  now: Optional[datetime] = None,
                  thermo_genomes: Optional[Mapping[str, Optional[str]]] = None) -> Dict:
    timestamp = (now or datetime.now(timezone.utc)).isoformat()
    fasta_mapping = thermo_genomes
    if fasta_mapping is None:
        candidate = params.get("thermo_genomes") if isinstance(params, dict) else None
        fasta_mapping = candidate if isinstance(candidate, dict) else {}
    return {
        "generated": timestamp,
        "tool_versions": tool_versions(),
        "databases": [db_fingerprint(database) for database in databases],
        "thermo_genomes": {
            database: fasta_fingerprint((fasta_mapping or {}).get(database))
            for database in databases
        },
        "parameters": params,
        "template": template_info or {},
        "host": (
            os.uname().nodename if hasattr(os, "uname")
            else os.environ.get("COMPUTERNAME", "")
        ),
    }
