"""Reproducibility manifests for tools, scientific models and reference inputs.

The default fingerprint policy is deliberately lightweight so an ordinary assay
does not reread multi-gigabase plant genomes just to build its report. Callers
can request strong_hashes=True to stream the complete referenced files through
SHA-256 for archival/publication-grade provenance.
"""
from __future__ import annotations

import hashlib
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from . import __version__


PROVENANCE_SCHEMA_VERSION = "2.0"
MODEL_VERSIONS = {
    "specificity": "full-primer-realignment-v1",
    "gel": "agarose-log-gap-v1",
    "risk": "heuristic-score-v1",
    "restriction": "rebase404-cleavage-v1",
}


def model_versions() -> Dict[str, str]:
    """Return version identifiers for scientific decision models."""
    return dict(MODEL_VERSIONS)


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


def git_revision() -> Optional[str]:
    """Return the source revision when discoverable without requiring Git."""
    for key in ("PRIMERBLAST_GIT_REVISION", "GITHUB_SHA"):
        value = os.environ.get(key)
        if value:
            return value.strip() or None
    git = shutil.which("git")
    if not git:
        return None
    try:
        root = Path(__file__).resolve().parents[1]
        process = subprocess.run(
            [git, "rev-parse", "HEAD"], cwd=str(root),
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, timeout=3)
        if process.returncode == 0:
            value = process.stdout.decode(errors="ignore").strip()
            return value or None
    except (OSError, subprocess.SubprocessError):
        pass
    return None


def tool_versions() -> Dict[str, str]:
    return {
        "primerblast_oss": __version__,
        "python": sys.version.split()[0],
        "primer3_core": _run_version(["primer3_core", "--version"]),
        "blastn": _run_version(["blastn", "-version"]),
        "makeblastdb": _run_version(["makeblastdb", "-version"]),
    }


def sha1_file(path: str, limit_bytes: int = 4_000_000) -> Optional[str]:
    """Compatibility quick fingerprint over only the leading file bytes."""
    file_path = Path(path)
    if not file_path.exists():
        return None
    digest = hashlib.sha1()
    with file_path.open("rb") as handle:
        digest.update(handle.read(limit_bytes))
    return digest.hexdigest()


def sha256_file(path: str, chunk_bytes: int = 4 * 1024 * 1024) -> Optional[str]:
    """Stream an entire file through SHA-256 without materializing it in RAM."""
    file_path = Path(path)
    if not file_path.is_file():
        return None
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        while True:
            block = handle.read(chunk_bytes)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _timestamp(stat) -> str:
    return datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()


def file_fingerprint(path: Optional[str], *, strong: bool = False) -> Dict:
    """Fingerprint an arbitrary reference file."""
    if not path:
        return {
            "path": None, "bytes": None, "mtime": None,
            "sha1_prefix": None, "sha256": None,
        }
    file_path = Path(path)
    if not file_path.is_file():
        return {
            "path": path, "bytes": None, "mtime": None,
            "sha1_prefix": None, "sha256": None,
        }
    stat = file_path.stat()
    return {
        "path": path,
        "bytes": str(stat.st_size),
        "mtime": _timestamp(stat),
        "sha1_prefix": sha1_file(path),
        "sha256": sha256_file(path) if strong else None,
    }


def fasta_fingerprint(path: Optional[str], *, strong: bool = False) -> Dict:
    """Compatibility FASTA-shaped wrapper around file_fingerprint."""
    result = file_fingerprint(path, strong=strong)
    return {
        "fasta": result["path"],
        "bytes": result["bytes"],
        "mtime": result["mtime"],
        "sha1_prefix": result["sha1_prefix"],
        "sha256": result["sha256"],
    }


def _blast_index_paths(db_path: str) -> List[Path]:
    """Return BLAST nucleotide index files for a database prefix."""
    prefix = Path(db_path)
    parent = prefix.parent if str(prefix.parent) else Path(".")
    stem = prefix.name
    matches = []
    for candidate in parent.glob(stem + ".*"):
        suffix = candidate.name[len(stem):]
        if candidate.is_file() and re.match(
                r"^(?:\.\d+)?\.n[a-z0-9]+$", suffix):
            matches.append(candidate)
    return sorted(matches, key=lambda item: item.name)


def db_fingerprint(db_path: str, *, strong: bool = False) -> Dict:
    """Fingerprint a BLAST database while retaining legacy summary fields."""
    index_files = _blast_index_paths(db_path)
    info: Dict = {
        "db": db_path,
        "index_bytes": None,
        "mtime": None,
        "files": [],
    }
    summary_path = None
    for extension in (".nin", ".nsq", ".ndb"):
        candidate = Path(db_path + extension)
        if candidate.is_file():
            summary_path = candidate
            break
    if summary_path is None and index_files:
        summary_path = index_files[0]
    if summary_path is not None:
        stat = summary_path.stat()
        info["index_bytes"] = str(stat.st_size)
        info["mtime"] = _timestamp(stat)
    for path in index_files:
        stat = path.stat()
        info["files"].append({
            "name": path.name,
            "bytes": str(stat.st_size),
            "mtime": _timestamp(stat),
            "sha256": sha256_file(str(path)) if strong else None,
        })
    return info


def make_manifest(params: Dict, databases: Sequence[str],
                  template_info: Optional[Dict] = None,
                  now: Optional[datetime] = None,
                  thermo_genomes: Optional[Mapping[str, Optional[str]]] = None,
                  reference_files: Optional[Mapping[str, Optional[str]]] = None,
                  strong_hashes: bool = False) -> Dict:
    """Build one versioned provenance manifest.

    Strong hashing is opt-in so normal interactive runs do not gain a full
    multi-gigabyte read of every reference file.
    """
    timestamp = (now or datetime.now(timezone.utc)).isoformat()
    fasta_mapping = thermo_genomes
    if fasta_mapping is None:
        candidate = params.get("thermo_genomes") if isinstance(params, dict) else None
        fasta_mapping = candidate if isinstance(candidate, dict) else {}
    references = {
        name: file_fingerprint(path, strong=strong_hashes)
        for name, path in (reference_files or {}).items()
    }
    return {
        "provenance_schema_version": PROVENANCE_SCHEMA_VERSION,
        "generated": timestamp,
        "source_revision": git_revision(),
        "model_versions": model_versions(),
        "fingerprint_policy": {
            "mode": "full_sha256" if strong_hashes else "lightweight",
            "full_sha256": bool(strong_hashes),
            "lightweight_sha1_prefix_bytes": 4_000_000,
        },
        "tool_versions": tool_versions(),
        "databases": [
            db_fingerprint(database, strong=strong_hashes)
            for database in databases
        ],
        "thermo_genomes": {
            database: fasta_fingerprint(
                (fasta_mapping or {}).get(database), strong=strong_hashes)
            for database in databases
        },
        "reference_files": references,
        "parameters": params,
        "template": template_info or {},
        "host": (
            os.uname().nodename if hasattr(os, "uname")
            else os.environ.get("COMPUTERNAME", "")
        ),
    }
