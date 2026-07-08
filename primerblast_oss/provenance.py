"""Reproducibility manifest: pin tool versions, parameters and databases.

Local analysis can be made fully reproducible (unlike NCBI, whose databases
change under you). This captures everything needed to reproduce a run for a
paper / thesis / lab notebook.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from . import __version__


def _run_version(cmd: List[str]) -> str:
    exe = shutil.which(cmd[0])
    if not exe:
        return "not found"
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=15)
        return p.stdout.decode(errors="ignore").strip().splitlines()[0] if p.stdout else ""
    except Exception as e:  # noqa: BLE001
        return f"error: {e}"


def tool_versions() -> Dict[str, str]:
    return {
        "primerblast_oss": __version__,
        "python": sys.version.split()[0],
        "primer3_core": _run_version(["primer3_core", "--version"]),
        "blastn": _run_version(["blastn", "-version"]),
        "makeblastdb": _run_version(["makeblastdb", "-version"]),
    }


def db_fingerprint(db_path: str) -> Dict[str, Optional[str]]:
    """Identify a BLAST db by its index files' size/mtime (cheap, stable)."""
    info: Dict[str, Optional[str]] = {"db": db_path, "nin_bytes": None, "mtime": None}
    for ext in (".nin", ".nsq", ".ndb"):
        p = Path(db_path + ext)
        if p.exists():
            st = p.stat()
            info["nin_bytes"] = str(st.st_size)
            info["mtime"] = datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat()
            break
    return info


def sha1_file(path: str, limit_bytes: int = 4_000_000) -> Optional[str]:
    p = Path(path)
    if not p.exists():
        return None
    h = hashlib.sha1()
    with p.open("rb") as fh:
        h.update(fh.read(limit_bytes))
    return h.hexdigest()


def make_manifest(params: Dict, databases: Sequence[str],
                  template_info: Optional[Dict] = None,
                  now: Optional[datetime] = None) -> Dict:
    ts = (now or datetime.now(timezone.utc)).isoformat()
    return {
        "generated": ts,
        "tool_versions": tool_versions(),
        "databases": [db_fingerprint(db) for db in databases],
        "parameters": params,
        "template": template_info or {},
        "host": os.uname().nodename if hasattr(os, "uname") else os.environ.get("COMPUTERNAME", ""),
    }
