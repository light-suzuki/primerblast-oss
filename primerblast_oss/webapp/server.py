"""Stdlib-only HTTP backend for the primerblast-oss web GUI.

No third-party dependencies: the whole server is built on ``http.server`` so it
inherits the package's "no runtime deps" contract. Long-running BLAST jobs run
in background threads; the browser submits a job and polls for the result.

Endpoints
---------
GET  /                     -> static UI (index.html)
GET  /<asset>              -> static asset (app.js, style.css, i18n.js, ...)
GET  /api/health           -> tool availability + versions
GET  /api/databases        -> discovered BLAST nucleotide databases
GET  /api/references       -> local reference paths and gene-ID format metadata
POST /api/run/<mode>       -> {job_id}          (mode = design|check|tile|sequence|assay|markers|makedb)
GET  /api/job/<job_id>     -> {status, result?, error?}
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import threading
import traceback
import uuid
import copy
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Callable, Dict, List, Optional
from urllib.parse import unquote, urlparse

from .. import __version__
from .. import progress

STATIC_DIR = Path(__file__).resolve().parent / "static"

# Directories scanned for pre-built BLAST databases. The first that exists wins
# for discovery but all are scanned; users can also type an absolute path.
DEFAULT_DB_DIRS = [
    Path.home() / ".codex" / "blast_databases",
    Path.home() / "blast_databases",
    Path.home() / "primerblast-oss" / "databases",
]

# --------------------------------------------------------------------------- #
# parameter builders (mirror the CLI defaults)
# --------------------------------------------------------------------------- #
# Compatibility aliases for earlier integrations; execution lives in workflows.
from ..workflows import (
    HANDLERS, execute, _associated_genomes, _find_gene_seqid, _run_sequence,
    _run_design, _run_check, _run_tile, _run_assay, _run_markers, _run_makedb,
)


# --------------------------------------------------------------------------- #
# job manager
# --------------------------------------------------------------------------- #
class JobManager:
    def __init__(self) -> None:
        self._jobs: Dict[str, Dict] = {}
        self._lock = threading.Lock()

    def submit(self, mode: str, params: Dict) -> str:
        if mode not in HANDLERS:
            raise ValueError(f"unknown mode: {mode}")
        if not isinstance(params, dict):
            raise ValueError("params must be a JSON object")
        handler = lambda payload: execute(mode, payload, allow_db_write=mode == "makedb")
        job_id = uuid.uuid4().hex
        with self._lock:
            self._jobs[job_id] = {"status": "running", "mode": mode,
                                  "started_at": time.time(), "partial_revision": 0}
        t = threading.Thread(target=self._work, args=(job_id, handler, params),
                             daemon=True)
        t.start()
        return job_id

    def _work(self, job_id: str, handler: Callable, params: Dict) -> None:
        def publish(event):
            with self._lock:
                job = self._jobs[job_id]
                if 'partial_result' in event:
                    job['partial_result'] = copy.deepcopy(event['partial_result'])
                    job['partial_revision'] += 1
                else:
                    job.update(event)
                    job['_progress_at'] = time.monotonic()
        try:
            with progress.observe(publish) as recorder:
                progress.report('prepare')
                result = handler(params)
                timing = recorder.finish()
            with self._lock:
                job = self._jobs[job_id]
                job.update(status="done", result=result, timing=timing,
                           finished_at=time.time())
                job.pop('partial_result', None)
        except Exception as exc:  # noqa: BLE001 - surface any engine error to UI
            timing = recorder.finish('error')
            with self._lock:
                self._jobs[job_id].update(
                    status="error", error=str(exc),
                    timing=timing, finished_at=time.time(),
                    trace=traceback.format_exc())

    def get(self, job_id: str) -> Optional[Dict]:
        with self._lock:
            job = self._jobs.get(job_id)
            if not job:
                return None
            snapshot = dict(job)
            snapshot.pop('_progress_at', None)
            if job['status'] == 'running' and 'progress' in job:
                observation = dict(job['progress'])
                delta = max(0, time.monotonic() - job['_progress_at'])
                observation['elapsed_seconds'] += delta
                observation['stage_elapsed_seconds'] += delta
                snapshot['progress'] = observation
            return snapshot


JOBS = JobManager()


# --------------------------------------------------------------------------- #
# environment probes
# --------------------------------------------------------------------------- #
def _tool_version(binary: str, args: List[str]) -> Optional[str]:
    path = shutil.which(binary)
    if not path:
        return None
    try:
        out = subprocess.run([path] + args, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, timeout=8)
        first = out.stdout.decode(errors="ignore").strip().splitlines()
        return first[0] if first else path
    except Exception:
        return path


def health() -> Dict:
    return {
        "app_version": __version__,
        "tools": {
            "primer3_core": _tool_version("primer3_core", ["-about"]) or None,
            "blastn": _tool_version("blastn", ["-version"]) or None,
            "makeblastdb": _tool_version("makeblastdb", ["-version"]) or None,
        },
        "ok": bool(shutil.which("primer3_core") and shutil.which("blastn")),
    }


def discover_databases() -> List[Dict]:
    seen = set()
    found: List[Dict] = []
    for d in DEFAULT_DB_DIRS:
        if not d.is_dir():
            continue
        for nin in sorted(d.glob("*.nin")):
            prefix = str(nin.with_suffix(""))
            if prefix in seen:
                continue
            seen.add(prefix)
            found.append({"name": nin.stem, "path": prefix, "dir": str(d)})
    return found


# --------------------------------------------------------------------------- #
# HTTP handler
# --------------------------------------------------------------------------- #
_CONTENT_TYPES = {
    ".html": "text/html; charset=utf-8",
    ".js": "application/javascript; charset=utf-8",
    ".css": "text/css; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".svg": "image/svg+xml",
    ".ico": "image/x-icon",
}


class Handler(BaseHTTPRequestHandler):
    server_version = "primerblast-oss-webapp"

    # -- helpers -----------------------------------------------------------
    def _send_json(self, obj, code: int = 200) -> None:
        body = json.dumps(obj, default=str).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _send_static(self, rel: str) -> None:
        if rel in ("", "/"):
            rel = "index.html"
        rel = rel.lstrip("/")
        target = (STATIC_DIR / rel).resolve()
        if not str(target).startswith(str(STATIC_DIR.resolve())) or not target.is_file():
            self.send_error(404, "Not found")
            return
        body = target.read_bytes()
        self.send_response(200)
        self.send_header("Content-Type",
                         _CONTENT_TYPES.get(target.suffix, "application/octet-stream"))
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _read_body(self) -> Dict:
        length = int(self.headers.get("Content-Length", 0) or 0)
        if not length:
            return {}
        raw = self.rfile.read(length)
        try:
            return json.loads(raw.decode("utf-8"))
        except Exception:
            return {}

    # -- routes ------------------------------------------------------------
    def do_GET(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        try:
            if path == "/api/health":
                self._send_json(health())
            elif path == "/api/databases":
                self._send_json({"databases": discover_databases()})
            elif path == "/api/references":
                from .references import reference_catalog
                self._send_json(reference_catalog(discover_databases()))
            elif path == "/api/enzymes":
                from ..restriction_catalog import catalog
                self._send_json(catalog())
            elif path.startswith("/api/job/"):
                job_id = unquote(path[len("/api/job/"):])
                job = JOBS.get(job_id)
                if job is None:
                    self._send_json({"error": "unknown job"}, 404)
                else:
                    self._send_json(job)
            else:
                self._send_static(path)
        except Exception as exc:  # noqa: BLE001
            self._send_json({"error": str(exc)}, 500)

    def do_POST(self) -> None:  # noqa: N802
        path = urlparse(self.path).path
        try:
            if path.startswith("/api/run/"):
                mode = unquote(path[len("/api/run/"):])
                body = self._read_body()
                params = body.get("params", body)
                job_id = JOBS.submit(mode, params)
                self._send_json({"job_id": job_id})
            else:
                self._send_json({"error": "not found"}, 404)
        except Exception as exc:  # noqa: BLE001
            self._send_json({"error": str(exc)}, 400)

    def log_message(self, fmt, *args):  # quiet by default
        if os.environ.get("PRIMERBLAST_WEB_VERBOSE"):
            super().log_message(fmt, *args)


def serve(host: str = "127.0.0.1", port: int = 8799) -> ThreadingHTTPServer:
    httpd = ThreadingHTTPServer((host, port), Handler)
    return httpd
