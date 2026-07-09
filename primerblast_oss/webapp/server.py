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
POST /api/run/<mode>       -> {job_id}          (mode = design|check|tile|assay|markers|makedb)
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
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple
from urllib.parse import unquote, urlparse

from ..design import DesignParams, clean_sequence, read_fasta
from ..specificity import SpecParams, in_silico_pcr
from ..pipeline import run_pipeline
from ..tiling import design_tiling
from ..tools import make_blastdb
from .. import __version__
from .. import report as R
from .. import outputs as OUT

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
    )


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
                name = line[1:].split()[0] if len(line) > 1 else "seq"
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
    """Fast pre-scan: return the chromosome (col 1) of the first GFF3 line
    mentioning the gene, so the full parse can be bounded to one seqid.

    All of a gene's features share a chromosome, so a substring match on any
    of its lines (gene/mRNA/exon/CDS) yields the right seqid."""
    import gzip
    opener = gzip.open if gff3_path.endswith(".gz") else open
    id_tag, name_tag = f"ID={gene}", f"Name={gene}"
    try:
        with opener(gff3_path, "rt", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line and line[0] != "#" and (id_tag in line or name_tag in line):
                    return line.split("\t", 1)[0] or None
    except OSError:
        return None
    return None


def _gene_to_template(p: Dict) -> Tuple[str, str]:
    """Resolve a gene name/ID (+GFF3 +genome) into a (template_id, sequence)."""
    from ..genome import Genome
    from ..regions import resolve_gene, extract_template

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
    return (tmpl.id, tmpl.seq)


def _run_design(p: Dict) -> Dict:
    dp = _design_params(p)
    sp = _spec_params(p)
    dbs = _databases(p)
    size_tol = int(_f(p, "size_tolerance", 10))
    if _f(p, "source", "sequence") == "gene":
        templates = [_gene_to_template(p)]
    else:
        templates = _templates(p)
    results = []
    for tid, seq in templates:
        res = run_pipeline(tid, seq, dbs, design_params=dp, spec_params=sp,
                           size_tolerance=size_tol)
        d = R.to_dict(res)
        d["tsv"] = R.to_tsv(res)
        results.append(d)
    return {"mode": "design", "templates": results}


def _run_check(p: Dict) -> Dict:
    sp = _spec_params(p)
    dbs = _databases(p)
    primers: Dict[str, str] = {}
    if p.get("forward"):
        primers["F"] = str(p["forward"]).upper().strip()
    if p.get("reverse"):
        primers["R"] = str(p["reverse"]).upper().strip()
    for i, spec in enumerate(p.get("primers") or [], 1):
        spec = str(spec).strip()
        if not spec:
            continue
        if "=" in spec:
            name, seq = spec.split("=", 1)
        else:
            name, seq = f"P{i}", spec
        primers[name.strip()] = seq.upper().strip()
    if not primers:
        raise ValueError("Provide at least a forward/reverse primer or a primer list.")
    results = [in_silico_pcr(primers, db, sp=sp) for db in dbs]
    return {"mode": "check", **R.insilico_to_dict(results, primers)}


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
        )
        reg = region or (0, len(clean_sequence(seq)) - 1)
        out.append(R.tiling_to_dict(tiles, tid, reg, dbs))
    return {"mode": "tile", "templates": out}


def _run_assay(p: Dict) -> Dict:
    from ..genome import Genome
    from ..regions import resolve_gene, resolve_interval, resolve_snp
    from ..assay import run_assay
    from ..provenance import make_manifest
    from ..vcf import parse_vcf

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
                              feature=_f(p, "gene_feature", "cds"), flank=0)
    elif _f(p, "interval", None):
        chrom, span = str(p["interval"]).split(":")
        s, e = span.split("-")
        region = resolve_interval(chrom, int(s), int(e),
                                  name=_f(p, "name", p["interval"]))
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
                       spec_params=sp, variants=variants, caps_snp=caps_snp)
    prov = make_manifest({"design": dp.__dict__, "spec": sp.__dict__, "flank": flank},
                         dbs, template_info=result["target"])
    result["provenance"] = prov
    pairs = result.get("pairs", [])
    result["exports"] = {
        "csv": _safe(OUT.pairs_to_csv, pairs),
        "order": _safe(OUT.order_table, pairs),
        "bed": _safe(OUT.products_to_bed, pairs),
        "ascii": "\n\n".join(_safe(OUT.ascii_offtarget_map, pr) or "" for pr in pairs),
    }
    result["mode"] = "assay"
    return result


def _run_markers(p: Dict) -> Dict:
    from ..genome import Genome
    from ..regions import resolve_interval
    from ..assay import design_qtl_markers

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
    "design": _run_design,
    "check": _run_check,
    "tile": _run_tile,
    "assay": _run_assay,
    "markers": _run_markers,
    "makedb": _run_makedb,
}


# --------------------------------------------------------------------------- #
# job manager
# --------------------------------------------------------------------------- #
class JobManager:
    def __init__(self) -> None:
        self._jobs: Dict[str, Dict] = {}
        self._lock = threading.Lock()

    def submit(self, mode: str, params: Dict) -> str:
        handler = HANDLERS.get(mode)
        if handler is None:
            raise ValueError(f"unknown mode: {mode}")
        job_id = uuid.uuid4().hex
        with self._lock:
            self._jobs[job_id] = {"status": "running", "mode": mode}
        t = threading.Thread(target=self._work, args=(job_id, handler, params),
                             daemon=True)
        t.start()
        return job_id

    def _work(self, job_id: str, handler: Callable, params: Dict) -> None:
        try:
            result = handler(params)
            with self._lock:
                self._jobs[job_id].update(status="done", result=result)
        except Exception as exc:  # noqa: BLE001 - surface any engine error to UI
            with self._lock:
                self._jobs[job_id].update(
                    status="error", error=str(exc),
                    trace=traceback.format_exc())

    def get(self, job_id: str) -> Optional[Dict]:
        with self._lock:
            job = self._jobs.get(job_id)
            return dict(job) if job else None


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
