"""Machine-readable local CLI over the same workflows used by the GUI."""
import contextlib
import json
from pathlib import Path
import sys

from .workflows import HANDLERS, execute

SCHEMA_VERSION = "primerblast-agent/1"


def _reject_constant(value):
    raise ValueError("Non-finite JSON number: " + value)

COMMON_INPUT = {
    "db": {"type": "array", "items": {"type": "string"}},
    "template": {"type": "string", "description": "DNA or pasted FASTA"},
    "template_id": {"type": "string"},
    "gene": {"type": "string"}, "gff3": {"type": "string"},
    "genome": {"type": "string", "description": "Local FASTA with an existing .fai"},
    "interval": {"type": "string", "description": "chrom:start-end (1-based inclusive)"},
    "db_genomes": {"type": ["object", "string"], "description": "Explicit DB-to-FASTA associations"},
    "num_threads": {"type": "integer", "minimum": 1},
}
MODE_INPUT = {
    "design": {"product_size": {"type": "string"}, "num_return": {"type": "integer", "minimum": 1}},
    "check": {"forward": {"type": "string"}, "reverse": {"type": "string"}, "primers": {"type": "array", "items": {"type": "string"}}},
    "tile": {"amplicon_min": {"type": "integer"}, "amplicon_max": {"type": "integer"}, "overlap": {"type": "integer"}},
    "sequence": {"source": {"enum": ["sequence", "gene", "interval"]}, "amplicon_size": {"type": "string"}, "m13_tails": {"type": "boolean"}},
    "assay": {"snp": {"type": "string"}, "alt": {"type": "string"}},
    "markers": {"n_markers": {"type": "integer"}, "spacing": {"type": "integer"}},
    "makedb": {"infile": {"type": "string"}, "out_db": {"type": "string"}, "title": {"type": "string"}},
}


def schema():
    return {
        "schema_version": SCHEMA_VERSION,
        "operations": [{"name": name, "writes_database": name == "makedb",
                        "input_schema": {"type": "object", "properties": {
                            **(COMMON_INPUT if name != "makedb" else {}), **MODE_INPUT[name]},
                            "additionalProperties": True},
                        "input": "GUI parameters; primary fields described above; full options in native CLI help"}
                       for name in HANDLERS],
        "request_schema": {"type": "object", "required": ["operation", "params"],
                           "additionalProperties": False, "properties": {
            "operation": {"enum": list(HANDLERS)}, "params": {"type": "object"},
            "request_id": {"type": "string"}}},
        "execution": "local", "starts_gui": False, "uploads_data": False,
    }


def command(args):
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    if args.action == "schema":
        print(json.dumps(schema(), ensure_ascii=False))
        return 0
    request = {}
    try:
        raw = sys.stdin.read() if args.input == "-" else Path(args.input).read_text(encoding="utf-8-sig")
        request = json.loads(raw.lstrip("\ufeff"), parse_constant=_reject_constant)
        if not isinstance(request, dict):
            raise ValueError("Request must be a JSON object")
        if set(request) - {"operation", "params", "request_id"}:
            raise ValueError("Unknown request fields")
        if not isinstance(request.get("operation"), str) or not isinstance(request.get("params"), dict):
            raise ValueError("operation and params are required")
        if "request_id" in request and not isinstance(request["request_id"], str):
            raise ValueError("request_id must be a string")
        # Native adapters may print diagnostics; stdout stays exactly one JSON value.
        with contextlib.redirect_stdout(sys.stderr):
            result = execute(request["operation"], request["params"], allow_db_write=args.allow_db_write)
        output = {"schema_version": SCHEMA_VERSION, "ok": True,
                  "operation": request["operation"], "result": result}
        if "request_id" in request:
            output["request_id"] = request["request_id"]
        try:
            encoded = json.dumps(output, ensure_ascii=False, allow_nan=False)
        except (TypeError, ValueError) as exc:
            raise RuntimeError("Operation returned a result that cannot be encoded as JSON") from exc
        code = 0
    except Exception as exc:
        output = {"schema_version": SCHEMA_VERSION, "ok": False,
                  "error": {"code": "invalid_request" if isinstance(exc, (ValueError, KeyError, json.JSONDecodeError)) else "execution_failed",
                            "message": str(exc)}}
        code = 2 if output["error"]["code"] == "invalid_request" else 1
        if isinstance(request, dict) and isinstance(request.get("request_id"), str):
            output["request_id"] = request["request_id"]
        encoded = json.dumps(output, ensure_ascii=False, allow_nan=False)
    print(encoded)
    return code
