"""Published reports must not contain personal filesystem or host metadata."""
import json
import re
from pathlib import Path


def test_published_benchmark_reports_have_no_personal_paths():
    root = Path(__file__).resolve().parents[1] / "benchmarks"
    pattern = re.compile(r"/home/[^/\s]+/|[A-Za-z]:[\\/]+Users[\\/]+[^\\/\s]+[\\/]", re.I)
    for path in root.rglob("*"):
        if path.suffix in (".json", ".html", ".csv", ".bed", ".md"):
            assert not pattern.search(path.read_text(encoding="utf-8")), path.name


def test_published_host_metadata_is_redacted():
    root = Path(__file__).resolve().parents[1] / "benchmarks"
    def check(value):
        if isinstance(value, dict):
            if "host" in value:
                assert value["host"] in ("local-workstation", "redacted")
            for item in value.values():
                check(item)
        elif isinstance(value, list):
            for item in value:
                check(item)
    for path in root.rglob("*.json"):
        check(json.loads(path.read_text(encoding="utf-8")))
