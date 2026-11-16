"""Guard against truncated/corrupted Python sources (e.g. a lossy transfer
that replaces content with U+2026), which previously broke CI with
``SyntaxError: invalid character '\u2026'``."""
import ast
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ELLIPSIS = "\u2026"


def _python_files():
    for sub in ("primerblast_oss", "tests", "benchmarks", "scripts", "tools"):
        base = ROOT / sub
        if base.is_dir():
            yield from base.rglob("*.py")


def test_python_sources_parse_and_have_no_ellipsis_char():
    problems = []
    for path in _python_files():
        text = path.read_text(encoding="utf-8")
        if ELLIPSIS in text:
            problems.append(f"{path.relative_to(ROOT)}: contains U+2026")
        try:
            ast.parse(text, filename=str(path))
        except SyntaxError as exc:
            problems.append(f"{path.relative_to(ROOT)}: {exc}")
    assert not problems, "\n".join(problems)
