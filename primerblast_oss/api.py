"""Stable public library API for downstream applications (e.g. SnapyGene).

This module is the versioned contract described in ``docs/PUBLIC_API.md``.
Functions here return JSON-safe plain dicts (no dataclass instances) and
every result carries an ``api_version`` key so persisted output can be
validated. The lower-level functions SnapyGene's adapter already calls
(``design.DesignParams``, ``pipeline.run_pipeline``,
``specificity.spec_params_for_profile`` / ``pair_specificity`` /
``in_silico_pcr``) are unchanged and remain importable from their original
modules.

Coordinate and strand conventions (documented in ``docs/PUBLIC_API.md``):
template ``left_start``/``right_start`` are 0-based template positions, all
genomic product coordinates are 1-based inclusive, and ``strand`` is
``+``/``-`` (extension direction is derived per product).
"""
from __future__ import annotations

import shutil
import subprocess
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence

from . import __version__
from .design import PrimerPair
from .errors import (
    CancelledError,
    SearchIncompleteError,
)
from .pipeline import run_pipeline
from .specificity import (
    SEARCH_COMPLETE,
    PrimingSite,
    SpecParams,
    Amplicon,
    in_silico_pcr,
    pair_specificity,
    spec_params_for_profile,
)
from .tools import make_blastdb

API_VERSION = "1.0"
"""Schema version of the JSON-safe result models in this module.

Persisted results carry this value in their ``api_version`` key; bump it
when a result field is added, removed, or redefined in a way that is not a
strict superset.
"""


def json_safe(value: Any) -> Any:
    """Deeply convert dataclass instances to plain dicts for JSON output.

    Nested containers are converted recursively; everything else (str, int,
    float, bool, None) is passed through unchanged. The result of this
    function can be passed to ``json.dumps`` directly.
    """
    if isinstance(value, (PrimingSite, Amplicon, PrimerPair)):
        return json_safe(value.__dict__)
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def capabilities() -> Dict[str, Any]:
    """Report what this installed package can do.

    ``multiplex`` and ``thermodynamic_filtering`` depend on the optional
    ``primer3-py`` package, so they are reported dynamically.
    """
    from . import thermo as thermo_module

    thermo = bool(thermo_module.available())
    return {
        "api_version": API_VERSION,
        "package_version": __version__,
        "capabilities": {
            "design": True,
            "pair_specificity": True,
            "pool_in_silico_pcr": True,
            "tiling": True,
            "multiplex": thermo,
            "qpcr_probe_design": False,
            "blast_database_creation": True,
            "thermodynamic_filtering": thermo,
        },
    }


def discover_tools() -> Dict[str, Any]:
    """Locate the external executables and report actionable setup messages.

    Returns ``{"tools": {name: {"path", "version", "available"}},
    "missing": [...], "messages": {name: hint}, "complete": bool}``.
    """
    tools: Dict[str, Dict[str, Any]] = {}
    for name, version_args in (
        ("primer3_core", ["--version"]),
        ("blastn", ["-version"]),
        ("makeblastdb", ["-version"]),
    ):
        path = shutil.which(name)
        version: Optional[str] = None
        if path:
            try:
                process = subprocess.run(
                    [path] + version_args,
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=15,
                )
                first_line = (
                    process.stdout.decode(errors="ignore").strip().splitlines()
                    if process.stdout else []
                )
                version = first_line[0] if first_line else ""
            except (OSError, subprocess.SubprocessError):
                version = None
        tools[name] = {"path": path, "version": version, "available": bool(path)}
    missing = [name for name, info in tools.items() if not info["available"]]
    messages = {
        "primer3_core": (
            "Install the 'primer3' package (e.g. `apt install primer3`) or "
            "pass primer3_bin to the API call."),
        "blastn": (
            "Install BLAST+ (e.g. `apt install ncbi-blast+`) or pass "
            "blastn_bin to the API call."),
        "makeblastdb": (
            "Install BLAST+ (e.g. `apt install ncbi-blast+`) or pass "
            "makeblastdb_bin to the API call."),
    }
    return {
        "tools": tools,
        "missing": missing,
        "messages": {name: messages[name] for name in missing},
        "complete": not missing,
    }


def _check_cancelled(cancel_check: Optional[Callable[[], bool]]) -> None:
    if cancel_check is not None and cancel_check():
        raise CancelledError("operation cancelled by caller")


def design_and_screen(
    template_id: str,
    sequence: str,
    databases: Sequence[str],
    design_params=None,
    spec_params: Optional[SpecParams] = None,
    primer3_bin: Optional[str] = None,
    blastn_bin: Optional[str] = None,
    size_tolerance: int = 10,
    genome=None,
    genomes_by_db: Optional[Mapping[str, object]] = None,
    thermo_params=None,
    thermo_gate: bool = True,
    dimer_params=None,
    strict_search: bool = False,
    progress_callback: Optional[Callable[[str, float], None]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """Design primers and screen them against one or more databases.

    Returns a JSON-safe dict with ``api_version``, template metadata, and a
    ``pairs`` list (each pair carries its specificity summary per database).
    With ``strict_search`` True, a search whose evidence is not ``complete``
    raises :class:`~primerblast_oss.errors.SearchIncompleteError` instead of
    returning an indeterminate result.

    ``progress_callback(stage, fraction)`` is called between pairs;
    ``cancel_check()`` returning True aborts with
    :class:`~primerblast_oss.errors.CancelledError`. Cancellation is checked
    between pairs, so it is coarse-grained.
    """
    _check_cancelled(cancel_check)
    result = run_pipeline(
        template_id,
        sequence,
        databases,
        design_params=design_params,
        spec_params=spec_params,
        primer3_bin=primer3_bin,
        blastn_bin=blastn_bin,
        size_tolerance=size_tolerance,
        genome=genome,
        genomes_by_db=genomes_by_db,
        thermo_params=thermo_params,
        thermo_gate=thermo_gate,
        dimer_params=dimer_params,
        progress_callback=progress_callback,
        cancel_check=cancel_check,
    )
    if strict_search:
        incomplete = [
            pair.specificity.get("search_completeness")
            for pair in result.pairs
            if pair.specificity.get("search_completeness") != SEARCH_COMPLETE
        ]
        if incomplete:
            raise SearchIncompleteError(
                "search evidence incomplete for %d pair(s): %s"
                % (len(incomplete), ", ".join(sorted(set(map(str, incomplete))))))
    return {
        "api_version": API_VERSION,
        "template_id": result.template_id,
        "template_len": result.template_len,
        "primer3_explain": result.primer3_explain,
        "databases": result.databases,
        "params": json_safe(result.params),
        "pairs": [json_safe(pair.__dict__) for pair in result.pairs],
    }


def pair_specificity_result(
    forward: str,
    reverse: str,
    db: str,
    designed_size: Optional[int] = None,
    sp: Optional[SpecParams] = None,
    blastn_bin: Optional[str] = None,
    size_tolerance: int = 10,
    genome=None,
    thermo_params=None,
    thermo_gate: bool = True,
    allowed_primer_mismatches: Optional[Mapping[str, int]] = None,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """JSON-safe, versioned wrapper around ``specificity.pair_specificity``."""
    _check_cancelled(cancel_check)
    result = pair_specificity(
        forward,
        reverse,
        db,
        designed_size=designed_size,
        sp=sp,
        blastn_bin=blastn_bin,
        size_tolerance=size_tolerance,
        genome=genome,
        thermo_params=thermo_params,
        thermo_gate=thermo_gate,
        allowed_primer_mismatches=allowed_primer_mismatches,
        cancel_check=cancel_check,
    )
    return {"api_version": API_VERSION, **json_safe(result)}


def pool_in_silico_pcr(
    primers: Dict[str, str],
    db: str,
    sp: Optional[SpecParams] = None,
    blastn_bin: Optional[str] = None,
    genome=None,
    thermo_params=None,
    thermo_gate: bool = True,
    cancel_check: Optional[Callable[[], bool]] = None,
) -> Dict[str, Any]:
    """JSON-safe, versioned wrapper around ``specificity.in_silico_pcr``."""
    _check_cancelled(cancel_check)
    result = in_silico_pcr(
        primers,
        db,
        sp=sp,
        blastn_bin=blastn_bin,
        genome=genome,
        thermo_params=thermo_params,
        thermo_gate=thermo_gate,
        cancel_check=cancel_check,
    )
    return {"api_version": API_VERSION, **json_safe(result)}


def create_database(
    fasta: str,
    out: Optional[str] = None,
    title: Optional[str] = None,
    parse_seqids: bool = True,
    makeblastdb_bin: Optional[str] = None,
) -> Dict[str, Any]:
    """Build a nucleotide BLAST database; returns a JSON-safe result dict."""
    prefix = make_blastdb(
        fasta,
        out=out,
        title=title,
        parse_seqids=parse_seqids,
        makeblastdb_bin=makeblastdb_bin,
    )
    return {
        "api_version": API_VERSION,
        "fasta": fasta,
        "db": prefix,
        "parse_seqids": parse_seqids,
    }


__all__ = [
    "API_VERSION",
    "capabilities",
    "discover_tools",
    "json_safe",
    "design_and_screen",
    "pair_specificity_result",
    "pool_in_silico_pcr",
    "create_database",
    "run_pipeline",
    "pair_specificity",
    "in_silico_pcr",
    "spec_params_for_profile",
    "SpecParams",
    "PrimingSite",
    "Amplicon",
    "PrimerPair",
]
