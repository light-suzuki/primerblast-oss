# Public library API

This document defines the **stable library contract** for applications that
embed `primerblast-oss` (notably SnapyGene). Everything described here is
versioned; everything else in the package is internal and may change without
notice.

## Versioning

- Package version: `primerblast_oss.__version__` (SemVer).
- API schema version: `primerblast_oss.API_VERSION` (currently `"1.0"`).
  **Every result dict from the `api` module carries an `api_version` key** so
  persisted output can be validated against the schema that produced it.
- Breaking changes to the public API require a **major** version bump.
- When a public result field is renamed, the old name is kept as a
  compatibility alias for at least **one minor release** before removal.

## What SnapyGene already calls (kept stable)

These symbols are unchanged and remain importable from their original
modules:

| Symbol | Module |
|---|---|
| `DesignParams` | `primerblast_oss.design` |
| `run_pipeline` | `primerblast_oss.pipeline` |
| `spec_params_for_profile` | `primerblast_oss.specificity` |
| `pair_specificity` | `primerblast_oss.specificity` |
| `in_silico_pcr` | `primerblast_oss.specificity` |

`pair_specificity` / `in_silico_pcr` gained an optional `cancel_check`
keyword (additive; existing callers are unaffected).

## Recommended entry point: `primerblast_oss.api`

New integrations should use the `api` module, which returns **JSON-safe plain
dicts** (no dataclass instances; pass them straight to `json.dumps`):

| Function | Returns |
|---|---|
| `capabilities()` | what this install can do (design, pair specificity, pool in-silico PCR, tiling, multiplex, qPCR/probe design, BLAST DB creation, thermodynamic filtering) |
| `discover_tools()` | `primer3_core` / `blastn` / `makeblastdb` path + version, missing-component messages, `complete` flag |
| `design_and_screen(template_id, sequence, databases, ...)` | versioned design-and-screen result |
| `pair_specificity_result(forward, reverse, db, ...)` | versioned pair specificity result |
| `pool_in_silico_pcr(primers, db, ...)` | versioned pool in-silico PCR result |
| `create_database(fasta, out=None, ...)` | versioned database creation result |
| `json_safe(value)` | deep-convert dataclass results to plain dicts |

### Result fields

`design_and_screen` returns `api_version`, `template_id`, `template_len`,
`primer3_explain`, `databases`, `params`, and `pairs` (each pair carries
`forward`/`reverse`, `left_start`/`left_len`/`right_start`/`right_len`
(0-based template positions), `product_size`, `tm_f`/`tm_r`, `gc_f`/`gc_r`,
`penalty`, and `specificity`).

The `specificity` summary per pair includes, per screened database:
template coordinates, primer sequences and coordinates, product
coordinates/sizes, total and 3′-terminal mismatch counts (`fwd_mismatch`,
`rev_mismatch`, `fwd_tp5`, `rev_tp5`), on-target/off-target classification,
search-completeness state, high-copy and target-limit warnings, nearest
off-target size gap, gel distinguishability, score/rank, and database
identity (`db`) plus the parameters used (`blast_limits`).

## Coordinate and strand conventions

| Quantity | Convention |
|---|---|
| Template `left_start` / `right_start` | 0-based positions on the template sequence |
| Genomic product `start` / `end` | **1-based inclusive** |
| `strand` on a priming site | `+` = matches template orientation, `-` = reverse complement |
| BED output | 0-based half-open (converted at the boundary) |
| FASTA / VCF / GFF3 positions | 1-based inclusive |

## Error hierarchy

`primerblast_oss.errors`:

- `PrimerblastError` (base, subclass of `RuntimeError` — legacy `except
  RuntimeError` still works)
  - `ToolMissingError` — external executable not found (raised by design /
    specificity when the binary is absent)
  - `InvalidDatabaseError` — database or index unusable
  - `Primer3Error` — `primer3_core` process failed
  - `BlastError` — `blastn` / `makeblastdb` process failed
  - `MalformedInputError` — malformed input
  - `SearchIncompleteError` — raised by `design_and_screen(strict_search=True)`
    when search evidence is not `complete`
  - `CancelledError` — raised when `cancel_check()` returns True

## Cancellation and progress

Long-running calls accept:

- `progress_callback(stage: str, fraction: float)` — called between pairs of
  `run_pipeline` / `design_and_screen`. Fraction is `0.0`–`1.0`.
- `cancel_check() -> bool` — checked between pairs (coarse-grained); raising
  `CancelledError` aborts the call. SnapyGene runs these calls off the GUI
  thread; the callback may raise from another thread's state without
  deadlock (no locks are held during the check).

## Testing the contract

`tests/test_public_api.py` verifies that results serialize to JSON and
reload, that missing-tool errors are raised, that cancellation works, and
that the wrapper results carry `api_version`. A full end-to-end run against a
real BLAST database is exercised by `benchmarks/` (not part of the unit
suite).
