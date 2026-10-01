# Full-primer realignment and site identity

BLAST nominates candidate loci; when a matching genome is provided, each
candidate is verified against the complete primer sequence. If the projected
full-primer footprint is an exact match in that genome, its coordinates are
preserved. This matters in tandem repeats: several adjacent placements can be
equally optimal, and choosing the first fitting-alignment tie can move a known
exact site onto its neighbor. The check uses the actual genome sequence, not
BLAST's reported percentage identity. Partial HSPs are eligible when their
projected full-primer footprint is exact. Indel-bearing or imperfect footprints
continue through the existing fitting aligner and mismatch filters.

Repeated HSP evidence is coalesced only when it describes the same binding
alignment. Identity includes primer, reference, strand, both terminal
coordinates, primer length, alignment strings, mismatch counts and available
thermodynamic annotations. Distinct 3′ endpoints or distinct gapped alignments
remain separate evidence, even when their products have the same outer
coordinates and length. Product coordinates alone are not a deduplication key.
HSP-only sites retain their query/target alignment strings for this purpose.

## Counts and completeness

- `raw_blast_hits` and `raw_hits_per_primer` retain every input HSP, including
  duplicate or malformed rows
- `PrimerHitStats.priming_sites` counts unique accepted binding alignments
- `PrimerHitStats.realigned_sites` and
  `full_length_realignment.accepted_sites` count unique accepted genome-verified
  full-primer alignments
- Binding-site counts and thermodynamic evaluations use those unique sites
- Raw-hit high-copy warnings, malformed-row reporting and unresolved genome
  candidate reporting remain active; site-based high-copy thresholds use the
  unique site count

This prevents repeated HSPs for one site from inflating product counts. It does
not make BLAST exhaustive or eliminate all alignment ambiguity. The existing
imperfect-alignment scoring/tie policy, search-completeness heuristics,
thermodynamic model and default mismatch thresholds are unchanged.

## Synthetic regression

`tests/test_realign_site_identity.py` creates a synthetic reference containing
15 `AC` repeats, a forward primer of ten `AC` repeats, and a downstream reverse
primer. Exhaustive substring matching independently establishes six forward
sites and six products of 130, 132, 134, 136, 138 and 140 bp. Complete exact HSPs
are supplied directly, separating post-BLAST realignment correctness from
BLAST seed/reporting completeness. The reverse-complemented reference is also
tested, alongside partial HSPs, contig edges, insertions/deletions, nearby
imperfect candidates, converging HSPs and distinct alignment evidence.

Run the regression without external tools:

```sh
python -m pytest -q tests/test_realign_site_identity.py
```

For the existing end-to-end synthetic checks, install BLAST+, Primer3 and the
optional `thermo` dependency, set a writable reference cache, and run:

```sh
export BLAST_USAGE_REPORT=false
export SEQWB_REFERENCE_CACHE_DIR="$(mktemp -d)"
python benchmarks/continuous_benchmark.py --max-seconds 30
```

These are correctness regressions, not experimental PCR validation or proof of
whole-genome sensitivity.
