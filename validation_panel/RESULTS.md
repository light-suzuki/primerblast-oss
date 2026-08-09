# Wet-lab validation panel — results

**Status: panel defined, no assessed assays yet.**

`validation_panel/manifest.example.json` defines the versioned manifest
schema (inputs pinned by `sequence_sha256`, protocol metadata, observed
outcomes, exclusion provenance) and the scoring protocol is documented in
[`validation_panel/README.md`](README.md). Metrics cannot be reported until
real wet-lab outcomes are recorded; until then every cell below reads
`n/a` by design.

## Planned report

| Metric | Value | n |
|---|---|---|
| Off-target detection sensitivity | n/a | — |
| Off-target false-positive rate | n/a | — |
| Intended-product amplification success | n/a | — |
| Rank calibration (A/B/C/D/I) | n/a | — |
| Risk calibration (low/medium/high) | n/a | — |
| CAPS digest concordance | n/a | — |

This table is produced by:

```bash
python benchmarks/validation_analysis.py \
  --manifest validation_panel/manifest.example.json
```

Once assays are recorded, the numbers above are filled in and the
`RESULTS.md` link from the README points here.
