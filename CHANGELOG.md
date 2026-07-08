# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/).

## [0.2.0] - 2026-07-08

Plant-breeding pipeline: turns the specificity engine into a breeding-oriented
in-silico PCR toolkit addressing the well-known NCBI Primer-BLAST pain points.

### Added
- `assay` subcommand: full pipeline from a **gene (GFF3)**, **interval**, or
  **SNP** — design, multi-reference specificity, SNPs under primers (VCF),
  amplicon conservation, CAPS/dCAPS enzyme scan, and an experimenter
  **risk (low/medium/high)** with reasons.
- `markers` subcommand: evenly spaced markers across a **QTL interval**.
- Local annotation input: **GFF3**, **VCF**, **BED** parsers and a region
  resolver; strand-aware genome fetch via `.fai`.
- **CAPS/dCAPS** module: 40 restriction enzymes, differential digest between two
  alleles, gel-resolvability, enzymes gained/lost.
- Explicit **3'-terminal 5 bp / 10 bp** mismatch counts per priming site.
- **Multi-reference conservation** scoring and **SNP-under-primer** detection
  (with a 3'-end flag).
- Output formats: **CSV**, **BED**, self-contained **HTML report**, oligo
  **order table**, and an **ASCII off-target map**.
- **Provenance manifest** (tool versions, parameters, DB fingerprints) for
  reproducibility.
- `blastn -num_threads` for throughput (design cue from PrimerServer2).

### Fixed
- Anchor-aware reclassification: only the product at the template's own genomic
  locus is "intended"; perfect-match products elsewhere are surfaced as
  **paralog off-targets** (the generic heuristic hid multi-copy genes).
- CAPS assays now target the SNP so every designed amplicon spans it.

## [0.1.0] - 2026-07-08

Initial local Primer-BLAST equivalent.

### Added
- `design`: Primer3 design + BLAST specificity by pairing hits into predicted
  amplicons (F/R, R/F, F/F, R/R) with a 3'-anchored priming model.
- `check`: orientation-free in-silico PCR from pasted primers.
- `tile`: cover a whole region with overlapping amplicons.
- `makedb`: build a BLAST database (with `-parse_seqids`).
- Multi-database screening; gel-resolvability of off-targets; A–D ranking.
