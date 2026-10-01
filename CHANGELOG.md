# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/), and the project aims to follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Native `blast` and `primer3` CLI commands with FASTA file/stdin inputs, text /
  JSON / TSV reports, custom executable paths and file exports. These require no
  GUI server. Primer3 runs without a BLAST database and exports candidate oligos.
- Native `check` supports separate reverse-complement hypotheses, input-sequence
  evidence, TSV output, optional reference-product FASTA export and per-reference
  GFF3 gene annotation. Existing profiles and thermodynamic options are preserved;
  sequence extraction runs only when requested. CLI streams/files use UTF-8.
- PCR-check evidence distinguishes amplification with the supplied oligos from
  reverse-complement alternatives, shows actual F/R extension directions and
  separates label swaps and single-primer products. Binding-site directions are
  available even without products; incomplete searches and input-error diagnosis
  remain unresolved. Literal and changed inputs have separate result sections;
  comparison tables show input, reverse, complement and reverse complement.
  Site tables include up to 200 candidates with full counts.
- Standalone local BLAST nucleotide search and Primer3 GUI tabs, also available
  as `blast` / `primer3` operations through the shared agent API. BLAST reports
  both-strand HSP alignments without PCR filters; Primer3 requires no database
  and labels specificity as not evaluated, with candidate transfer to PCR check.
- Explicit automatic reverse-complement hypotheses for up to two copied
  genomic primer inputs, preserving ordered 5'-3' oligos and separate searches.
- Predicted-product reference FASTA downloads and overlapping GFF3 gene names
  in PCR check. Extraction uses an associated indexed FASTA or `blastdbcmd`;
  missing sequence/annotation remains unresolved. Existing design, sequencing
  and assay maps also export their reference amplification span.
- Allele-specific PCR / ARMS design for SNP assays: ref/alt primers terminate
  at the SNP, optional deliberate -2/-3 mismatches are enumerated, best
  classical AS-PCR and tetra-primer ARMS sets are generated, AA/AB/BB gel
  patterns are scored, and top sets are re-screened against local BLAST
  databases. Allele discrimination and genome-wide off-target evidence remain
  separate outputs.
- Gel-aware CAPS/dCAPS ranking based on digest fragment sizes, log-size
  separation, ladder coverage, agarose size ranges, heterozygote band patterns,
  and deterministic M/AA/AB/BB virtual-gel SVG output. Built-in 100 bp and 1 kb
  ladder presets can be replaced by user-supplied band sizes.
- Target-aware `sequence` workflow for Sanger/amplicon resequencing. It
  reuses the tiling engine, accepts gene/interval/FASTA inputs, exposes
  configurable amplicon length and overlap, reports coverage gaps, maps
  minus-strand targets back to genomic coordinates, and can add M13 tails to
  order oligos without changing specificity evaluation.

### Changed / Fixed
- Primer-dimer / hairpin analysis now runs in **`design` and `tile`** too (not
  only `assay` / `multiplex`), affecting rank, and dimer ΔG/Tm are shown in the
  text, TSV and CSV outputs.
- The **high-copy-primer warning** now triggers on the number of *priming sites*
  (`--high-copy-site-threshold`, default 500) instead of raw BLAST HSP count, so
  a specific primer on a large genome is no longer falsely flagged
  repeat-sensitive.

### Added
- **`multiplex-design` subcommand**: designs primers for several targets (a
  multi-record template FASTA) and picks **one mutually-compatible pair per
  target** so no forward/reverse primer forms a concerning cross-dimer with a
  primer of another target (backtracking selector + greedy partial fallback,
  `--candidates-per-target`, `--require-specific`). NCBI Primer-BLAST designs
  each amplicon independently and does not do this.
- **Primer-dimer / hairpin / multiplex analysis** (`primer3-py`): per-pair
  hairpin, self-dimer and forward×reverse cross-dimer scoring (ΔG + Tm), folded
  into risk. New `multiplex` subcommand checks every-primer-vs-every-primer
  cross-dimers across a pool to pick compatible sets — something NCBI
  Primer-BLAST does not do.
- **NCBI Primer-BLAST head-to-head benchmark** (Arabidopsis TAIR10, published):
  on the same template, primerblast-oss's top de-novo pair and its specificity
  verdict match the live NCBI web service exactly (`benchmarks/RESULTS.md` §8),
  extended to a **six-locus panel** (§8b) run against the live NCBI service and
  both local tools — primerblast-oss stays within the NCBI / PrimerServer2 range
  on all six, matching NCBI in rejecting a non-3'-anchored off-target that
  PrimerServer2 keeps.
- **Automated multi-locus PrimerServer2 benchmark** (`benchmarks/head_to_head_ps2.py`):
  designs a pair in each of N genome windows and compares both tools' predicted
  amplicons under matched parameters. Across 40 TAIR10 loci the tools agree on the
  exact amplicon set for 92% of non-repetitive loci; residual differences are the
  3'-anchor-vs-Tm-window model difference and repeat-copy enumeration, none an
  implementation error (`benchmarks/RESULTS.md` §9).
- **Optional thermodynamic scoring** (`primer3-py`): each priming site gets a
  duplex Tm and 3'-end stability ΔG; thermodynamically non-viable sites are gated
  out of amplicon prediction. Enabled via `--genome-fasta` (automatic in
  `assay`), with `--min-anneal-tm` / `--max-3p-dg` / `--no-thermo` /
  `--no-thermo-gate`. Falls back to the mismatch/3'-anchor model when primer3-py
  is absent. Tm/ΔG are reported per product, and because non-viable sites are
  gated, off-target counts and risk reflect only thermodynamically-plausible
  products.

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
