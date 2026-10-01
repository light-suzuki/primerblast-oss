# primerblast-oss

[![Release](https://img.shields.io/github/v/release/light-suzuki/primerblast-oss?sort=semver)](https://github.com/light-suzuki/primerblast-oss/releases)
[![CI](https://github.com/light-suzuki/primerblast-oss/actions/workflows/ci.yml/badge.svg)](https://github.com/light-suzuki/primerblast-oss/actions/workflows/ci.yml)
[![Benchmark](https://github.com/light-suzuki/primerblast-oss/actions/workflows/benchmark.yml/badge.svg)](https://github.com/light-suzuki/primerblast-oss/actions/workflows/benchmark.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)

**English** | [日本語](README.ja.md)

This repository is the PCR engine. [Sequence Workbench](https://github.com/light-suzuki/Gene-research)
embeds it as part of broader breeding and gene-research workflows. See
[the embedding boundary](docs/EMBEDDING.md). The standalone GUI defaults to PCR
tools; legacy research views remain available through the legacy-tools switch.

The GUI (`python -m primerblast_oss.webapp`) includes standalone **BLAST** and
**Primer3** tabs. `blastn`, `makeblastdb` and `primer3_core` are required external
tools. BLAST accepts DNA or multi-record FASTA, searches both strands and exports
HSP alignments as TSV. Primer3 designs without a database, labels specificity as
not evaluated, and transfers candidate oligos to PCR check.

PCR check accepts ordered 5'-3' oligos without F/R swapping or strand settings.
For sequences copied from a genomic display, select the copied-sequence input
mode to test original/reverse-complement combinations separately (up to two
primers). Review the candidate oligos; products from separate hypotheses are not
products of a single mixed PCR reaction. Results distinguish literal-input
products, F/R label swaps (R extends right and F extends left), single-primer
products, and alternative oligos requiring reverse-complement changes. Each
product shows input versus candidate sequences, extension directions, coordinates
and mismatches. A comparison table displays input, reverse, complement (3'-5')
and reverse complement (5'-3'); supplied inputs are never overwritten. Literal
results and changed-oligo hypotheses appear in separate sections. Binding-site
tables show directions even without products (up to
200 sites, with full counts). A missing product under an incomplete search is
unresolved; search cannot establish a typo or guarantee experimental amplification.
Predicted products export reference
FASTA from the associated indexed genome or `blastdbcmd`, and display overlapping
gene names/IDs from a matching GFF3. Reference FASTA is the genomic plus strand,
without primer mismatches or 5' tails incorporated. The FASTA/GFF3 fields apply to
the first selected database; other references use `db_genomes` and `db_gff3`
mappings in the shared API. Missing extraction or annotation is unresolved.
Design, sequencing and assay maps also export their reference amplification span.

Standalone BLAST and Primer3 are GUI tools. The native and JSON agent CLIs
provide integrated PCR workflows.

A local, open-source, Primer-BLAST-style **command-line tool** for plant breeding
and genetics. It designs PCR primers with **Primer3** and verifies their
**specificity** entirely offline against local BLAST+ databases — including
unpublished genomes and several cultivars at once — and adds in-silico PCR,
whole-region tiling, SNP-under-primer detection, amplicon conservation analysis,
CAPS/dCAPS and allele-specific PCR (AS-PCR / tetra-ARMS) marker design, and an experimenter-facing risk score.

The browser's **Enzymes / 制限酵素** tab includes an offline 1,088-name catalog,
related enzyme and verified product names, and both-strand cleavage diagrams.
See [catalog coverage, provenance and prediction limits](docs/restriction-catalog.md).

> The core is pure Python (standard library only) and calls out to `primer3_core`
> and BLAST+. The unit tests require **no external tools or data**.

## Why local and open source

NCBI Primer-BLAST is excellent, but it is **not open source** and exists only as a
hosted web service. You therefore cannot audit, fork, or self-host it, cannot run
it next to your data, and cannot submit unpublished or embargoed genomes to it.

Depending on an external service also binds a pipeline to that service's
availability and policies — rate limits, maintenance windows, occasional outages.
Running the workflow **locally and offline** removes that dependency and makes a
run fully reproducible: pinned FASTA, GFF3, VCF, BLAST database, and tool
versions, with no queue and no login. This matters most in precisely the case the
hosted tool cannot serve: local, unpublished, multi-cultivar genomes.

primerblast-oss is **MIT-licensed**, so anyone can read, run, and build on it.

## Why this exists

NCBI Primer-BLAST performs two steps: (1) Primer3 designs candidate primers, and
(2) BLAST screens each primer against a database and **pairs the hits into
predicted amplicons**, flagging unintended PCR products. Most local "primer +
BLAST" scripts perform only the per-primer BLAST and omit step (2) — the step
that actually detects off-target products.

`primerblast-oss` is an independent implementation of step (2) — using a
BLAST-alignment-based priming model, not NCBI's exact algorithm — with a focus
on things a local, breeding-oriented workflow needs:

| | NCBI Primer-BLAST | PrimerServer2 | primerblast-oss |
|---|---|---|---|
| Primer3 design | ✅ | ✅ | ✅ |
| Pairs BLAST hits into predicted amplicons | ✅ | ✅ | ✅ |
| Off-target products from either primer as F/F, R/R, F/R | ✅ | ✅ | ✅ |
| 3'-end-aware priming model | ✅ | ✅ (`--use-3-end`) | ✅ |
| Runs offline on **unpublished / local** genomes | hard | ✅ | ✅ |
| Thermodynamic off-target scoring (Tm-based) | ✅ | ✅ (core model) | optional (primer3-py) |
| **Multiplex** primer-dimer *checking* of a pool | — | ✅ | ✅ |
| **Multiplex** compatible-set *design* (one pair/target) | — | — | ✅ |
| Screen against **multiple databases** in one run | — | partial | ✅ |
| **In-silico PCR** from pasted primers (orientation-free) | — | ✅ | ✅ |
| **Tile a whole region** with overlapping amplicons | — | — | ✅ |
| Gel-resolvability of off-targets (size-gap aware) | — | — | ✅ |
| Breeding assay: GFF3/VCF/CAPS/QTL + risk | — | — | ✅ |
| Scriptable CLI + library, no queue/login | limited | ✅ | ✅ |
| Curated, continuously-updated databases | ✅ | — | — |
| Mature hosted web server | ✅ | ✅ | local GUI only |

This is about **fit for a local, offline workflow**, not a claim of being better
overall. NCBI Primer-BLAST has real advantages this tool does not: curated,
continuously updated databases, a mature thermodynamic model, and deeper
primer-dimer / hairpin analysis. [PrimerServer2](https://github.com/billzt/PrimerServer2)
is also a strong local tool and shares much of the core recipe; primerblast-oss's
additions over it are whole-region tiling, gel-resolvability, one-run multi-database
screening, the breeding assay (GFF3/VCF/CAPS/QTL/risk), and multiplex-set *design*
rather than only dimer *checking*. A ✅ in more than one column means the capability
exists on each side — not that the underlying models are identical or that outputs
will match.

**Benchmarks (summary):** across **40 randomly-placed Arabidopsis TAIR10 loci**,
primerblast-oss and PrimerServer2 predict the same amplicon set (count, size,
coordinates) on **92 % of non-repetitive loci**; three hand-checked *Lotus
japonicus* pairs matched PrimerServer2 exactly; and across **six loci** run
against the **live NCBI Primer-BLAST** service, primerblast-oss stays within the
NCBI / PrimerServer2 range on every one — matching NCBI in rejecting a
non-3'-anchored off-target that PrimerServer2 keeps. Full method, numbers, and an
analysis of the residual disagreements are in
[`benchmarks/RESULTS.md`](benchmarks/RESULTS.md) §7–§9.

## Validation status

The short-term goal is to be a **superset of PrimerServer2** for local,
scriptable primer work: matching its specificity behaviour on local genomes while
adding multi-database screening, tiling, marker design, breeding-assay outputs,
and offline reproducibility. The evidence to date:

- **PrimerServer2, 40-locus automated head-to-head (Arabidopsis TAIR10).** With
  matched parameters, the two tools agree on the exact predicted amplicon set for
  **33 / 36 (92 %) non-repetitive loci**. Every residual disagreement is
  accounted for: repetitive loci where both tools call the primer non-specific but
  enumerate repeat copies differently, and marginal sites where a primer's 3' end
  is not fully aligned — which primerblast-oss rejects as non-priming and
  PrimerServer2 keeps on duplex Tm. None trace to an implementation error
  ([`benchmarks/RESULTS.md`](benchmarks/RESULTS.md) §9).
- **PrimerServer2, *Lotus japonicus*.** Three hand-checked pairs matched exactly
  on amplicon count, size, and coordinates (§7).
- **NCBI Primer-BLAST (six loci).** Six pairs were run through the live NCBI
  service and both local tools. All three agree exactly on the three clean loci;
  on the three borderline loci primerblast-oss sits **within the NCBI /
  PrimerServer2 range** — on one it matches NCBI in rejecting an off-target that
  PrimerServer2 keeps (a non-3'-anchored site), on another it matches
  PrimerServer2. Every difference is a single borderline product near a Tm or
  3'-alignment threshold, not an error (§8b). Still a modest sample, and NCBI
  screens its own *Arabidopsis* assembly rather than the local FASTA.
- **Continuous regression benchmark.** A separate scheduled GitHub Actions workflow
  builds a synthetic FASTA/BLAST database and exercises Primer3 design, BLAST
  amplicon pairing, duplicate/off-target classification, thermodynamic gating,
  and multiplex dimer checks weekly (and on manual dispatch).
- **Wet-lab validation panel.** A versioned panel
  ([`validation_panel/`](validation_panel/README.md)) records prospective
  experimental outcomes (protocol metadata, observed bands, exclusion
  provenance) against immutable in-silico predictions, and
  [`benchmarks/validation_analysis.py`](benchmarks/validation_analysis.py)
  reports off-target sensitivity / false-positive rate, amplification failure,
  rank/risk calibration, and CAPS digest concordance. Until real assays are
  recorded, the panel is defined but uncalibrated
  ([`validation_panel/RESULTS.md`](validation_panel/RESULTS.md)).

Against **NCBI Primer-BLAST**, no claim of drop-in equivalence is made: NCBI
retains the advantage in curated, continuously-updated databases, hosted UX, and a
private, long-matured specificity model. primerblast-oss is the stronger choice
when the data that matter are local, unpublished, multi-reference, or must run
reproducibly in scripts.

## How specificity is judged

For each primer pair and each database:

1. `blastn -task blastn-short` finds every near-full-length hit of each primer.
2. A hit becomes a **priming site** only if the primer's **3' end is aligned**,
   its **3'-terminal base matches**, mismatches within the 3' window
   (`--three-prime-window`, default 5) are `≤ --max-3prime-mismatch` (default 1),
   and total mismatches over the full primer are `≤ --max-total-mismatch`
   (default 4). Unaligned 5' bases count as mismatches.
3. On each subject, a plus-strand priming site is paired with every downstream
   minus-strand priming site within the product-size window
   (`--min-product`..`--max-product`). The product span is measured 5'→5'
   (true PCR amplicon length).
4. A product is **on-target** if both primers anneal perfectly (0 mismatch),
   in F/R orientation, at the designed size (± `--size-tolerance`). Everything
   else is **off-target**. A pair is **specific** when exactly one product (the
   intended one) is predicted in every screened database.

Pairs are scored and ranked A–D by specificity, Tm balance, GC, and 3'-dimer
strength.

Genome-verified exact loci retain their nominated coordinates. Repeated HSPs
for the same binding alignment count once; distinct 3' endpoints or gapped
alignments remain separate evidence. See [site identity and count semantics](
docs/priming-site-identity.md).

Use `--specificity-profile ncbi` to switch the mismatch thresholds to a
NCBI-Primer-BLAST-like stringency profile: up to 5 total mismatches are kept as
candidate priming sites, up to 1 mismatch is allowed within the 3'-terminal 5 bp,
and a terminal-base mismatch is counted rather than rejected outright. This is a
compatibility profile for threshold behavior; it is still not NCBI's private
algorithm or database.

## Requirements

- Python ≥ 3.8 (standard library only)
- `primer3_core` (Debian/Ubuntu: `apt install primer3`)
- BLAST+ `blastn` / `makeblastdb` (`apt install ncbi-blast+`)
- A nucleotide BLAST database (see below)

## Install

```bash
pip install -e .              # provides the `primerblast-oss` command
pip install -e '.[thermo]'    # + optional primer3-py (thermodynamics & dimers)
# or run without installing:
python -m primerblast_oss --help
```

Both `primerblast-oss <subcommand>` (after install) and
`python -m primerblast_oss <subcommand>` are equivalent; this README uses the
`python -m` form so the examples work without installing.

## Quick start

```bash
# 1. build a BLAST database from a genome FASTA (once)
python -m primerblast_oss makedb genome.fa --out-db mydb

# 2. design primers on a template and check them against that genome
python -m primerblast_oss design \
  --template-fasta my_gene.fa --product-size 150-500 --db mydb

# 3. or just in-silico-PCR a pair you already have
python -m primerblast_oss check \
  --forward GACAAGGAATCAGCGGCTCT --reverse GCAGCGTTTTGTAGTGGGTG --db mydb
```

A local browser GUI wrapping these subcommands also exists, but it is an
optional extra — see [Web GUI (optional)](#web-gui-optional) near the end.

## Usage

primerblast-oss is a **CLI tool**. Subcommands: **design**, **check**,
**multiplex**, **multiplex-design**, **tile**, **sequence**, **assay**,
**markers**, **makedb**. Run `python -m primerblast_oss <subcommand> --help` for the full
option list of any one.

The native CLI focuses on integrated PCR workflows: design plus specificity
screening, primer checks, multiplex compatibility and marker assays.
Standalone BLAST and Primer3 remain in the GUI.

In PCR check, automatic reverse-complement alternatives are OFF by default.
Choose ON/OFF in the GUI, or use native `check --no-auto-orientation`
(`--input-orientation as_supplied`) to explicitly turn them OFF.
Use `--input-orientation auto` for ON. Both modes search both strands and all
F/R, R/F, F/F and R/R combinations under the selected search conditions.
OFF does not remove predictions for the supplied oligos.

`multiplex` checks primer-dimer compatibility across a pool of primers (needs
`primer3-py`) — every primer against every other, to pick sets you can run
together:

```bash
python -m primerblast_oss multiplex \
  --primer A_F=... --primer A_R=... --primer B_F=... --primer B_R=...
```

`multiplex-design` goes a step further: give it several targets (a multi-record
template FASTA) and it designs candidates for each, then picks **one
mutually-compatible pair per target** so no two primers form a concerning
cross-dimer. NCBI Primer-BLAST designs each amplicon independently and cannot do
this.

```bash
python -m primerblast_oss multiplex-design \
  --template-fasta targets.fa --db $DB --genome-fasta genome.fa \
  --product-size 80-300 --candidates-per-target 5 --require-specific
```

### `design` — region + product size → primer pairs

```bash
python -m primerblast_oss design \
  --template-fasta my_gene.fa \
  --db /path/to/genome_db \
  --product-size 150-500 --format text
```

Screen against several cultivar genomes at once (specific in *all* of them):

```bash
python -m primerblast_oss design \
  --template "ACGT..." --template-id MyLocus \
  --db /data/blastdb/cultivarA --db /data/blastdb/cultivarB \
  --product-size 200-800 --format tsv
```

### `check` — primer sequences → all predicted PCR products (in-silico PCR)

Paste primers; orientation is **not** constrained (any primer may act as
forward or reverse). Every product is listed with its size and the size gap to
the nearest other product, so you can judge whether extra bands are resolvable.

`--input-orientation auto` separately tests reverse-complement alternatives for
up to two primers; the default `as_supplied` preserves the entered oligos. Text
and JSON distinguish literal-input products, label swaps and changed-oligo
candidates. `--show-sequence-forms` includes reverse/complement forms in text;
JSON always contains them. `--primers-fasta -` reads a pool from stdin.

```bash
primerblast-oss check --primers-fasta pair.fa --db mydb \
  --input-orientation auto --genome-fasta genome.fa --gff3 genes.gff3 \
  --products-fasta products.fa --format json --out check.json
```

The product FASTA is reference plus-strand sequence, with separate IDs for DBs,
result groups and products. It does not substitute mismatched oligo bases or
5' tails; incomplete extraction is reported on stderr and in JSON. Extraction
uses the associated indexed FASTA or `blastdbcmd` and runs only for
`--products-fasta`. `--gff3` annotates the first DB; repeated `--db-gff3 DB=GFF3`
maps other references. Annotation alone does not require sequence extraction.
`--format tsv` produces product rows with DB, changed inputs, completeness,
participating oligos/directions and genes; zero-product groups and search state
are reported on stderr. JSON retains all groups, binding-site samples and
unresolved states. Existing specificity profiles and thermodynamic options apply.

```bash
python -m primerblast_oss check \
  --forward GCACTCTAGAGGTTCAAGGCC --reverse TGGTACGTGTGGTTCAGTTTCA \
  --db /path/to/genome_db
# or a pool of primers by name:
python -m primerblast_oss check \
  --primer F1=ACGT... --primer F2=TTGC... --primer R1=GGCA... \
  --db /path/to/genome_db --format json
```

### `tile` — region + amplicon length → overlapping amplicons covering the whole region

Primer-BLAST designs one amplicon around a target; `tile` instead walks the
entire region with overlapping amplicons (e.g. to sequence a whole gene).

```bash
python -m primerblast_oss tile \
  --template-fasta gene.fa \
  --amplicon-min 400 --amplicon-max 700 --overlap 60 \
  --db /path/to/genome_db
```

### `sequence` — overlapping amplicons for Sanger / resequencing

`sequence` is a target-aware wrapper around the tiling engine. Give a gene,
genomic interval, raw sequence, or FASTA and choose both the desired amplicon
length and overlap between neighboring products.

```bash
python -m primerblast_oss sequence \
  --gene Psat.cameor.v2.1g00050 --gff3 genome.gff3 --genome genome.fa \
  --db $DB/cameor_v2 --amplicon-size 600-800 --overlap 120
```

For genomic targets, `--flank` adds sequence outside the requested region for
primer placement while coverage is still scored against the target itself.
Coverage gaps are reported explicitly. Optional `--m13-tails` prepends the
universal M13 sequences to the **order oligos only**; genome specificity is
evaluated on the annealing portion of each primer.

```bash
python -m primerblast_oss sequence \
  --interval chr1:100000-104000 --genome genome.fa --db $DB/genome \
  --amplicon-size 500-700 --overlap 100 --m13-tails --format order
```

### `assay` — full breeding assay from a gene / interval / SNP

Resolves a target from a **local genome + GFF3/VCF**, designs primers, checks
specificity across **several reference genomes**, flags **SNPs under primers**,
scores **amplicon conservation**, runs an optional **CAPS/dCAPS** enzyme scan,
and assigns an experimenter **risk (low/medium/high)**. Outputs text, JSON,
CSV, BED, an oligo **order table**, or a self-contained **HTML** report.

```bash
DB=/path/to/blastdb
# a gene, screened across three cultivars, with a VCF and CDS feature
python -m primerblast_oss assay \
  --gene Psat.cameor.v2.1g00050 --gene-feature cds --gff3 genome.gff3 \
  --genome genome.fa \
  --db $DB/cameor_v2 --db $DB/unpublished_cultivar --db $DB/ZW6 \
  --vcf variants.vcf --flank 100 --product-size 150-600 --format html --out report.html

# a CAPS/dCAPS marker spanning a SNP (alt allele given)
# auto-rank digest patterns by actual fragment sizes and write an M/AA/AB/BB gel
python -m primerblast_oss assay --snp chr1:6385 --alt A \
  --genome genome.fa --db $DB/cameor_v2 --flank 250 \
  --ladder auto --gel-percent auto --virtual-gel marker.svg --format text

# force a 100-bp ladder / 2% agarose model
python -m primerblast_oss assay --snp chr1:6385 --alt A \
  --genome genome.fa --db $DB/cameor_v2 \
  --ladder 100bp --gel-percent 2.0 --format json
```

For CAPS/dCAPS, biological specificity and experimental readability are
reported separately. A pair can therefore be non-specific at the sequence
level but still receive
`marker_verdict=gel_scorable_with_separated_offtargets` when predicted
off-target products remain clearly separated from the AA/AB/BB diagnostic band
patterns in every screened database. Off-target bands are overlaid per
database, not pooled across references.

Very short or very long diagnostic fragments are **not automatic rejection
criteria**. They remain in the candidate set and are reported as soft warnings
when they fall outside the selected ladder span or recommended agarose range.
This lets the experimenter keep an otherwise highly discriminating marker and
adjust the gel/ladder conditions manually.

### Allele-specific PCR / ARMS from the same SNP assay

When `assay --snp ... --alt ...` is used, the SNP workflow now also builds
allele-specific PCR candidates in parallel with CAPS/dCAPS. Ref- and alt-specific
primers place the SNP at the **3' terminal base**. Candidate deliberate
destabilizing mismatches at the **-2 and -3 positions from the 3' end** are
enumerated and ranked by the difference between intended- and non-target-allele
terminal mismatch patterns.

The same parent pair is also used to build **tetra-primer ARMS-PCR** candidates:
the two outer primers form a control amplicon, while two allele-specific inner
primers generate different-size ref and alt bands. The expected AA / AB / BB
patterns are checked for gel separation.

```bash
python -m primerblast_oss assay \
  --snp chr1:6385 --alt A \
  --genome genome.fa --db $DB/cameor_v2 \
  --aspcr-candidates-per-allele 4 \
  --aspcr-tetra-candidates 6 \
  --aspcr-pairs-to-screen 2 \
  --format json
```

For the best ref-AS, alt-AS, and tetra-ARMS sets, the tool re-runs genome-wide
in-silico PCR against the selected databases. The allele-discrimination score
  and genome off-target screen are reported separately.

Automatic AS-PCR/tetra-ARMS preference requires an anchored, complete genome
screen without unexpected products and an evaluated primer-structure check
without concerns. Size alone does not establish an intended product. Other
assemblies without mapped target anchors remain `unverified_size_only`, and
unscreened or unresolved designs remain candidates rather than preferred modes.

**Important:** a 3'-terminal mismatch does not guarantee complete allele
rejection. ARMS designs commonly add a deliberate near-3' mismatch to increase
discrimination, but polymerase, annealing conditions, and the exact mismatch
combination matter. These scores rank designs; they are not probabilities of
successful genotyping and should be validated experimentally.

### `markers` — evenly spaced markers across a QTL interval

```bash
python -m primerblast_oss markers --interval chr1:80000000-90000000 \
  --genome genome.fa --db $DB/cameor_v2 --n-markers 20 --format json
```

### `makedb` — build a database (with `-parse_seqids`)

```bash
python -m primerblast_oss makedb genome.fa --out-db genome_db
```

## How it approaches NCBI Primer-BLAST's common pain points

These are the pain points the tool is designed around. Coverage varies and some
items are partial — see [Limitations](#limitations).

| # | Common pain point | primerblast-oss approach |
|---|---|---|
| 1 | weak at batch / many regions | `markers`, `assay` over BED/gene lists; CLI + library, scriptable |
| 2 | poor fit for local / custom assemblies | everything runs on local FASTA + BLAST DB; `makedb` helper |
| 3 | weak multi-reference comparison | `--db` repeatable; amplicon **conservation** scored per reference |
| 4 | whole-chromosome design is clumsy | `tile` + `markers` generate primers across whole regions/intervals |
| 5 | primer strand/orientation unclear | every binding site reports strand, 5'/3' coords, extension direction |
| 6 | unexpected side products hard to read | BLAST hits **paired into predicted amplicons** with sizes |
| 7 | F-F / R-R products hard to see | enumerated explicitly and shown in the ASCII map & tables |
| 8 | 3'-end mismatch not visible | explicit 3'-terminal **5 bp / 10 bp** mismatch counts per hit |
| 9 | paralogs / repeats / duplications | genome-wide pairing surfaces duplicated priming sites (no dedicated repeat mask; bounded by BLAST `-max_target_seqs`) |
| 10 | not built for CAPS/dCAPS | `caps` scan: enzymes that digest two alleles differently, gel gap |
| 11 | weak GFF3 / VCF / QTL integration | `--gene`/`--gff3`, `--vcf`, `--interval`, BED input |
| 12 | opaque empty results | Primer3 explain string surfaced; per-stage diagnostics |
| 13 | weak reproducibility | provenance manifest pins tool versions, params, DB fingerprints |
| 14 | weak experimenter-facing scoring | `risk` rolls up every signal into low/medium/high with reasons |
| 15 | side products not visualized | ASCII off-target map + BED track for a genome browser |

Design cues were taken from
[PrimerServer2](https://github.com/billzt/PrimerServer2) (strand-aware BLAST-hit
pairing, multi-threaded `blastn`, coordinate input) and NCBI Primer-BLAST.

## Limitations

Honest scope, so you know what it does *not* do:

- **Only spot-checked against NCBI Primer-BLAST.** It is an independent
  implementation. One published-Arabidopsis locus matched NCBI exactly (same top
  pair, same specificity verdict — `benchmarks/RESULTS.md` §8), but that is a
  single locus, not systematic validation; results elsewhere are plausible and
  internally checked, not guaranteed to match NCBI's output.
- **Thermodynamic scoring is optional** (needs `pip install primer3-py`). When a
  `--genome-fasta` is supplied (automatic in `assay`), each site gets a duplex Tm
  and 3'-end ΔG via primer3 and thermodynamically non-viable sites are gated out;
  without it, priming falls back to the BLAST-alignment mismatch / 3'-anchor rule.
- **Off-target discovery is bounded by BLAST `-max_target_seqs`** (default 5000).
  In extremely repetitive regions some hits can be missed; there is no dedicated
  repeat mask.
- **dCAPS support is best-effort** and CAPS calls depend on the enzyme table
  (~40 common enzymes), not an exhaustive REBASE set.
- **Batch/QTL modes work but are not benchmarked at large scale**; each pair
  costs a BLAST search, so wide sweeps are IO/CPU bound.
- **primer-dimer / hairpin analysis needs primer3-py** (optional). With it, each
  pair gets hairpin / self-dimer / cross-dimer scoring and the `multiplex`
  subcommand checks a whole pool; without it, only Primer3's design-time limits
  apply.

Contributions toward any of these are welcome.

### Key options (shared by design/check/tile)

| Option | Meaning | Default |
|---|---|---|
| `--product-size` (design) | one or more ranges, e.g. `150-500,500-1000` | `70-1000` |
| `--amplicon-min/--amplicon-max/--overlap` (tile) | tiling geometry | 400/800/40 |
| `--opt-tm/--min-tm/--max-tm` | primer melting temperature window | 60/57/63 |
| `--specificity-profile` | mismatch preset: `local-strict` or `ncbi` | `local-strict` |
| `--max-total-mismatch` | mismatches allowed for an off-target to still prime | 4 |
| `--max-3prime-mismatch` | mismatches allowed in the 3' window | 1 |
| `--three-prime-window` | size of the 3' window | 5 |
| `--max-product` | largest off-target amplicon considered amplifiable | 4000 |
| `--max-target-seqs` | BLAST hit cap; raise for repetitive genomes | 5000 |
| `--exhaustive` | convenience mode using a higher BLAST hit cap | off |
| `--num-threads` | `blastn` worker threads | 4 |
| `--high-copy-hit-threshold` | raw BLAST HSP count that triggers a repeat-sensitivity warning | 10000 |
| `--gel-min-gap` | size gap (bp) to resolve two products on a gel | 50 |
| `--no-3prime-terminal` | allow a mismatched 3' terminal base | off |
| `--genome-fasta` | enable primer3-py thermodynamic site scoring for design/check/tile | off |
| `--min-anneal-tm` / `--max-3p-dg` | thermodynamic off-target gate thresholds | 40 / -5 |
| `--dimer-dg-warn` / `--dimer-tm-warn` | primer-dimer/hairpin warning thresholds for assay/multiplex | -8 / 45 |
| `--format` | `text` \| `json` \| `tsv` (design only) | text |

### Verdicts & scoring

Design/tile rank pairs **A–D**: **A** = single product in every database;
**B** = extra products but all far enough in size to resolve on a gel;
**C/D** = one or more off-targets that **co-migrate** with the intended band.
Co-migrating off-targets are penalized heavily; size-resolvable ones only
lightly — matching how such pairs are used in practice.

## JSON output (scripting-friendly)

`--format json` emits a stable schema for piping into other tools (or a GUI).
For native `design` and `tile`, a single template retains the flat result;
multiple templates produce one document with `mode` and a `templates` array.
PCR-check binding-site tables retain thermodynamically rejected candidates;
these are excluded from product predictions, with evaluation state reported
separately. Failed thermodynamic calculations remain unresolved.
Every object carries what a consumer needs without recomputation:

- **primer**: `forward`/`reverse`, `tm_f`/`tm_r`, `gc_f`/`gc_r`, `left_start`/
  `right_start` (0-based template positions), `product_size`, `penalty`.
- **product** (design/check): `subject`, `start`, `end`, `size`, `orientation`
  (e.g. `F/R`, `R/R`), `fwd_mismatch`+`rev_mismatch`, `nearest_gap` (bp to the
  closest other product — drives gel-resolvability shading).
- **pair verdict**: `rank`, `score`, `specific_all_db`, `gel_distinguishable`,
  `total_on_target`/`total_off_target`/`total_comigrating`, and `per_db`
  breakdown.
- **tile**: `index`, `covers` `[start,end]`, `gap_to_prev` (overlap>0 / gap<0),
  plus the full pair object — enough to draw the amplicon track over the region.

## Library API

```python
from primerblast_oss import run_pipeline, in_silico_pcr, design_tiling

# design + specificity
result = run_pipeline("MyLocus", template_seq, ["/data/db/genome"])
for pair in result.pairs:
    print(pair.forward, pair.reverse, pair.specificity["rank"])

# in-silico PCR from arbitrary primers
res = in_silico_pcr({"F": "ACGT...", "R": "TTGC..."}, "/data/db/genome")
for a in res["products"]:
    print(a.size, a.subject, a.orientation)

# tile a whole region
tiles = design_tiling("gene", template_seq, ["/data/db/genome"],
                      amplicon_min=400, amplicon_max=700, overlap=60)
```

## Web GUI (optional)

The CLI is the primary interface. As a convenience, a **local browser front end**
wraps the main design workflows — no cloud, no third-party Python dependencies (it is built
on the standard-library `http.server`). It is an extra, not the main tool. Run it
on the machine where `primer3_core`, `blastn`, and your BLAST databases live (e.g.
inside WSL):

```bash
python -m primerblast_oss.webapp             # serves http://127.0.0.1:8799, opens a browser
python -m primerblast_oss.webapp --port 9000 --no-browser
```

It exposes tabs for design / in-silico PCR / tiling / sequencing / assay / QTL markers / build
DB, auto-discovers BLAST databases under `~/.codex/blast_databases`,
`~/blast_databases` and `./databases`, runs each job in a background thread, and
offers one-click TSV/CSV/BED/JSON downloads. English / 日本語 toggle in the header.
It binds to loopback (`127.0.0.1`) only; on WSL2 Windows browsers reach it via the
default localhost forwarding. Integrated PCR workflows are also available from
the CLI; standalone BLAST and Primer3 are GUI tools.

The purpose-first start page, beginner/advanced controls, inline input guidance,
and responsive layout work in both languages. Load a small FASTA/text file into
the sequence editor, or enter server-side genome/GFF3 paths. Optional `DB=FASTA`
associations enable full-length verification per database; other assemblies are
never silently assigned the design-reference FASTA. Sequencing results include
PCR coverage gaps and optional M13 ordering sequences. Assay results expose
CAPS/dCAPS bands, virtual gels, and AS-PCR/tetra-ARMS oligos with raw evidence in
expandable details. Unfinished searches remain unresolved. All results are
computational candidates, rather than wet-validated assays.

Choose a reference first to fill FASTA, GFF3 and search database paths. The GUI
recognizes an indexed FASTA and a uniquely named sibling GFF3 beside a database;
it does not guess between multiple annotations. Explicit profiles can be listed
in the server-side `~/.codex/primerblast-oss/references.json`, or in the JSON file
named by `PRIMERBLAST_REFERENCES`:

```json
[{"name":"Reference name","genome":"/data/genome.fa","gff3":"/data/genes.gff3","database":"/data/blast/genome"}]
```

Use absolute server-readable paths and provide the FASTA `.fai` index. Gene
lookup accepts full-width characters, surrounding whitespace, case differences
and `gene:` prefixes, and searches GFF3 ID, Name, Alias, gene_id and locus_tag.
Ambiguous aliases report candidate IDs instead of selecting a gene. Version
and transcript suffixes are preserved because removing them can change the target.

The numeric gene-ID helper uses a bounded sample of the selected reference's
annotation: fixed text and leading zeros are supplied for formats such as
Arabidopsis `AT[chromosome]G[number]` and pea gene IDs. Mixed formats retain
manual entry. Profiles can explicitly set `gene_id_format`, for example
`{"prefix":"Psat.cameor.v2.","separator":"g","digits":5,"chromosomes":["1","2","3","4","5","6","7"]}`.
The same local `GET /api/references` information is available for a future
agent layer. Do not commit private inputs, local catalogs or raw private reports.

### Binding maps and saved projects / 結合位置とデータ保存

Design, sequencing and breeding-assay results show forward/reverse binding
positions on the input sequence and the amplified interval. Displayed coordinates
are 1-based and inclusive; genomic maps retain chromosome and strand, including
descending coordinates on minus-strand templates. Expand the sequence view to
inspect the annealing oligos and highlighted differences. dCAPS maps use the
**modified** primers and their own screened coordinates, separately from the
parent pair. CAPS/dCAPS views include recognition sequences, top-strand cut
boundaries, fragment lengths and illustrative gels. Off-target tables show
product coordinates, bp and size differences from the target. A size difference
alone does not establish that bands can be resolved experimentally.

The sequence viewer displays 60-base rows, pages of 1200 bases and a direct
start-index control; bases within the amplicon, F/R footprints, engineered
differences and annotated exons are distinct. Supply a GFF3 to genomic assay or
sequencing jobs to show overlapping genes and each transcript's exon/CDS shape,
with the amplified part and overlapping exon numbers. Gene-based design retains
the same genomic context. Missing annotation and a mismatched chromosome name
are kept separate from a loaded annotation with no overlapping genes.

Enzymes with complete cleavage sites can be expanded to see the **actual AA/BB
product sequence**, shaded recognition bases, and cuts between bases on aligned
5′→3′ and 3′→5′ strands. This includes cuts outside the recognition motif for
Type IIS enzymes. dCAPS uses primer-incorporated products, not the unchanged
reference. Candidate-only gel illustrations explicitly exclude off-target
background and preserve unresolved specificity.

「どこを増幅する？」でF/Rの結合位置・増幅範囲を確認し、配列表示を開くと
設計オリゴと参照配列の違いが見えます。CAPS/dCAPSは酵素・認識配列・切断位置・
断片長とゲル模式図を表示します。オフターゲットは座標・bp・標的とのサイズ差を
一覧表示します。未完了の探索とWet未検証の区別は維持します。

Save inputs and results together with **Save in this browser** (five recent
projects), or export/import a versioned project JSON. Restoring displays stored
results without running a calculation. Browser saves belong to this browser
profile and origin: clearing browser data or changing the server port can make
them unavailable. JSON is the portable backup; referenced genome/DB files remain
external and must be available to recalculate. Saving is explicit, not automatic.

入力条件と結果は「このブラウザーに保存」で直近5件を保持できます。
長期保管・他端末への移動にはプロジェクトJSONを出力してください。
再読込みでは計算しません。ゲノム・DB本体はJSONに含まれず、再計算には元の
ファイルが必要です。ブラウザーのデータ削除やポート変更でブラウザー保存が
見えなくなる場合があります。

配列は60塩基ずつ表示し、1200塩基単位の移動や開始位置の指定ができます。
GFF3を指定すると遺伝子・転写産物ごとのエクソン/CDSの形と、どのエクソンを
増幅するかが重なって見えます。制限酵素は切断部位があるものを開き、実際の
AA/BB配列・認識塩基・上下鎖の切断境界を確認できます。

The visual sequence and project-retention workflow was informed by
[Gene-research](https://github.com/light-suzuki/Gene-research); the current OSS
analysis engine remains authoritative. No additional frontend dependencies are needed.

## Tests

Unit tests are pure Python and need **no external tools or data**:

```bash
pip install -e ".[dev]"
pytest                                # or run the files directly:
python tests/test_specificity.py      # specificity / amplicon pairing
python tests/test_integration.py      # variants, conservation, risk, CAPS
node tests/web_locus.test.js          # optional: coordinate/strand and display regressions
```

## Benchmark

`benchmarks/run_benchmark.py` extracts a real region from an `.fai`-indexed
genome, designs primers, and screens specificity — reporting timing and the
number of predicted products per pair. Point it at your own data with
`export PBO_DBDIR=/path/to/blastdb`. Real runs on a pea genome (design,
multi-cultivar, in-silico PCR, tiling, full assay, CAPS) are written up in
[`benchmarks/RESULTS.md`](benchmarks/RESULTS.md).

`benchmarks/head_to_head_ps2.py` is the automated concordance benchmark against
[PrimerServer2](https://github.com/billzt/PrimerServer2): it designs a pair in
each of N windows across a genome and compares both tools' predicted amplicons
under matched parameters (see [`benchmarks/RESULTS.md`](benchmarks/RESULTS.md) §9).

```bash
python benchmarks/head_to_head_ps2.py --genome tair10.fa --db tair10.fa \
  --primertool /path/to/primertool --n-loci 40 --out h2h.json
```

`benchmarks/continuous_benchmark.py` is the CI-friendly regression benchmark:
it builds a tiny synthetic FASTA/BLAST database and exercises Primer3 design,
BLAST amplicon pairing, duplicate/off-target classification, optional
thermodynamic gating, and multiplex dimer checks.

```bash
python benchmarks/continuous_benchmark.py --max-seconds 30
```

## Contributing

Contributions welcome — see [CONTRIBUTING.md](CONTRIBUTING.md). Changes are
tracked in [CHANGELOG.md](CHANGELOG.md).

## Citation

If you use this in research, please cite it (see [CITATION.cff](CITATION.cff)).

## License

[MIT](LICENSE) © primerblast-oss contributors
