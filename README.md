# primerblast-oss

[![CI](https://github.com/light-suzuki/primerblast-oss/actions/workflows/ci.yml/badge.svg)](https://github.com/light-suzuki/primerblast-oss/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)

A local, open-source reimplementation of [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
for plant breeding and genetics. Design PCR primers with **Primer3** and check
their **specificity** entirely offline against local BLAST+ databases —
including unpublished genomes and several cultivars at once — plus in-silico
PCR, whole-region tiling, SNP-under-primer detection, amplicon conservation,
CAPS/dCAPS marker design, and experimenter risk scoring.

> Pure-Python core (standard library only); shells out to `primer3_core` and
> BLAST+. Unit tests need **no external tools or data**.

## Why local and open source

NCBI Primer-BLAST is excellent, but it is **not open source** and runs only as a
hosted web service. In practice that means you can't audit, fork, or self-host
it, you can't run it next to your data, and unpublished or embargoed genomes
can't be submitted to it.

Relying on any external service also ties your pipeline to its availability and
policies (rate limits, maintenance, occasional outages). Keeping the workflow
**local and offline** removes that dependency and makes runs fully
reproducible — pinned FASTA, GFF3, VCF, BLAST DB, and tool versions, with no
queue or login. That matters most for the case the web tool can't serve anyway:
local, unpublished, multi-cultivar genomes.

primerblast-oss is **MIT-licensed**, so anyone can read, run, and build on it.

## Why this exists

NCBI Primer-BLAST does two things: (1) Primer3 designs candidate primers, and
(2) BLAST screens each primer against a database and **pairs the hits into
predicted amplicons** to flag unintended PCR products. Most local "primer +
BLAST" scripts only do a per-primer BLAST and miss step (2) — the part that
actually detects off-target products.

`primerblast-oss` implements step (2) faithfully and adds capabilities NCBI's
web tool cannot offer:

| | NCBI Primer-BLAST | primerblast-oss |
|---|---|---|
| Primer3 design | ✅ | ✅ |
| Pair BLAST hits into predicted amplicons | ✅ | ✅ |
| Off-target products from either primer as F **or** R (F/F, R/R, F/R) | ✅ | ✅ |
| 3'-end anchored priming model (terminal match + 3'-window mismatch limit) | ✅ | ✅ |
| Runs offline on **unpublished / local** genomes | ❌ | ✅ |
| Screen against **multiple databases** in one run | ❌ | ✅ |
| **In-silico PCR** from pasted primers (orientation-free) | ❌ | ✅ |
| **Tile a whole region** with overlapping amplicons | ❌ (clusters on one side) | ✅ |
| Gel-resolvability of off-targets (size-gap aware) | ❌ | ✅ |
| Full control of every threshold, scriptable, no queue | partial | ✅ |
| No login, no rate limit | ❌ | ✅ |

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

## Requirements

- Python ≥ 3.8 (standard library only)
- `primer3_core` (Debian/Ubuntu: `apt install primer3`)
- BLAST+ `blastn` / `makeblastdb` (`apt install ncbi-blast+`)
- A nucleotide BLAST database (see below)

## Install

```bash
pip install -e .        # provides the `primerblast-oss` command
# or run without installing:
python -m primerblast_oss --help
```

## Usage

Four subcommands: **design**, **check**, **tile**, **makedb**.

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

Fixes Primer-BLAST's habit of clustering primers on one side: this walks the
entire region with overlapping amplicons (e.g. to sequence a whole gene).

```bash
python -m primerblast_oss tile \
  --template-fasta gene.fa \
  --amplicon-min 400 --amplicon-max 700 --overlap 60 \
  --db /path/to/genome_db
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
  --db $DB/cameor_v2 --db $DB/JI2694 --db $DB/ZW6 \
  --vcf variants.vcf --flank 100 --product-size 150-600 --format html --out report.html

# a CAPS marker spanning a SNP (alt allele given)
python -m primerblast_oss assay --snp chr1:6385 --alt A \
  --genome genome.fa --db $DB/cameor_v2 --flank 250 --format text
```

### `markers` — evenly spaced markers across a QTL interval

```bash
python -m primerblast_oss markers --interval chr1:80000000-90000000 \
  --genome genome.fa --db $DB/cameor_v2 --n-markers 20 --format json
```

### `makedb` — build a database (with `-parse_seqids`)

```bash
python -m primerblast_oss makedb genome.fa --out-db genome_db
```

## How this addresses NCBI Primer-BLAST's weak points

| # | NCBI Primer-BLAST weakness | primerblast-oss |
|---|---|---|
| 1 | weak at batch / many regions | `markers`, `assay` over BED/gene lists; CLI + library, scriptable |
| 2 | poor fit for local / custom assemblies | everything runs on local FASTA + BLAST DB; `makedb` helper |
| 3 | weak multi-reference comparison | `--db` repeatable; amplicon **conservation** scored per reference |
| 4 | whole-chromosome design is clumsy | `tile` + `markers` generate primers across whole regions/intervals |
| 5 | primer strand/orientation unclear | every binding site reports strand, 5'/3' coords, extension direction |
| 6 | unexpected side products hard to read | BLAST hits **paired into predicted amplicons** with sizes |
| 7 | F-F / R-R products hard to see | enumerated explicitly and shown in the ASCII map & tables |
| 8 | 3'-end mismatch not visible | explicit 3'-terminal **5 bp / 10 bp** mismatch counts per hit |
| 9 | paralogs / repeats / duplications | genome-wide pairing surfaces every duplicated priming site |
| 10 | not built for CAPS/dCAPS | `caps` scan: enzymes that digest two alleles differently, gel gap |
| 11 | weak GFF3 / VCF / QTL integration | `--gene`/`--gff3`, `--vcf`, `--interval`, BED input |
| 12 | opaque empty results | Primer3 explain string surfaced; per-stage diagnostics |
| 13 | weak reproducibility | provenance manifest pins tool versions, params, DB fingerprints |
| 14 | weak experimenter-facing scoring | `risk` rolls up every signal into low/medium/high with reasons |
| 15 | side products not visualized | ASCII off-target map + BED track for a genome browser |

Design cues were taken from
[PrimerServer2](https://github.com/billzt/PrimerServer2) (strand-aware BLAST-hit
pairing, multi-threaded `blastn`, coordinate input) and NCBI Primer-BLAST.

### Key options (shared by design/check/tile)

| Option | Meaning | Default |
|---|---|---|
| `--product-size` (design) | one or more ranges, e.g. `150-500,500-1000` | `70-1000` |
| `--amplicon-min/--amplicon-max/--overlap` (tile) | tiling geometry | 400/800/40 |
| `--opt-tm/--min-tm/--max-tm` | primer melting temperature window | 60/57/63 |
| `--max-total-mismatch` | mismatches allowed for an off-target to still prime | 4 |
| `--max-3prime-mismatch` | mismatches allowed in the 3' window | 1 |
| `--three-prime-window` | size of the 3' window | 5 |
| `--max-product` | largest off-target amplicon considered amplifiable | 4000 |
| `--gel-min-gap` | size gap (bp) to resolve two products on a gel | 50 |
| `--no-3prime-terminal` | allow a mismatched 3' terminal base | off |
| `--format` | `text` \| `json` \| `tsv` (design only) | text |

### Verdicts & scoring

Design/tile rank pairs **A–D**: **A** = single product in every database;
**B** = extra products but all far enough in size to resolve on a gel;
**C/D** = one or more off-targets that **co-migrate** with the intended band.
Co-migrating off-targets are penalized heavily; size-resolvable ones only
lightly — matching how such pairs are used in practice.

## GUI-ready output

`--format json` emits a stable schema intended to back a GUI. Every object
carries what a UI needs to render without recomputation:

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

## Tests

Unit tests are pure Python and need **no external tools or data**:

```bash
pip install -e ".[dev]"
pytest                                # or run the files directly:
python tests/test_specificity.py      # specificity / amplicon pairing
python tests/test_integration.py      # variants, conservation, risk, CAPS
```

## Benchmark

`benchmarks/run_benchmark.py` extracts a real region from an `.fai`-indexed
genome, designs primers, and screens specificity — reporting timing and the
number of predicted products per pair. Point it at your own data with
`export PBO_DBDIR=/path/to/blastdb`. Real runs on a pea genome (design,
multi-cultivar, in-silico PCR, tiling, full assay, CAPS) are written up in
[`benchmarks/RESULTS.md`](benchmarks/RESULTS.md).

## Contributing

Contributions welcome — see [CONTRIBUTING.md](CONTRIBUTING.md). Changes are
tracked in [CHANGELOG.md](CHANGELOG.md).

## Citation

If you use this in research, please cite it (see [CITATION.cff](CITATION.cff)).

## License

[MIT](LICENSE) © primerblast-oss contributors
