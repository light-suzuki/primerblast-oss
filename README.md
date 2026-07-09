# primerblast-oss

[![Release](https://img.shields.io/github/v/release/light-suzuki/primerblast-oss?sort=semver)](https://github.com/light-suzuki/primerblast-oss/releases)
[![CI](https://github.com/light-suzuki/primerblast-oss/actions/workflows/ci.yml/badge.svg)](https://github.com/light-suzuki/primerblast-oss/actions/workflows/ci.yml)
[![Benchmark](https://github.com/light-suzuki/primerblast-oss/actions/workflows/benchmark.yml/badge.svg)](https://github.com/light-suzuki/primerblast-oss/actions/workflows/benchmark.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)

**English** | [日本語](README.ja.md)

A local, open-source, Primer-BLAST-style **command-line tool** for plant breeding
and genetics. It designs PCR primers with **Primer3** and verifies their
**specificity** entirely offline against local BLAST+ databases — including
unpublished genomes and several cultivars at once — and adds in-silico PCR,
whole-region tiling, SNP-under-primer detection, amplicon conservation analysis,
CAPS/dCAPS marker design, and an experimenter-facing risk score.

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
- **Continuous regression benchmark.** CI builds a synthetic FASTA/BLAST database
  and exercises Primer3 design, BLAST amplicon pairing, duplicate/off-target
  classification, thermodynamic gating, and multiplex dimer checks on every push.

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
**multiplex**, **multiplex-design**, **tile**, **assay**, **markers**,
**makedb**. Run `python -m primerblast_oss <subcommand> --help` for the full
option list of any one.

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
wraps every subcommand — no cloud, no third-party Python dependencies (it is built
on the standard-library `http.server`). It is an extra, not the main tool. Run it
on the machine where `primer3_core`, `blastn`, and your BLAST databases live (e.g.
inside WSL):

```bash
python -m primerblast_oss.webapp             # serves http://127.0.0.1:8799, opens a browser
python -m primerblast_oss.webapp --port 9000 --no-browser
```

It exposes tabs for design / in-silico PCR / tiling / assay / QTL markers / build
DB, auto-discovers BLAST databases under `~/.codex/blast_databases`,
`~/blast_databases` and `./databases`, runs each job in a background thread, and
offers one-click TSV/CSV/BED/JSON downloads. English / 日本語 toggle in the header.
It binds to loopback (`127.0.0.1`) only; on WSL2 Windows browsers reach it via the
default localhost forwarding. Everything it does is also available from the CLI.

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
