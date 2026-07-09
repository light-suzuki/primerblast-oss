# Benchmark results

Machine: EVO-X1 (WSL2 Ubuntu). Tools: `primer3_core` 2.6.x, BLAST+ 2.12.0+.
Genomes: local pea (*Pisum sativum*) cultivar assemblies, ~3.8–3.9 Gbp each,
built as BLASTDB v5.

## 1. Single-genome self-search (correctness)

Template: `chr1:85,337,500–85,339,200` (1,701 bp) extracted from Cameor v2
(`pisum_v2`). Design: product 300–800 bp. Screen: `pisum_v2` (self).

- Design + specificity for 5 pairs: **58.7 s** (~11.7 s/pair; dominated by
  `blastn-short` against a 3.9 Gbp database).
- Top 3 pairs: **rank A, SPECIFIC** — the intended product is recovered as the
  single predicted amplicon (on-target 1, off-target 0).
- Two deliberately weaker pairs are correctly flagged **NON-SPECIFIC**: they
  keep the intended on-target product but also predict 16 off-target products
  (mostly reverse-primer R/R amplicons at degenerate sites with 4+4 mismatches).
  A per-primer BLAST check would have reported both primers as "found" and
  missed these unintended products entirely.

Coordinate check (independent): the predicted on-target amplicon for pair 1,
`chr1:85,338,374–85,338,741` (368 bp), has a 5′ end equal to the forward primer
and a 3′ end equal to the reverse-complement of the reverse primer — exact
match, confirming the genomic mapping.

## 2. Multi-cultivar screening (local, multiple databases)

Same template, screened simultaneously against three cultivar genomes
(`pisum_v2`, `unpublished_cultivar`, `pisum_zw6`) in a single command. The NCBI web tool does
not cover this case: these are local assemblies and it screens one database at a
time. A machine-readable example report is saved at
`benchmarks/example_report.json`; the summary is below.

Design + specificity for 3 pairs across 3 databases: **106 s** (~35 s/pair,
i.e. ~12 s per pair per database).

Highlighted pair (368 bp, `F=GCACTCTAGAGGTTCAAGGCC`, `R=TGGTACGTGTGGTTCAGTTTCA`):

| Database | Result | Note |
|---|---|---|
| `pisum_v2` (Cameor) | **specific** | single intended product |
| `pisum_zw6` (RefSeq ZW6) | **specific** | perfect ortholog, stored on the opposite strand — correctly recognized (R/F) |
| `unpublished_cultivar` (unpublished assembly) | 1 product, **forward primer over a cultivar-specific SNP** (mm 1+0); coordinates withheld (unpublished genome) | allele-specific behavior |

This is the payoff of local multi-database screening: the same primer pair is
clean in two cultivars but overlaps a cultivar-specific SNP in a third — which
requires screening `unpublished_cultivar`, a local unpublished assembly, directly. The two weaker pairs are again flagged with many off-target products
(R/R and F/R amplicons at degenerate sites) in every genome.

## 3. In-silico PCR (`check` mode)

Pasting the pair 1 primers and screening `pisum_v2` (no design step):

```
# pisum_v2: 1 predicted product(s)   sites/primer: {'F': 2425, 'R': 32}
    size   subject:start-end            primers    mm      Δsize
      368  chr1:85338374-85338741      F/R        0+0     -
```

The forward primer seeds at thousands of loci, but only one pairs with a
convergent reverse site inside the size window — a single 368 bp product at the
same coordinates verified in §1. Orientation is not constrained, so F/F and R/R
products would be listed too if present, each with its size gap to neighbours.

## 4. Whole-region tiling (`tile` mode)

Region `chr1:85,337,500–85,339,200` (1,701 bp, `benchmarks/example_template.fa`),
amplicon 400–600 bp, overlap 60: **4 overlapping amplicons** walk the region,
covering `3..1563` (positions 1-based on the template):

| amplicon | covers | size | rank | overlap w/ prev |
|---|---|---|---|---|
| 1 | 3–405 | 403 | A | — |
| 2 | 377–894 | 518 | B (size-resolvable) | 29 bp |
| 3 | 841–1246 | 406 | A | 54 bp |
| 4 | 1137–1563 | 427 | A | 110 bp |

Ranks come from the per-amplicon specificity check (A = unique product,
B = extra products but resolvable by size). The extreme 3' end is not covered
because a primer needs landing sequence *downstream* of what it amplifies —
extract the target with flanking padding for full CDS coverage. This whole-
region walk complements NCBI Primer-BLAST, which designs one amplicon around a
single target rather than tiling a region.

## 5. Full breeding assay (`assay`) — gene + multi-reference + VCF + risk

Target `Psat.cameor.v2.1g00050` CDS (resolved from `pisum_v2.gff3`,
chr1:5,039–6,385, +), template ±100 bp flank, screened across **3 cultivar
genomes** (pisum_v2 / unpublished_cultivar / pisum_zw6) with a synthetic VCF.

**Finding: this gene is multi-copy.** Every one of the 6 designed pairs is
flagged **risk HIGH** — each amplifies its intended locus *plus* 4–8 additional
**perfect-match, same-size products** at paralogous loci. For the 582 bp pair:

```
INTENDED  chr1:5,721-6,302     582 bp  F/R
off       chr1:26,590-27,171   582 bp  F/R  0+0 mismatch
off       chr1:68,561-69,142   582 bp  F/R  0+0 mismatch
off       chr2:366,172,812-... 581 bp  F/R  0+0 mismatch
off       scaffold1040:...     582 bp  F/R  0+0 mismatch
```

The pipeline (a) resolves the gene from GFF3, (b) uses the template's genomic
anchor to mark **only** the on-locus product as intended — the paralogs are
perfect-match products a size/identity heuristic would mislabel as "the intended
one", (c) confirms the amplicon is conserved in all 3 cultivars, and (d) detects
that one pair's reverse primer sits over the synthetic SNP at chr1:6,000. The
verdict — *no specific primer exists for this gene; it has ≥4 genome-wide
copies* — comes from screening the local assembly directly, and illustrates
pain point #9 (paralogs/duplications).

## 6. CAPS/dCAPS + output formats

A SNP-anchored assay (`--snp chr1:5878 --alt C`) forces the amplicon to span the
SNP (Primer3 target), builds the two allele amplicons, and scans 40 restriction
enzymes for a differential digest. The SNP creates an **MseI** site: the
recommended pair's product cuts to **[293, 137] bp for one allele vs [430] bp
(uncut) for the other** — a clean, gel-scorable CAPS marker. (All pairs here are
flagged risk HIGH because this locus has paralogs — the tool warns rather than
hiding it.)

Every output format renders from a single run:
**text / JSON / CSV / BED / HTML report / oligo order table / ASCII off-target
map**. See `benchmarks/assay_caps.json`, `.csv`, `.bed`, and `.html`.

A note on gel resolution: an earlier candidate SNP (chr1:5608) only shifted a
HinfI band by 23 bp — below the 25 bp default — and was correctly reported as
*not* distinguishable, illustrating that the CAPS scan judges real gel
resolvability rather than just "a site changed".

## 7. Head-to-head vs PrimerServer2

Genome: **Lotus japonicus MG20** (`lotja.MG20.gnm3`, ~455 Mbp, published), fresh
BLAST db. Three primer pairs designed by primerblast-oss on `Lj1:20,000,000+`
were run through **both** tools' specificity check against the same genome.
[PrimerServer2](https://github.com/billzt/PrimerServer2) 2.0.0b19 uses the same
core recipe (`blastn -task blastn-short -evalue 30000 -word_size 7`, then a
primer3 Tm check; off-target if a site's Tm is within `--Tm-diff` 20 °C of the
primer — matching our `--min-anneal-tm 40` floor).

| pair | product | primerblast-oss (+thermo) | PrimerServer2 |
|---|---|---|---|
| 1 | 573 bp | 1 amplicon @ `Lj1:20,000,521-20,001,093` | 1 amplicon @ same |
| 2 | 569 bp | 1 amplicon @ `Lj1:20,000,521-20,001,089` | 1 amplicon @ same |
| 3 | 570 bp | 1 amplicon @ `Lj1:20,000,521-20,001,090` | 1 amplicon @ same |

**Exact concordance** on amplicon count, coordinates and size for all three
pairs. Pair 3 is the informative case: primerblast-oss's *mismatch-only* model
flags 4 candidate off-targets, but **both tools' thermodynamic check drop them**
as non-amplifying — i.e. the optional thermo gate brings primerblast-oss into
agreement with PrimerServer2. Runtime for the check was a few seconds per tool.

primerblast-oss adds, on top of this parity: multi-database screening, CAPS/dCAPS
design, GFF3/VCF/QTL integration, whole-region tiling, experimenter risk scoring,
and — like PrimerServer2 — multiplex primer-dimer checking (`multiplex`) plus a
compatible-set *designer* (`multiplex-design`). PrimerServer2 still has a mature
hosted web server, which this does not.

*Reproducibility note:* PrimerServer2 needs `samtools`; where it was unavailable
we supplied a tiny `samtools faidx` shim backed by the same `.fai` reader. The
BLAST/primer3 computations are PrimerServer2's own.

## 8. Three-way head-to-head — NCBI Primer-BLAST **and** PrimerServer2 (Arabidopsis)

Genome: **Arabidopsis thaliana** (published TAIR10). Template `At1_5Mb`
([`benchmarks/at1_5mb.fa`](at1_5mb.fa)) = a 1,501 bp window at `chr1:5,000,000`.
The **same template / primers** were run through three tools on 2026-07-09:

- the live [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
  web service (database *Genomes for selected eukaryotic organisms*, organism
  limited to *Arabidopsis thaliana* — taxid 3702; job `TUeQPRoBF6kwk4eWivajpPDtspbd_qmL3A`),
- [PrimerServer2](https://github.com/billzt/PrimerServer2) 2.0.0b19
  `primertool check` against a local TAIR10 db, and
- primerblast-oss `design` / `check` against the same local TAIR10 db.

**Top de-novo pair — byte-for-byte identical across all three tools:**

| | forward | reverse | product | Tm (F/R) | specificity |
|---|---|---|---|---|---|
| NCBI Primer-BLAST | `GACAAGGAATCAGCGGCTCT` | `GCAGCGTTTTGTAGTGGGTG` | 342 bp | 60.11 / 60.04 | specific (no other targets) |
| PrimerServer2 | `GACAAGGAATCAGCGGCTCT` | `GCAGCGTTTTGTAGTGGGTG` | 342 bp | 60.11 / 60.04 | 1 amplicon `1:5000421-5000762` |
| primerblast-oss | `GACAAGGAATCAGCGGCTCT` | `GCAGCGTTTTGTAGTGGGTG` | 342 bp | 60.1 / 60.0 | rank A, 1 product `1:5000421-5000762` |

**Specificity concordance across three primer pairs** (primerblast-oss `check` vs
PrimerServer2 `check`, same TAIR10 db, both with their thermodynamic model):

| pair | forward / reverse | PrimerServer2 | primerblast-oss |
|---|---|---|---|
| atA | `GACAAGGAATCAGCGGCTCT` / `GCAGCGTTTTGTAGTGGGTG` | 1× 342 bp `1:5000421-5000762` | 1× 342 bp `1:5000421-5000762` |
| atB | `GGACGAAGCAGGAGATGGAG` / `GCAGCGTTTTGTAGTGGGTG` | 1× 317 bp `1:5000446-5000762` | 1× 317 bp `1:5000446-5000762` |
| atC | `GGACGAAGCAGGAGATGGAG` / `TGATCCTCCTTACACGCAGC` | 2× (332 bp `…777`, 356 bp `…801`) | 2× (332 bp `…777`, 356 bp `…801`) |

**Exact concordance** on amplicon count, size, and genomic coordinates for every
pair — including the two-product `atC` case, where both tools independently find
the same pair of co-migrating-ish amplicons (332 / 356 bp). PrimerServer2's Tm
values match NCBI's to the second decimal (both use primer3 `oligotm`); the
predicted amplicon coordinates match primerblast-oss exactly.

- **Design parity:** primerblast-oss's #1 de-novo pair is the pair NCBI returns as
  its Primer pair 1 (template positions 422..441 / 744..763; product 342 bp).
- **Beyond both:** for the same pair primerblast-oss additionally reports the
  forward×reverse primer-dimer / hairpin ΔG (F×R ΔG −2.89 kcal/mol, worst −4.82 —
  OK), and can pick a mutually compatible *multiplex* set across several targets
  (`multiplex-design`) — PrimerServer2 offers multiplex-dimer *checking*, and NCBI
  neither.

On this published Arabidopsis locus primerblast-oss **matches both NCBI
Primer-BLAST and PrimerServer2** on the primer it picks and on the specificity
verdict, while adding dimer thermodynamics and multiplex-set selection on top.

### 8b. Six-locus NCBI Primer-BLAST panel

To go beyond one locus, six pairs drawn from the automated benchmark below
(§9) — three "clean" single-product loci, one duplicated locus, and the two loci
where primerblast-oss and PrimerServer2 disagreed — were each submitted to the
live NCBI Primer-BLAST service (organism *Arabidopsis thaliana*, taxid 3702) and
compared with both local tools on the same primers (2026-07-09):

| locus | product | NCBI Primer-BLAST | PrimerServer2 | primerblast-oss |
|---|---|---|---|---|
| 1_2767379 | 448 bp | specific (1) | 1 | 1 |
| 2_7387233 | 172 bp | specific (1) | 1 | 1 |
| 3_2607814 | 383 bp | specific (1) | 1 | 1 |
| 1_5533258 | 433 bp | specific (1) | 2 | 2 |
| 1_13830895 | 303 bp | specific (1) | 2 (+365) | **1** |
| 4_18583553 | 350 bp | not specific (+3530 Mt) | 2 (+1980) | 1 |

On the three clean loci **all three tools agree exactly**. The other three are the
informative cases, and each disagreement is a *single* off-target sitting on one
tool's acceptance threshold — no tool is a strict superset of another:

- **1_13830895** — PrimerServer2 reports an extra 365 bp off-target; **NCBI and
  primerblast-oss both reject it** (the site's reverse primer is not 3'-anchored).
  Here NCBI confirms primerblast-oss's stricter 3'-end rule against PrimerServer2.
- **1_5533258** — primerblast-oss and PrimerServer2 both report a second 433 bp
  product ~7 kb from the target whose reverse primer anneals weakly (Tm ≈ 45 °C,
  2 mismatches); NCBI drops it. Here the two local tools agree and NCBI is the
  stricter one on a low-Tm site.
- **4_18583553** — NCBI flags a 3,530 bp forward/forward product on the
  *mitochondrion* (primer 3' ends mismatched) that primerblast-oss's 3'-anchor
  rule rejects; PrimerServer2 instead flags a different 1,980 bp nuclear site that
  NCBI does not. All three differ by one borderline product.

So on real Arabidopsis loci primerblast-oss sits **squarely within the NCBI /
PrimerServer2 range** — matching NCBI against PrimerServer2 on one locus and
PrimerServer2 against NCBI on another, with no systematic over- or under-calling.
Remaining differences are borderline off-targets near a Tm or 3'-alignment
threshold, plus the fact that NCBI screens its own *Arabidopsis* representative
assembly rather than the local TAIR10 FASTA — not algorithmic errors.

## 9. Automated multi-locus concordance vs PrimerServer2 (40 Arabidopsis loci)

The head-to-heads above use a handful of hand-checked pairs. To measure
concordance at scale without cherry-picking, `benchmarks/head_to_head_ps2.py`
selects **40 low-N 1,500 bp windows spread across all five TAIR10 chromosomes**,
designs a Primer3 pair in each, and predicts that pair's amplicons with **both**
tools against the same TAIR10 BLAST database. Parameters were matched on both
sides: off-target size window **50–2,000 bp**, PrimerServer2 `--Tm-diff 20`
against primerblast-oss's thermodynamic gate (primer3-py), and PrimerServer2's
amplicon cap raised so neither tool truncates. Concordance is scored on the exact
set of predicted amplicons (count, product size, and genomic coordinates).

| locus set | n | exact-coordinate concordance |
|---|---|---|
| all designed loci | 40 | **33 / 40 (82 %)** |
| non-repetitive loci (≤ 3 predicted amplicons) | 36 | **33 / 36 (92 %)** |
| repetitive / multi-copy loci (> 3 amplicons) | 4 | 0 / 4 (see below) |

The residual disagreements are **understood, not random**, and fall into two
groups:

- **Repetitive / multi-copy loci (4/40).** Both tools agree the primer is
  non-specific (many products); they differ only in *how many* individual repeat
  copies or tandem-repeat "ladder" rungs each enumerates. Any workflow would
  reject these primers regardless, so the exact count is moot.
- **Marginal 3'-end sites (3/36 non-repetitive loci).** In every case
  PrimerServer2 reports **one extra off-target** that primerblast-oss does not,
  and in every case that off-target depends on a primer binding site whose **3'
  end is not fully aligned** — e.g. locus `4_18583553`, where the reverse primer
  aligns only 18 of 20 bases (3'-terminal 2 nt unaligned), or `1_13830895`, where
  the forward primer's 3' end is unaligned (both confirmed by `blastn-short`).
  primerblast-oss's priming model **requires the 3' end to anchor** (terminal
  match + limited 3'-window mismatch), so it treats such sites as non-amplifying;
  PrimerServer2 keeps them because the overall duplex Tm is high enough. This is a
  deliberate model difference, not a defect: primerblast-oss is the stricter side
  on 3'-priming competence (biochemically the dominant factor for extension),
  while PrimerServer2 is more conservative and flags more candidate off-targets.

So on realistic, non-repetitive primers the two tools agree on **92 %** of loci
exactly, and the ~8 % that differ do so for one explainable reason — the 3'-anchor
vs Tm-window priming threshold — with primerblast-oss taking the stricter,
extension-competent-only interpretation. No discordance traced to an
implementation error.

## 10. Live NCBI direct-primer regression

The external service is deliberately not called from CI. Instead,
[`ncbi_live_regression.json`](ncbi_live_regression.json) records one polite,
browser-submitted smoke comparison with fixed public inputs. On 2026-07-10 JST,
the official NCBI Primer-BLAST job
`fXegCnrAd2hQVnJTfzNWYQUoR1MoO1xOKQ` checked
`GGTTCACCGCTCTCACTCAA` / `CTCCCTCAGTTGCAACCCAT` against
`NC_045512.2` with the RefSeq representative-genomes database restricted to
SARS-CoV-2. NCBI reported one specific 258 bp product at `28428-28685`; the
current local `check` run against the same accession FASTA also returned exactly
one 258 bp product at the same coordinates (0+0 mismatches).

This is an external integration regression only, not evidence of plant-genome
or NCBI-wide equivalence. The JSON records the source FASTA SHA-256, inputs,
job ID, and expected comparison result so a future browser run can detect drift.

## Reproduce

```bash
DB=/home/user/.codex/blast_databases
# design (single- and multi-db)
python benchmarks/run_benchmark.py
python benchmarks/run_benchmark.py --db $DB/pisum_v2 --db $DB/unpublished_cultivar --db $DB/pisum_zw6
# in-silico PCR
python -m primerblast_oss check \
  --forward GCACTCTAGAGGTTCAAGGCC --reverse TGGTACGTGTGGTTCAGTTTCA --db $DB/pisum_v2
# whole-region tiling
python -m primerblast_oss tile --template-fasta benchmarks/example_template.fa \
  --amplicon-min 400 --amplicon-max 600 --overlap 60 --candidates-per-tile 4 --db $DB/pisum_v2
# full breeding assay (gene + 3 refs + VCF)
python -m primerblast_oss assay --gene Psat.cameor.v2.1g00050 --gene-feature cds \
  --gff3 $DB/pisum_v2/pisum_v2.gff3 --genome $DB/pisum_v2.fa \
  --db $DB/pisum_v2 --db $DB/unpublished_cultivar --db $DB/pisum_zw6 \
  --vcf benchmarks/synthetic_variants.vcf --flank 100 --product-size 150-600
# CAPS marker across a SNP
python -m primerblast_oss assay --snp chr1:6385 --alt A --genome $DB/pisum_v2.fa \
  --db $DB/pisum_v2 --flank 250 --format text
# NCBI head-to-head (Arabidopsis): de-novo design + specificity vs a local TAIR10 db
python -m primerblast_oss design --template-fasta benchmarks/at1_5mb.fa \
  --template-id At1_5Mb --product-size 300-400 --num-return 5 \
  --db $DB/tair10 --genome-fasta $DB/tair10.fa
# multiplex-set design across several targets (beyond NCBI Primer-BLAST)
python -m primerblast_oss multiplex-design --template-fasta benchmarks/multi_targets.fa \
  --product-size 80-300 --candidates-per-target 5 --require-specific \
  --db $DB/tair10 --genome-fasta $DB/tair10.fa
# PrimerServer2 side of the Arabidopsis head-to-head (section 8), same TAIR10 db:
#   printf "atA GACAAGGAATCAGCGGCTCT GCAGCGTTTTGTAGTGGGTG\n" > q.txt
#   primertool check q.txt $DB/tair10.fa -t q.tsv        # PrimerServer2 2.0.0b19
# Automated 40-locus concordance vs PrimerServer2 (section 9):
python benchmarks/head_to_head_ps2.py \
  --genome $DB/tair10.fa --db $DB/tair10.fa \
  --primertool /path/to/primertool \
  --n-loci 40 --check-size-min 50 --check-size-max 2000 --out h2h.json
```
