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
```
