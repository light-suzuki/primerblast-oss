# Reproduce the README screenshots

Captured in the Codex in-app browser on 2026-10-02 from a source checkout.
These are actual GUI screenshots, cropped where useful. No output text or
scientific measurements were edited. All sequences, references, identifiers
and coordinates are synthetic, from `benchmarks/continuous_benchmark.py`.
`ExampleGene` is a synthetic GFF3 annotation, not an experimental gene finding.

From a source checkout, with BLAST+, Primer3 and `.[thermo]` installed:

```bash
python scripts/readme_demo.py --port 8880
```

Open `http://localhost:8880/`. This server discovers only its generated
synthetic database and reference profile. It does not scan the normal local
reference catalog. Generated FASTA/BLAST/GFF3 files stay in ignored
`.local/readme-demo/`; no large or private data belongs in these screenshots.

1. Capture the task selection screen in Japanese and English using the language
   button: `start-ja.jpg`, `start-en.jpg`.
2. Open PCR check. Select **Synthetic PCR demo** to fill matching paths and the DB.
3. Paste F=`ACGTTGCAAGTCCGATCGTA`, R=`CGTAAGCTAGCATACGGTCA`.
   Leave automatic reverse-complement alternatives **OFF** for the input shot:
   `input-ja.jpg`, `input-en.jpg`.
4. Turn alternatives **ON** and calculate. The original oligos yield a 120 bp
   F/F product. Changing R to its reverse complement adds two 103 bp F/R products.
   Collapse the sequence-form comparison and capture the two result sections:
   `results-ja.jpg`, `results-en.jpg`.
5. Expand the `chr_target:1–103` candidate for input/changed oligos, actual
   extension directions, ExampleGene and reference FASTA:
   `product-ja.jpg`, `product-en.jpg`.

Language changes rebuild the result details; reopen/collapse the relevant
sections after switching languages. Screenshots illustrate interface behavior,
not wet-PCR success or exhaustive off-target discovery.
