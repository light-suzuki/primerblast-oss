# Wet-lab validation panel

Computational concordance benchmarks ([benchmarks/RESULTS.md](../benchmarks/RESULTS.md))
establish that the tool implements its own model consistently. They do **not**
measure the outcome that matters for `rank A` / `specific` / `risk low`:
whether a primer pair produces the predicted bands in a real PCR. This panel
is the bridge between in-silico predictions and experimental outcomes.

## Layout

| File | Purpose |
|---|---|
| `manifest.example.json` | Versioned, machine-readable panel: immutable inputs (`sequence_sha256`), expected observations, protocol metadata, observed outcomes, exclusion provenance |
| `RESULTS.md` | Concise benchmark report linked from the README (currently empty of wet-lab data) |
| `benchmarks/validation_analysis.py` | Reconciles observed vs predicted and computes the panel metrics |

The manifest is **schema 1.x** (key `schema_version`). The analysis script
rejects other schemas.

## Panel composition (target)

Cover, per the design in issue #15:

- unique loci, recent paralogs, tandem duplications, repeats, organellar
  off-targets;
- clean pairs and deliberately marginal pairs spanning the A–D / risk range;
- primer-site SNPs/indels, especially in the terminal 5 bases;
- CAPS markers with observed digest patterns;
- multiple plant genome sizes and at least several cultivars/accessions
  (private cultivar genomes may use hashed identifiers, keeping a fully
  public subset reproducible);
- negative cases: no amplification, weak amplification, multiple bands,
  unexpected band sizes.

## Protocol and scoring rules

- **Prediction blocks are immutable.** Fill `prediction` once from a design
  run at panel-definition time; never edit them afterwards. `sequence_sha256`
  pins the primer inputs — the analysis script flags any mismatch as an
  integrity issue.
- **Detection limit:** absence of a visible gel band is *not* proof that a
  sequence-level off-target does not exist. Record the gel percentage, ladder,
  and the smallest resolved band; flag off-targets predicted near the
  detection limit rather than counting them as experimentally absent.
- **Replicates:** an assay needs ≥ 2 independent replicates to be assessed.
  Record `replicate_agreement`; disagreement blocks the outcome until
  resolved.
- **Statuses:**
  - `confirmed` — intended band at the expected size, no extra bands;
  - `weak_amplification` — intended band faint but reproducible;
  - `failed_amplification` — no product;
  - `extra_bands` — additional bands beyond the intended one;
  - `unexpected_sizes` — band(s) present but sizes do not match prediction;
  - `pending` — not yet assessed;
  - `excluded` — never used in metrics; the `exclusion_reason` is
    mandatory and retained for provenance (e.g. failed QC, contaminated run).
- **Blinding:** record `observed` without looking at `prediction` where
  practical; the manifest allows filling observed first and prediction later,
  as long as prediction is finalised before the analysis run that includes
  the assay.
- **Band identity:** when an extra band could be the predicted off-target or
  something else, sequence the amplicon and record `amplicon_sequenced:
  true`.

## Metrics

`benchmarks/validation_analysis.py` reports:

- off-target **sensitivity** (observed extra bands that were predicted);
- off-target **false-positive rate** of predicted off-targets;
- **amplification failure rate** among predicted intended products;
- **calibration** of `rank` and `risk_level` against success probability;
- success stratified by repetitive vs non-repetitive locus class;
- CAPS predicted-vs-observed digest **concordance**.

Denominators are printed explicitly; a cell with no data reads `n/a`.

## How to use

```bash
# after recording observed outcomes in a manifest:
python benchmarks/validation_analysis.py \
  --manifest validation_panel/manifest.example.json \
  --csv /tmp/panel-summary.csv
```

New assays are added to the manifest as JSON entries (copy one from the
example, replace the inputs, and compute the primer hash with:

```bash
python -c "import hashlib,sys; print(hashlib.sha256(('%s|%s' % (sys.argv[1], sys.argv[2])).encode()).hexdigest())" FWD REV
```

When an assay's observed outcome contradicts a prediction, first confirm the
protocol metadata (annealing temperature, Mg2+, cycles) matches the intended
conditions, then keep the entry and record the discrepancy — contradictory
data calibrates the model.

## Excluded / inconclusive provenance

Never delete an assay from the manifest. Mark it `excluded` with an explicit
reason; the analysis lists every exclusion. This keeps the panel auditable.
