# Contributing to primerblast-oss

Thanks for your interest! This is a small, dependency-light tool and
contributions are welcome.

## Development setup

```bash
git clone https://github.com/light-suzuki/primerblast-oss
cd primerblast-oss
pip install -e ".[dev]"
```

Runtime needs the external executables `primer3_core` and BLAST+
(`blastn`, `makeblastdb`) on `PATH`:

```bash
# Debian/Ubuntu
sudo apt install primer3 ncbi-blast+
```

## Running tests

The unit tests are **pure Python and need no external tools or data** — they
cover the specificity pairing, variant/conservation/risk logic, and CAPS:

```bash
pytest                      # or:
python tests/test_specificity.py
python tests/test_integration.py
```

The `benchmarks/` scripts exercise the full pipeline against a real BLAST
database; they require a local genome + BLAST DB and are not run in CI.

## Guidelines

- Standard library only in the core package (the tool shells out to
  `primer3_core` / `blastn`); avoid adding hard Python dependencies.
- Keep coordinate conventions explicit in docstrings (GFF3/VCF are 1-based
  inclusive; BED is 0-based half-open; BLAST/primer positions are documented per
  function).
- Add a test for new logic that can be tested without external tools.
- Keep `CHANGELOG.md` up to date.

## Reporting issues

Please include the command you ran, the `primer3_core` / `blastn` versions
(`primerblast-oss` prints them in the provenance manifest with `--format json`),
and a minimal example. For wrong-result reports, the intended vs predicted
product coordinates are the most useful detail.
