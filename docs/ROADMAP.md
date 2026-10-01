# Next priorities / 次に進める課題

Triage checked against open GitHub issues on 2026-10-02. This is a priority
guide, not a claim that the acceptance criteria have been completed.

| Priority | Issue | Next useful step | Completion boundary |
|---|---|---|---|
| First code task / まず進める実装 | [#46: stream/index large VCF](https://github.com/light-suzuki/primerblast-oss/issues/46) | Stream plain/gzip VCF, query indexed intervals when available, push region selection into assay/marker workflows | Preserve variant selection, keep a dependency-free fallback, measure memory and wall time on bounded synthetic inputs before larger runs |
| Follow-up acceptance / 実装済み機能の受入 | [#34: batch BLAST](https://github.com/light-suzuki/primerblast-oss/issues/34) | Check remaining malformed-row attribution and continuous batch coverage, then isolate real-genome memory/runtime measurements | Batching already exists. The recorded ten-primer real-genome check did not demonstrate speed or memory improvement; 100/500-primer acceptance remains open |
| Experimental work / 実験待ち | [#33: prospective wet-PCR panel](https://github.com/light-suzuki/primerblast-oss/issues/33) | Freeze predictions before measuring PCR bands, failures and protocol metadata | Requires assessed laboratory assays; synthetic or literature-sequence comparisons cannot close it |
| After #33 / 実験結果が集まってから | [#35: rank/risk calibration](https://github.com/light-suzuki/primerblast-oss/issues/35) | Evaluate existing score classes against independent outcomes; use hold-out/fresh assays for revised thresholds | Do not fit thresholds to the tool's own predictions or present heuristic classes as calibrated probabilities |

The nine unadjudicated TAIR10/PrimerServer2 disagreements remain a correctness
follow-up for #34. [Benchmark sections 11–12](../benchmarks/RESULTS.md)
distinguish that comparison, synthetic regressions and the reported literature
panels. A smaller independent review can advance disagreement attribution
without claiming full real-genome or wet acceptance.

GUI usability follow-ups include pair-by-pair CSV/TSV import, clearer error
locations and reference DB/FASTA/GFF3 association checks. CLI follow-ups include
independent `pair_id,F,R` batch reactions and resumable/per-pair results.
Standalone BLAST and Primer3 entrypoints stay in the GUI; integrated PCR
workflows remain available through the CLI.

まず#46を小さく進め、#34は残る低コストな正しさの確認から着手するのが妥当です。
#33と#35は実験が必要なので、計算結果だけで完了扱いにはしません。
