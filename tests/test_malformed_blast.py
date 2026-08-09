"""Regression tests: malformed BLAST outfmt rows must not be silently skipped.

An unparseable row (wrong field count or non-numeric coordinate) previously
incremented the raw hit count but was dropped from the result, biasing the
search toward "no off-targets". These tests pin the new behaviour: malformed
rows are counted, given a short reason, and downgrade search completeness so a
definitive ``specific=True`` can never be produced.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from primerblast_oss.specificity import (  # noqa: E402
    SEARCH_COMPLETE,
    SEARCH_POSSIBLY_TRUNCATED,
    SpecParams,
    _search_metadata,
    _specificity_verdict,
    priming_sites_with_stats,
)


def _blast_output(lines):
    return "\n".join(lines) + "\n"


def _valid_line():
    q = "ACGTACGTACGTACGTACGT"
    return "\t".join([
        "primer", "chr1", "100.0", "20", "0", "0", "1", "20",
        "500", "519", "1e-5", "40", "plus", q, q, "20",
    ])


def _run_with_output(output):
    """Run priming_sites_with_stats against a fixed BLAST output."""
    import primerblast_oss.specificity as specificity

    old = specificity._run_blast
    specificity._run_blast = lambda *_args, **_kwargs: output
    try:
        return priming_sites_with_stats(
            "ACGTACGTACGTACGTACGT", "F", "db",
            SpecParams(max_target_seqs=5000), "blastn")
    finally:
        specificity._run_blast = old


def test_short_row_is_not_silently_skipped():
    short = "\t".join(["primer", "chr1", "100.0", "20", "0", "0", "1", "20",
                       "500", "519", "1e-5", "40", "plus", "ACGT", "ACGT"])
    sites, stats = _run_with_output(_blast_output([short]))
    assert len(sites) == 0
    assert stats.raw_blast_hits == 1
    assert stats.malformed_rows == 1
    assert "15" in stats.malformed_row_reason
    assert stats.completeness != SEARCH_COMPLETE
    assert stats.completeness == SEARCH_POSSIBLY_TRUNCATED


def test_non_numeric_coordinate_is_counted_as_malformed():
    bad = _valid_line().split("\t")
    bad[7] = "not-a-number"  # qend
    sites, stats = _run_with_output(_blast_output(["\t".join(bad)]))
    assert len(sites) == 0
    assert stats.malformed_rows == 1
    assert "non-numeric" in stats.malformed_row_reason
    assert stats.completeness == SEARCH_POSSIBLY_TRUNCATED


def test_clean_16_column_output_is_unchanged():
    lines = [_valid_line(), _valid_line()]
    sites, stats = _run_with_output(_blast_output(lines))
    assert len(sites) == 2
    assert stats.raw_blast_hits == 2
    assert stats.malformed_rows == 0
    assert stats.malformed_row_reason is None
    assert stats.completeness == SEARCH_COMPLETE
    assert stats.priming_sites == 2


def test_mixed_clean_and_malformed_rows_are_not_complete():
    lines = [_valid_line(), "\t".join(["primer", "chr1"])]
    sites, stats = _run_with_output(_blast_output(lines))
    assert len(sites) == 1
    assert stats.raw_blast_hits == 2
    assert stats.malformed_rows == 1
    assert stats.completeness == SEARCH_POSSIBLY_TRUNCATED


def test_malformed_rows_cannot_confirm_specific_true():
    metadata = _search_metadata(
        {
            "F": _run_with_output(_blast_output(["\t".join(["primer", "chr1"])]))[1],
        },
        SpecParams(max_target_seqs=5000),
    )
    assert metadata["search_complete"] is False
    assert metadata["search_completeness"] == SEARCH_POSSIBLY_TRUNCATED
    assert metadata["malformed_rows_per_primer"] == {"F": 1}
    assert metadata["malformed_row_reason_per_primer"]["F"]
    assert _specificity_verdict(True, metadata["search_completeness"]) is None
    assert "unparseable" in metadata["completeness_recommendation"]


def test_row_snippet_is_truncated():
    from primerblast_oss.specificity import _row_snippet

    long_row = "x" * 500
    snippet = _row_snippet(long_row)
    assert len(snippet) == 123  # 120 + "..."
    assert snippet.endswith("...")


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception:
            failed += 1
            print(f"FAIL {fn.__name__}")
            traceback.print_exc()
    print(f"\n{len(fns)-failed}/{len(fns)} passed")
    raise SystemExit(1 if failed else 0)
