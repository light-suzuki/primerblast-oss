from primerblast_oss.assay import reclassify_by_anchor
from primerblast_oss.design import PrimerPair
from primerblast_oss.pipeline import PipelineResult, _score_pair
from primerblast_oss.report import to_text, to_tsv
from primerblast_oss.risk import assess_risk
from primerblast_oss.specificity import (
    Amplicon,
    SEARCH_COMPLETE,
    SEARCH_POSSIBLY_TRUNCATED,
    SEARCH_REPEAT_LIMITED,
    SpecParams,
    _classify_hit_list,
)


def _pair():
    return PrimerPair(
        index=0,
        template_id="t",
        forward="ACGTACGTACGTACGTACGT",
        reverse="TGCATGCATGCATGCATGCA",
        left_start=0,
        left_len=20,
        right_start=99,
        right_len=20,
        product_size=100,
        tm_f=60.0,
        tm_r=60.0,
        gc_f=50.0,
        gc_r=50.0,
    )


def _db_result(completeness=SEARCH_COMPLETE, specific=True, db="db"):
    amp = Amplicon("chr1", 1, 100, 100, "F", "R", 0, 0, on_target=True)
    return {
        "db": db,
        "n_products": 1,
        "n_on_target": 1,
        "n_off_target": 0,
        "n_comigrating": 0,
        "gel_distinguishable": True,
        "on_target": [amp],
        "off_target": [],
        "specific_observed": True,
        "specific": specific,
        "specificity_status": (
            "specific" if specific is True
            else "indeterminate" if specific is None
            else "non_specific"
        ),
        "search_completeness": completeness,
        "primer_search_completeness": {
            "F": completeness,
            "R": SEARCH_COMPLETE,
        },
        "raw_hits_per_primer": {"F": 95, "R": 1},
        "unique_subjects_per_primer": {"F": 95, "R": 1},
        "completeness_recommendation": "rerun exhaustively",
    }


def test_hit_list_completeness_thresholds():
    params = SpecParams(
        max_target_seqs=100,
        high_copy_hit_threshold=1000,
        high_copy_site_threshold=500,
    )
    assert _classify_hit_list(10, 2, 10, params)[2] == SEARCH_COMPLETE
    assert _classify_hit_list(95, 2, 95, params)[2] == SEARCH_POSSIBLY_TRUNCATED
    assert _classify_hit_list(100, 2, 100, params)[2] == SEARCH_REPEAT_LIMITED
    assert _classify_hit_list(10, 500, 10, params)[2] == SEARCH_REPEAT_LIMITED


def test_complete_clean_search_can_receive_rank_a():
    pair = _pair()
    _score_pair(pair, [_db_result()])
    assert pair.specificity["specificity_status"] == "specific"
    assert pair.specificity["specific_all_db"] is True
    assert pair.specificity["rank"] == "A"


def test_one_incomplete_database_makes_summary_indeterminate():
    pair = _pair()
    complete = _db_result(db="complete_db")
    incomplete = _db_result(
        completeness=SEARCH_POSSIBLY_TRUNCATED,
        specific=None,
        db="limited_db",
    )
    _score_pair(pair, [complete, incomplete])
    assert pair.specificity["specificity_status"] == "indeterminate"
    assert pair.specificity["specific_all_db"] is False
    assert pair.specificity["specific_observed_all_db"] is True
    assert pair.specificity["search_complete_all_db"] is False
    assert pair.specificity["rank"] == "I"
    assert pair.specificity["incomplete_databases"] == ["limited_db"]


def test_observed_offtarget_remains_definitively_non_specific():
    pair = _pair()
    result = _db_result(
        completeness=SEARCH_REPEAT_LIMITED,
        specific=False,
    )
    result["specific_observed"] = False
    result["n_off_target"] = 1
    result["n_comigrating"] = 1
    _score_pair(pair, [result])
    assert pair.specificity["specificity_status"] == "non_specific"
    assert pair.specificity["rank"] != "A"
    assert pair.specificity["rank"] != "I"


def test_anchor_does_not_turn_incomplete_search_into_specific_true():
    amp = Amplicon("chr1", 100, 468, 369, "F", "R", 0, 0)
    result = reclassify_by_anchor(
        {
            "db": "db",
            "on_target": [amp],
            "off_target": [],
            "search_completeness": SEARCH_POSSIBLY_TRUNCATED,
        },
        "chr1",
        50,
        600,
        369,
    )
    assert result["specific_observed"] is True
    assert result["specific"] is None
    assert result["specificity_status"] == "indeterminate"


def test_incomplete_search_can_never_be_low_risk():
    risk = assess_risk(search_completeness=SEARCH_POSSIBLY_TRUNCATED)
    assert risk.level != "low"
    assert any("not exhaustive" in reason for reason in risk.reasons)


def test_text_and_tsv_expose_completeness():
    pair = _pair()
    incomplete = _db_result(
        completeness=SEARCH_POSSIBLY_TRUNCATED,
        specific=None,
    )
    _score_pair(pair, [incomplete])
    result = PipelineResult(
        template_id="t",
        template_len=100,
        pairs=[pair],
        primer3_explain="",
        databases=["db"],
    )
    text = to_text(result)
    tsv = to_tsv(result)
    assert "INDETERMINATE" in text
    assert "rerun exhaustively" in text
    assert "search_completeness" in tsv.splitlines()[0]
    assert SEARCH_POSSIBLY_TRUNCATED in tsv
