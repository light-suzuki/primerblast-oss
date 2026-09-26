"""Gel-aware CAPS/dCAPS scoring tests."""
from primerblast_oss.caps import ENZYME_METADATA, caps_scan
from primerblast_oss.cli import build_parser
from primerblast_oss.gel import (
    analyze_digest_patterns,
    best_analysis_from_assay,
    virtual_gel_svg,
)


def test_cut_uncut_pattern_prefers_100bp_ladder_and_common_dense_gel():
    analysis = analyze_digest_patterns([430], [293, 137])
    assert analysis["ladder"] == "100bp"
    assert analysis["gel_percent"] == 2.0
    assert analysis["diagnostic_a"] == [430]
    assert analysis["diagnostic_b"] == [293, 137]
    assert [band["size"] for band in analysis["heterozygote_bands"]] == [
        430, 293, 137]
    assert analysis["score"] >= 70


def test_large_digest_uses_1kb_ladder():
    analysis = analyze_digest_patterns([5000], [3000, 2000])
    assert analysis["ladder"] == "1kb"
    assert analysis["gel_percent"] in (1.0, 1.2)


def test_tiny_diagnostic_fragment_is_penalized_and_reported():
    easy = analyze_digest_patterns([500], [300, 200])
    tiny = analyze_digest_patterns([500], [470, 30])
    assert tiny["score"] < easy["score"]
    assert any("below 50 bp" in reason for reason in tiny["reasons"])


def test_custom_ladder_bands_are_used_exactly():
    analysis = analyze_digest_patterns(
        [900], [600, 300],
        ladder="custom",
        custom_ladder_bands=[50, 300, 600, 900, 1200],
        gel_percent=2.0,
    )
    assert analysis["ladder"] == "custom"
    assert analysis["ladder_bands"] == [50, 300, 600, 900, 1200]
    assert analysis["gel_percent"] == 2.0


def test_virtual_gel_has_marker_and_three_genotype_lanes():
    analysis = analyze_digest_patterns([430], [293, 137])
    svg = virtual_gel_svg(analysis)
    assert svg.startswith("<svg")
    assert "100%%" not in svg
    assert 'width="100%"' in svg
    for lane in ("M", "AA", "AB", "BB"):
        assert ">%s</text>" % lane in svg
    assert "gel score" in svg


def test_caps_scan_attaches_gel_analysis():
    ref = "A" * 100 + "GAATTC" + "A" * 100
    alt = "A" * 100 + "GACTTC" + "A" * 100
    results = caps_scan(
        ref,
        alt,
        enzymes={"EcoRI": ENZYME_METADATA["EcoRI"]},
        gel_min_gap=25,
    )
    assert len(results) == 1
    result = results[0]
    assert result.gel_analysis is not None
    assert result.gel_analysis["lanes"]["AA"]
    assert result.gel_analysis["lanes"]["AB"]
    assert result.gel_analysis["lanes"]["BB"]


def test_best_analysis_from_assay_handles_natural_and_dcaps():
    analysis = analyze_digest_patterns([430], [293, 137])
    natural = {
        "pairs": [{
            "caps": {
                "best_marker_type": "CAPS",
                "best_result": {"gel_analysis": analysis},
            }
        }]
    }
    assert best_analysis_from_assay(natural) is analysis

    dcaps = {
        "pairs": [{
            "caps": {
                "best_marker_type": "dCAPS",
                "dcaps": {"best": {"digest": {"gel_analysis": analysis}}},
            }
        }]
    }
    assert best_analysis_from_assay(dcaps) is analysis


def test_assay_cli_accepts_gel_options():
    parser = build_parser()
    args = parser.parse_args([
        "assay",
        "--snp", "chr1:1000",
        "--alt", "G",
        "--genome", "genome.fa",
        "--db", "db",
        "--ladder", "100bp",
        "--gel-percent", "2.0",
        "--virtual-gel", "marker.svg",
    ])
    assert args.ladder == "100bp"
    assert args.gel_percent == "2.0"
    assert args.virtual_gel == "marker.svg"
