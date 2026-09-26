"""Pure-Python tests for sequencing / primer-walking orchestration."""
from primerblast_oss.cli import build_parser
from primerblast_oss.design import PrimerPair
from primerblast_oss.regions import GenomicRegion, Template
from primerblast_oss.sequencing import (
    M13_FORWARD_TAIL,
    M13_REVERSE_TAIL,
    build_sequence_plan,
    coverage_summary,
    local_region_for_template,
)


def _pair(index, left, right, product=600, rank="A"):
    pair = PrimerPair(
        index=index,
        template_id="target",
        forward="A" * 20,
        reverse="T" * 20,
        left_start=left,
        left_len=20,
        right_start=right,
        right_len=20,
        product_size=product,
        tm_f=60.0,
        tm_r=60.0,
        gc_f=50.0,
        gc_r=50.0,
    )
    pair.specificity = {
        "rank": rank,
        "score": 95.0,
        "specificity_status": "specific",
        "search_completeness": "complete",
        "total_off_target": 0,
    }
    return pair


def test_coverage_summary_reports_internal_and_terminal_gaps():
    tiles = [
        {"covers": (10, 99)},
        {"covers": (80, 149)},
        {"covers": (170, 219)},
    ]
    result = coverage_summary(tiles, (0, 249))
    assert result["covered_intervals"] == [[10, 149], [170, 219]]
    assert result["gaps"] == [[0, 9], [150, 169], [220, 249]]
    assert result["covered_bases"] == 190
    assert result["target_bases"] == 250
    assert result["full_coverage"] is False


def test_local_region_for_plus_and_minus_templates():
    plus_region = GenomicRegion("chr1", 100, 199, "+", "plus")
    plus = Template(
        id="plus", seq="A" * 140, region=plus_region,
        ext_start=80, ext_end=219, anchor_coord=80, anchor_strand="+",
        flank=20,
    )
    assert local_region_for_template(plus) == (20, 119)

    minus_region = GenomicRegion("chr1", 100, 199, "-", "minus")
    minus = Template(
        id="minus", seq="A" * 140, region=minus_region,
        ext_start=80, ext_end=219, anchor_coord=219, anchor_strand="-",
        flank=20,
    )
    assert local_region_for_template(minus) == (20, 119)


def test_sequence_plan_keeps_annealing_primers_and_adds_m13_only_to_order():
    tiles = [
        {
            "index": 1,
            "covers": (0, 599),
            "gap_to_prev": None,
            "pair": _pair(0, 0, 599),
        },
        {
            "index": 2,
            "covers": (500, 1099),
            "gap_to_prev": 100,
            "pair": _pair(1, 500, 1099),
        },
    ]
    plan = build_sequence_plan(
        tiles,
        "target",
        (0, 1099),
        ["db"],
        requested_overlap=100,
        amplicon_range=(500, 700),
        m13_tails=True,
    )
    first = plan["amplicons"][0]
    second = plan["amplicons"][1]
    assert first["forward"] == "A" * 20
    assert first["order_forward"] == M13_FORWARD_TAIL + ("A" * 20)
    assert first["order_reverse"] == M13_REVERSE_TAIL + ("T" * 20)
    assert first["overlap_to_next"] == 100
    assert second["overlap_to_prev"] == 100
    assert plan["coverage"]["full_coverage"] is True


def test_sequence_plan_maps_minus_strand_genomic_coordinates():
    region = GenomicRegion("chr5", 1000, 1999, "-", "geneX", source="gff3:cds")
    template = Template(
        id="geneX", seq="A" * 1200, region=region,
        ext_start=900, ext_end=2099, anchor_coord=2099, anchor_strand="-",
        flank=100,
    )
    tiles = [{
        "index": 1,
        "covers": (100, 699),
        "gap_to_prev": None,
        "pair": _pair(0, 100, 699),
    }]
    plan = build_sequence_plan(
        tiles, "geneX", (100, 1099), ["db"], genomic_template=template)
    genomic = plan["amplicons"][0]["genomic"]
    assert genomic == {
        "chrom": "chr5", "start": 1400, "end": 1999, "strand": "-"
    }


def test_sequence_cli_accepts_length_overlap_and_gene_target():
    parser = build_parser()
    args = parser.parse_args([
        "sequence",
        "--gene", "PsGene1",
        "--genome", "genome.fa",
        "--gff3", "genome.gff3",
        "--db", "db",
        "--amplicon-size", "600-800",
        "--overlap", "120",
        "--m13-tails",
    ])
    assert args.cmd == "sequence"
    assert args.gene == "PsGene1"
    assert args.amplicon_size == "600-800"
    assert args.overlap == 120
    assert args.m13_tails is True
