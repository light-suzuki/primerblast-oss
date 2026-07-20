from types import SimpleNamespace

from primerblast_oss.assay import reclassify_by_anchor
from primerblast_oss.cli import _thermo_setup, build_parser
from primerblast_oss.design import PrimerPair
from primerblast_oss.outputs import pairs_to_csv
from primerblast_oss.regions import GenomicRegion, Template
from primerblast_oss.specificity import (
    Amplicon, PrimerHitStats, PrimingSite, SpecParams,
    annotate_thermo, pair_specificity,
)


def _template():
    return Template(
        id="t", seq="A" * 500,
        region=GenomicRegion("chr1", 1000, 1499, "+", "t"),
        ext_start=1000, ext_end=1499,
        anchor_coord=1000, anchor_strand="+",
    )


def _pair():
    return PrimerPair(
        index=0, template_id="t", forward="A" * 20, reverse="T" * 20,
        left_start=100, left_len=20, right_start=299, right_len=20,
        product_size=200, tm_f=60.0, tm_r=60.0, gc_f=0.0, gc_r=0.0,
    )


def test_anchored_dcaps_allows_expected_intentional_mismatch():
    amp = Amplicon("chr1", 1100, 1299, 200, "F", "R", 1, 0)
    result = reclassify_by_anchor(
        {"db": "db", "on_target": [], "off_target": [amp],
         "search_completeness": "complete"},
        template=_template(), pair=_pair(),
        expected_primer_mismatches={"F": 1},
    )
    assert result["intended_status"] == "unique"
    assert result["specific"] is True
    assert result["on_target"] == [amp]


def test_pair_specificity_accepts_configured_engineered_mismatch():
    import primerblast_oss.specificity as specificity

    sites = [
        PrimingSite("F", "chr1", "+", 119, 1, 0, plen=20),
        PrimingSite("R", "chr1", "-", 280, 0, 0, plen=20),
    ]
    stats = {
        "F": PrimerHitStats("F", 1, 1, 1, False, False),
        "R": PrimerHitStats("R", 1, 1, 1, False, False),
    }
    old_detect = specificity._detect_blastn
    old_screen = specificity.screen_primers_with_stats
    try:
        specificity._detect_blastn = lambda value: "blastn"
        specificity.screen_primers_with_stats = lambda *args: (sites, stats)
        result = pair_specificity(
            "A" * 20, "T" * 20, "db", designed_size=200,
            sp=SpecParams(min_product=40, max_product=1000),
            allowed_primer_mismatches={"F": 1},
        )
    finally:
        specificity._detect_blastn = old_detect
        specificity.screen_primers_with_stats = old_screen
    assert result["specific"] is True
    assert result["n_on_target"] == 1


def test_no_thermo_sentinel_never_calls_evaluator():
    import primerblast_oss.thermo as thermo

    class Genome:
        fasta = "x.fa"
        def fetch(self, *args):
            return "A" * 20

    old_available = thermo.available
    old_evaluate = thermo.evaluate
    try:
        thermo.available = lambda: True
        thermo.evaluate = lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("evaluator must not be called"))
        kept, viable, stats = annotate_thermo(
            [PrimingSite("F", "chr1", "+", 20, 0, 0, plen=20)],
            {"F": "A" * 20}, Genome(), tp=False, gate=True)
    finally:
        thermo.available = old_available
        thermo.evaluate = old_evaluate
    assert len(kept) == 1 and viable == {}
    assert sum(stats["attempted_per_primer"].values()) == 0


def test_unresolved_contig_is_reported_not_silently_evaluated():
    import primerblast_oss.thermo as thermo

    class MissingGenome:
        fasta = "missing.fa"
        def fetch(self, *args):
            raise KeyError("missing contig")

    old_available = thermo.available
    try:
        thermo.available = lambda: True
        kept, viable, stats = annotate_thermo(
            [PrimingSite("F", "chr1", "+", 20, 0, 0, plen=20)],
            {"F": "A" * 20}, MissingGenome(), tp=None, gate=True)
    finally:
        thermo.available = old_available
    assert len(kept) == 1 and viable == {}
    assert stats["unresolved_per_primer"] == {"F": 1}
    assert stats["evaluated_per_primer"] == {}


def test_cli_no_thermo_uses_explicit_disable_sentinel():
    parser = build_parser()
    args = parser.parse_args([
        "check", "--forward", "A" * 20, "--db", "db", "--no-thermo"])
    genomes, params, gate = _thermo_setup(args)
    assert genomes == {}
    assert params is False and gate is True


def test_design_database_fasta_conflict_is_rejected(tmp_path):
    design = tmp_path / "design.fa"
    other = tmp_path / "other.fa"
    design.write_text(">chr1\nAAAA\n")
    other.write_text(">chr1\nCCCC\n")
    fake_genome = SimpleNamespace(fasta=str(design))
    args = SimpleNamespace(
        db=["db"], db_genome=["db=%s" % other], genome_fasta=None,
        no_thermo=True, min_anneal_tm=40.0, max_3p_dg=-5.0,
        no_thermo_gate=False,
    )
    try:
        _thermo_setup(args, design_genome=fake_genome)
    except ValueError as error:
        assert "must match --genome" in str(error)
    else:
        raise AssertionError("conflicting first-DB FASTA should be rejected")


def test_csv_keeps_historical_columns_and_appends_new_fields():
    text = pairs_to_csv([{
        "name": "P1", "forward": "AAAA", "reverse": "TTTT",
        "product_size": 100, "tm_f": 60.0, "tm_r": 60.0,
        "gc_f": 0.0, "gc_r": 0.0, "tp5_mismatch_min": 1,
        "snp_in_primer": True, "caps_enzyme": "EcoRI",
        "dimers": {"cross_dimer_dg": -2.0}, "risk": "medium",
    }])
    header = text.splitlines()[0].split(",")
    assert header[:20] == [
        "name", "forward", "reverse", "product_size", "tm_f", "tm_r",
        "tm_diff", "gc_f", "gc_r", "n_off_target", "n_ff", "n_rr",
        "n_fr_offtarget", "tp5_mismatch_min", "snp_in_primer",
        "conserved_refs", "caps_enzyme", "gel_distinguishable",
        "cross_dimer_dg", "risk",
    ]
    assert "marker_type" in header
    assert "specificity_status_all_db" in header
