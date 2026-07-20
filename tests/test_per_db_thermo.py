from primerblast_oss.cli import _parse_db_genome_specs, build_parser
from primerblast_oss.design import PrimerPair
from primerblast_oss.pipeline import (
    resolve_genome_for_database,
    run_pipeline,
    thermo_metadata,
)


class FakeGenome:
    def __init__(self, fasta, sequence):
        self.fasta = fasta
        self.sequence = sequence

    def fetch(self, subject, start, end, strand="+"):
        assert subject == "shared_contig"
        return self.sequence[start - 1:end]


def _pair():
    return PrimerPair(
        index=0,
        template_id="target",
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


def test_explicit_mapping_selects_correct_genome_for_same_contig_name():
    genome_a = FakeGenome("A.fa", "A" * 200)
    genome_b = FakeGenome("B.fa", "C" * 200)
    mapping = {"dbA": genome_a, "dbB": genome_b}
    selected_a, association_a = resolve_genome_for_database(
        "dbA", ["dbA", "dbB"], genomes_by_db=mapping)
    selected_b, association_b = resolve_genome_for_database(
        "dbB", ["dbA", "dbB"], genomes_by_db=mapping)
    assert selected_a is genome_a
    assert selected_b is genome_b
    assert selected_a.fetch("shared_contig", 1, 3) == "AAA"
    assert selected_b.fetch("shared_contig", 1, 3) == "CCC"
    assert association_a == association_b == "explicit_db_mapping"


def test_legacy_single_genome_is_never_reused_for_secondary_database():
    genome = FakeGenome("design.fa", "A" * 200)
    first, first_association = resolve_genome_for_database(
        "dbA", ["dbA", "dbB"], genome=genome)
    second, second_association = resolve_genome_for_database(
        "dbB", ["dbA", "dbB"], genome=genome)
    assert first is genome
    assert first_association == "legacy_design_database"
    assert second is None
    assert second_association == "unassociated"


def test_db_genome_cli_requires_exact_database_key():
    assert _parse_db_genome_specs(
        ["dbA=A.fa", "dbB=B.fa"], ["dbA", "dbB"]
    ) == {"dbA": "A.fa", "dbB": "B.fa"}
    try:
        _parse_db_genome_specs(["A=A.fa"], ["dbA"])
    except ValueError as error:
        assert "does not exactly match" in str(error)
    else:
        raise AssertionError("non-exact DB key should fail")


def test_parser_accepts_repeated_db_genome():
    parser = build_parser()
    arguments = parser.parse_args([
        "check",
        "--forward", "ACGTACGTACGTACGTACGT",
        "--db", "dbA",
        "--db", "dbB",
        "--db-genome", "dbA=A.fa",
        "--db-genome", "dbB=B.fa",
    ])
    assert arguments.db == ["dbA", "dbB"]
    assert arguments.db_genome == ["dbA=A.fa", "dbB=B.fa"]


def test_pipeline_passes_each_database_its_associated_genome():
    import primerblast_oss.pipeline as pipeline
    import primerblast_oss.dimers as dimers

    genome_a = FakeGenome("A.fa", "A" * 200)
    genome_b = FakeGenome("B.fa", "C" * 200)
    seen = {}

    original_design = pipeline.design_primers
    original_specificity = pipeline.pair_specificity
    original_available = dimers.available
    try:
        pipeline.design_primers = lambda *args, **kwargs: ([_pair()], "")

        def fake_specificity(forward, reverse, database, **kwargs):
            seen[database] = kwargs.get("genome")
            return {
                "db": database,
                "n_products": 1,
                "n_on_target": 1,
                "n_off_target": 0,
                "n_comigrating": 0,
                "gel_distinguishable": True,
                "on_target": [],
                "off_target": [],
                "specific_observed": True,
                "specific": True,
                "specificity_status": "specific",
                "search_completeness": "complete",
            }

        pipeline.pair_specificity = fake_specificity
        dimers.available = lambda: False
        result = run_pipeline(
            "target",
            "A" * 100,
            ["dbA", "dbB"],
            genomes_by_db={"dbA": genome_a, "dbB": genome_b},
            thermo_params=object(),
        )
    finally:
        pipeline.design_primers = original_design
        pipeline.pair_specificity = original_specificity
        dimers.available = original_available

    assert seen == {"dbA": genome_a, "dbB": genome_b}
    per_db = result.pairs[0].specificity["per_db"]
    assert per_db[0]["thermo_genome_fasta"] == "A.fa"
    assert per_db[1]["thermo_genome_fasta"] == "B.fa"


def test_thermo_metadata_marks_unassociated_database_as_skipped():
    import primerblast_oss.thermo as thermo

    original_available = thermo.available
    try:
        thermo.available = lambda: True
        metadata = thermo_metadata(
            None, object(), True, "unassociated")
    finally:
        thermo.available = original_available
    assert metadata["thermo_status"] == "skipped_no_associated_genome"
    assert metadata["thermo_evaluated"] is False
    assert metadata["thermo_genome_fasta"] is None
