"""Tests for sequence-aware digestion of predicted off-target PCR products."""
from primerblast_oss.amplicon_digest import (
    digest_offtarget_products,
    reconstruct_pcr_product,
)


class FakeGenome:
    def __init__(self, sequence):
        self.sequence = sequence

    def fetch(self, subject, start, end, strand="+"):
        assert subject == "chr1"
        seq = self.sequence[start - 1:end]
        if strand == "-":
            from primerblast_oss.genome import revcomp
            return revcomp(seq)
        return seq


def test_reconstruct_product_replaces_genomic_ends_with_actual_primers():
    genome = FakeGenome("A" * 40)
    product = {
        "subject": "chr1",
        "start": 1,
        "end": 30,
        "size": 30,
        "orientation": "F/R",
    }
    sequence, error = reconstruct_pcr_product(
        genome,
        product,
        {"F": "GAATTC", "R": "TTTTTT"},
    )
    assert error is None
    assert sequence is not None
    assert len(sequence) == 30
    # The reference genome is all A; EcoRI exists only because the actual
    # forward primer is incorporated into the PCR product.
    assert sequence.startswith("GAATTC")
    assert sequence.endswith("AAAAAA")


def test_digest_offtarget_uses_primer_introduced_restriction_site():
    genome = FakeGenome("A" * 40)
    per_db = [{
        "db": "dbA",
        "n_off_target": 1,
        "products": [{
            "subject": "chr1",
            "start": 1,
            "end": 30,
            "size": 30,
            "orientation": "F/R",
            "on_target": False,
        }],
    }]
    result = digest_offtarget_products(
        per_db,
        {"dbA": genome},
        {"F": "GAATTC", "R": "TTTTTT"},
        "EcoRI",
    )
    digest = result[0]["offtarget_digest"]
    assert digest["complete"] is True
    assert digest["n_digested_products"] == 1
    assert digest["products"][0]["fragments"] == [29, 1]
    assert digest["background_fragments"] == [29, 1]
    assert digest["model"] == (
        "reconstructed_pcr_product_with_primer_incorporation")


def test_missing_database_fasta_makes_offtarget_digest_indeterminate():
    per_db = [{
        "db": "dbA",
        "n_off_target": 1,
        "products": [{
            "subject": "chr1",
            "start": 1,
            "end": 30,
            "size": 30,
            "orientation": "F/R",
            "on_target": False,
        }],
    }]
    result = digest_offtarget_products(
        per_db,
        {},
        {"F": "GAATTC", "R": "TTTTTT"},
        "EcoRI",
    )
    digest = result[0]["offtarget_digest"]
    assert digest["complete"] is False
    assert digest["background_fragments"] == []
    assert digest["unresolved_products"][0]["reason"] == "no_associated_genome"


def test_clean_database_is_complete_even_without_fasta():
    result = digest_offtarget_products(
        [{"db": "dbA", "n_off_target": 0, "products": []}],
        {},
        {"F": "GAATTC", "R": "TTTTTT"},
        "EcoRI",
    )
    digest = result[0]["offtarget_digest"]
    assert digest["complete"] is True
    assert digest["background_fragments"] == []
