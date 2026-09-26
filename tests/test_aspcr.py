"""Tests for allele-specific PCR / ARMS design logic."""
from primerblast_oss.aspcr import (
    build_aspcr,
    generate_as_primers,
    screen_aspcr,
    tetra_arms_from_pair,
)
from primerblast_oss.assay import _template_allele_base
from primerblast_oss.regions import GenomicRegion, Template
from primerblast_oss.specificity import Amplicon
from primerblast_oss.design import PrimerPair


def _sequence():
    # Balanced synthetic target with enough flank on both sides.
    return ("ACGT" * 80)


def _pair():
    return PrimerPair(
        index=0,
        template_id="target",
        forward="ACGTACGTACGTACGTACGT",
        reverse="TGCATGCATGCATGCATGCA",
        left_start=50,
        left_len=20,
        right_start=250,
        right_len=20,
        product_size=201,
        tm_f=60.0,
        tm_r=60.0,
        gc_f=50.0,
        gc_r=50.0,
    )


def test_as_primers_end_at_snp_for_forward_and_reverse_roles():
    seq = _sequence()
    snp = 150
    ref = seq[snp]
    alt = next(base for base in "ACGT" if base != ref)
    candidates = generate_as_primers(
        seq,
        snp,
        ref,
        alt,
        min_length=20,
        max_length=20,
        deliberate_positions=(2,),
    )
    assert candidates
    for candidate in candidates:
        assert candidate["primer_3p"] == snp
        assert candidate["intended_profile"]["terminal_match"] is True
        assert candidate["non_target_profile"]["terminal_mismatch"] is True
        assert candidate["deliberate_mismatch"]["position_from_3prime"] == 2
        assert candidate["intended_profile"]["last3"] == 1
        assert candidate["non_target_profile"]["last3"] >= 2


def test_multiple_snp_substitution_classes_generate_both_alleles_and_roles():
    seq = list(_sequence())
    snp = 150
    for ref, alt in (("A", "G"), ("C", "T"), ("A", "C"), ("G", "T")):
        seq[snp] = ref
        candidates = generate_as_primers(
            "".join(seq),
            snp,
            ref,
            alt,
            min_length=20,
            max_length=20,
            deliberate_positions=(3,),
        )
        assert {(c["allele"], c["role"]) for c in candidates} == {
            ("ref", "F"), ("ref", "R"), ("alt", "F"), ("alt", "R")
        }


def test_classical_aspcr_build_has_ref_and_alt_reactions():
    seq = list(_sequence())
    snp = 150
    seq[snp] = "A"
    result = build_aspcr(
        "".join(seq),
        _pair(),
        snp,
        "G",
        max_classical_per_allele=2,
        max_tetra=2,
        min_length=20,
        max_length=20,
        opt_length=20,
    )
    assert result["status"] == "candidates_found"
    assert result["best_classical_ref"]["allele"] == "ref"
    assert result["best_classical_alt"]["allele"] == "alt"
    assert result["best_classical_ref"]["product_size"] > 0
    assert result["best_classical_alt"]["product_size"] > 0


def test_tetra_arms_reports_control_and_three_genotype_patterns():
    seq = list(_sequence())
    snp = 150
    seq[snp] = "A"
    tetra = tetra_arms_from_pair(
        "".join(seq),
        _pair(),
        snp,
        "A",
        "G",
        candidates_per_orientation=2,
        min_length=20,
        max_length=20,
        opt_length=20,
    )
    assert tetra
    best = tetra[0]
    assert best["control_product_size"] == _pair().product_size
    assert best["control_product_size"] in best["genotype_bands"]["AA_ref"]
    assert best["control_product_size"] in best["genotype_bands"]["AB"]
    assert best["control_product_size"] in best["genotype_bands"]["BB_alt"]
    assert best["ref_product_size"] in best["genotype_bands"]["AA_ref"]
    assert best["alt_product_size"] in best["genotype_bands"]["BB_alt"]
    assert isinstance(best["gel_scorable"], bool)
    assert best["ladder"] in ("100bp", "1kb")


def test_no_deliberate_mismatch_candidate_is_available_for_classical_screening():
    seq = list(_sequence())
    snp = 150
    seq[snp] = "C"
    candidates = generate_as_primers(
        "".join(seq),
        snp,
        "C",
        "T",
        min_length=20,
        max_length=20,
        deliberate_positions=(None,),
    )
    assert candidates
    assert all(c["deliberate_mismatch"] is None for c in candidates)
    assert all(c["intended_profile"]["total"] == 0 for c in candidates)
    assert all(c["non_target_profile"]["terminal_mismatch"] for c in candidates)


def test_minus_strand_genomic_alt_is_complemented_for_oriented_template():
    region = GenomicRegion("chr1", 100, 200, "-", "gene")
    template = Template(
        id="gene",
        seq="A" * 121,
        region=region,
        ext_start=90,
        ext_end=210,
        anchor_coord=210,
        anchor_strand="-",
        flank=10,
    )
    assert _template_allele_base(template, "A") == "T"
    assert _template_allele_base(template, "C") == "G"


def test_screen_aspcr_screens_best_sets_and_reports_unexpected_products():
    import primerblast_oss.specificity as specificity

    aspcr = {
        "best_classical_ref": {
            "primer": "A" * 20,
            "common_role": "R",
            "common_primer": "T" * 20,
            "product_size": 100,
        },
        "best_classical_alt": {
            "primer": "C" * 20,
            "common_role": "R",
            "common_primer": "T" * 20,
            "product_size": 100,
        },
        "best_tetra": {
            "outer_forward": "A" * 20,
            "outer_reverse": "T" * 20,
            "ref_inner": {"role": "F", "primer": "G" * 20},
            "alt_inner": {"role": "R", "primer": "C" * 20},
            "control_product_size": 200,
            "ref_product_size": 120,
            "alt_product_size": 150,
        },
    }

    original = specificity.in_silico_pcr
    calls = []

    def fake_insilico(primers, database, **kwargs):
        calls.append((database, tuple(sorted(primers))))
        products = [
            Amplicon("chr1", 1, 100, 100, "x", "y", 0, 0),
            Amplicon("chr2", 1, 333, 333, "x", "y", 0, 0),
        ]
        return {
            "db": database,
            "products": products,
            "n_products": len(products),
            "search_completeness": "complete",
            "search_complete": True,
            "thermo_site_stats": {
                "attempted_per_primer": {},
                "evaluated_per_primer": {},
                "unresolved_per_primer": {},
                "gated_per_primer": {},
            },
        }

    specificity.in_silico_pcr = fake_insilico
    try:
        screened = screen_aspcr(aspcr, ["dbA", "dbB"])
    finally:
        specificity.in_silico_pcr = original

    assert len(calls) == 6  # ref, alt, tetra across two databases
    sets = screened["specificity_screen"]["sets"]
    assert sets["classical_ref"]["search_complete_all_db"] is True
    assert sets["classical_ref"]["max_unexpected_products"] == 1
    assert sets["tetra_arms"]["max_unexpected_products"] >= 1
