from primerblast_oss.caps import (
    CapsResult,
    ENZYME_METADATA,
    caps_scan,
    cut_events,
    dcaps_candidates,
    digest_fragment_sizes,
    materialize_dcaps_primers,
)
from primerblast_oss.design import PrimerPair
from primerblast_oss.outputs import order_table


def test_ecori_uses_real_top_and_bottom_cut_positions():
    sequence = "GGGAATTCCC"
    events = cut_events(sequence, ENZYME_METADATA["EcoRI"])
    assert len(events) == 1
    assert events[0].site_pos == 2
    assert events[0].top_cut == 3
    assert events[0].bottom_cut == 7
    fragments = digest_fragment_sizes(sequence, ENZYME_METADATA["EcoRI"])
    assert fragments == [7, 3]
    assert sum(fragments) == len(sequence)


def test_bsai_type_iis_cuts_outside_recognition_sequence():
    sequence = "AAAAGGTCTCAAAAA"
    events = cut_events(sequence, ENZYME_METADATA["BsaI"])
    assert len(events) == 1
    assert events[0].site_pos == 4
    assert events[0].top_cut == 11       # 4 + len(GGTCTC) + 1
    assert events[0].bottom_cut == 15    # 4 + len(GGTCTC) + 5
    assert events[0].complete is True


def test_dcaps_search_finds_nonpalindromic_reverse_orientation():
    candidates = dcaps_candidates(
        sequence="GAGACT",
        snp_index=5,
        ref_base="C",
        alt_base="T",
        enzymes={"BsaI": ENZYME_METADATA["BsaI"]},
        max_primer_mismatch=0,
    )
    assert any(
        candidate["orientation"] == "-"
        and candidate["recognition_oriented"] == "GAGACC"
        and candidate["present_in"] == "ref"
        for candidate in candidates
    )


def test_dpni_is_not_recommended_for_unmethylated_pcr_product():
    results = caps_scan("AAAAGATCAAAA", "AAAAAATCAAAA")
    assert all(result.enzyme != "DpnI" for result in results)
    assert ENZYME_METADATA["DpnI"].pcr_compatible is False
    assert ENZYME_METADATA["MboI"].pcr_compatible is True


def test_materialized_dcaps_primer_leaves_snp_outside_oligo():
    sequence = "A" * 50
    candidate = {
        "enzyme": "EcoRI",
        "recognition": "GAATTC",
        "orientation": "+",
        "present_in": "alt",
        "mismatches": 1,
        "engineered_changes": [{"position": 18, "from": "A", "to": "G"}],
    }
    primers = materialize_dcaps_primers(
        sequence, snp_index=20, candidates=[candidate])
    assert len(primers) == 1
    primer = primers[0]
    assert primer["primer_role"] == "F"
    assert primer["primer_end"] == 19
    assert primer["three_prime_distance_to_snp"] == 1
    assert primer["primer_sequence"][18 - primer["primer_start"]] == "G"
    assert len(primer["primer_sequence"]) == 20


def test_end_to_end_dcaps_rechecks_modified_pair_against_each_database():
    import primerblast_oss.dcaps_workflow as workflow
    import primerblast_oss.assay as assay
    import primerblast_oss.dimers as dimers
    import primerblast_oss.amplicon_digest as amplicon_digest

    class Template:
        seq = "A" * 120

    original = PrimerPair(
        index=0,
        template_id="target",
        forward="A" * 20,
        reverse="T" * 20,
        left_start=0,
        left_len=20,
        right_start=99,
        right_len=20,
        product_size=100,
        tm_f=58.0,
        tm_r=58.0,
        gc_f=0.0,
        gc_r=0.0,
    )
    frame = {
        "enzyme": "EcoRI",
        "recognition": "GAATTC",
        "orientation": "+",
        "present_in": "alt",
        "mismatches": 1,
        "engineered_changes": [{"position": 49, "from": "A", "to": "G"}],
    }
    materialized = {
        **frame,
        "primer_role": "F",
        "primer_sequence": "A" * 19 + "G",
        "primer_start": 30,
        "primer_end": 49,
        "primer_length": 20,
        "orderable": True,
    }
    digest = CapsResult(
        enzyme="EcoRI",
        recognition="GAATTC",
        allele_a_fragments=[50, 20],
        allele_b_fragments=[70],
        distinguishable=True,
        min_gel_gap=50,
    )
    seen = []

    original_frames = workflow.dcaps_candidates
    original_materialize = workflow.materialize_dcaps_primers
    original_caps_scan = workflow.caps_scan
    original_specificity = workflow.pair_specificity
    original_analyze = assay.analyze_pair
    original_available = dimers.available
    original_digest = amplicon_digest.digest_offtarget_products
    try:
        workflow.dcaps_candidates = lambda *args, **kwargs: [frame]
        workflow.materialize_dcaps_primers = lambda *args, **kwargs: [materialized]
        workflow.caps_scan = lambda *args, **kwargs: [digest]

        def fake_specificity(forward, reverse, database, **kwargs):
            seen.append((database, forward, reverse, kwargs["designed_size"]))
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

        workflow.pair_specificity = fake_specificity
        assay.analyze_pair = lambda *args, **kwargs: {
            "specific": False,
            "specificity_status": "non_specific",
            "specific_all_db": False,
            "specificity_status_all_db": "non_specific",
            "intended_status": "unique",
            "risk": "high",
            "search_complete_all_db": True,
            "variant_in_primer_3prime": False,
            "dimers": None,
            "per_db_products": [
                {
                    "db": "dbA",
                    "n_off_target": 1,
                    "products": [
                        {"size": 70, "on_target": True},
                        {"size": 200, "on_target": False},
                    ],
                },
                {
                    "db": "dbB",
                    "n_off_target": 1,
                    "products": [
                        {"size": 70, "on_target": True},
                        {"size": 220, "on_target": False},
                    ],
                },
            ],
        }
        dimers.available = lambda: False

        def fake_digest(per_db_products, _genomes, _primers, enzyme):
            assert enzyme == "EcoRI"
            output = []
            for view in per_db_products:
                copied = dict(view)
                background = [200] if view["db"] == "dbA" else [220]
                copied["offtarget_digest"] = {
                    "enzyme": "EcoRI",
                    "complete": True,
                    "background_fragments": background,
                    "products": [],
                    "unresolved_products": [],
                }
                output.append(copied)
            return output

        amplicon_digest.digest_offtarget_products = fake_digest
        result = workflow.evaluate_dcaps_candidates(
            Template(),
            original,
            snp_local=50,
            alt_base="G",
            databases=["dbA", "dbB"],
            max_candidates_to_screen=2,
        )
    finally:
        workflow.dcaps_candidates = original_frames
        workflow.materialize_dcaps_primers = original_materialize
        workflow.caps_scan = original_caps_scan
        workflow.pair_specificity = original_specificity
        assay.analyze_pair = original_analyze
        dimers.available = original_available
        amplicon_digest.digest_offtarget_products = original_digest

    assert len(seen) == 2
    assert {entry[0] for entry in seen} == {"dbA", "dbB"}
    assert all(entry[1] == materialized["primer_sequence"] for entry in seen)
    assert result["n_orderable"] == 1
    assert result["best"]["specificity"]["specific_all_db"] is False
    assert result["best"]["gel_scorable_all_db"] is True
    assert result["best"]["recommendation_status"] == (
        "orderable_with_gel_separated_offtargets")


def test_order_sheet_uses_validated_engineered_primer():
    pair = {
        "name": "P1",
        "forward": "ORIGINALF",
        "reverse": "ORIGINALR",
        "caps_enzyme": "EcoRI",
        "marker_type": "dCAPS",
        "caps": {
            "best_marker_type": "dCAPS",
            "dcaps": {
                "best": {
                    "orderable": True,
                    "enzyme": "EcoRI",
                    "forward": "ENGINEEREDF",
                    "reverse": "ORIGINALR",
                    "tm_f": 60.0,
                    "tm_r": 59.0,
                    "gc_f": 50.0,
                    "gc_r": 50.0,
                    "product_size": 150,
                    "modified_primer_role": "F",
                    "engineered_mismatches": 1,
                }
            },
        },
    }
    output = order_table([pair])
    assert "ENGINEEREDF" in output
    assert "ORIGINALF" not in output
    assert "1 mismatch(es)" in output


def test_order_sheet_uses_tetra_arms_four_primer_set():
    pair = {
        "name": "SNP1",
        "forward": "OUTERF",
        "reverse": "OUTERR",
        "tm_f": 60.0,
        "tm_r": 60.0,
        "preferred_genotyping_mode": "tetra-ARMS",
        "aspcr": {
            "best_tetra": {
                "outer_forward": "OUTERF",
                "outer_reverse": "OUTERR",
                "ref_inner": {
                    "role": "F",
                    "primer": "REFINNER",
                    "tm": 59.0,
                    "deliberate_mismatch": {
                        "position_from_3prime": 2,
                        "from": "A",
                        "to": "C",
                    },
                },
                "alt_inner": {
                    "role": "R",
                    "primer": "ALTINNER",
                    "tm": 58.5,
                    "deliberate_mismatch": {
                        "position_from_3prime": 3,
                        "from": "G",
                        "to": "T",
                    },
                },
            }
        },
    }
    output = order_table([pair])
    for sequence in ("OUTERF", "OUTERR", "REFINNER", "ALTINNER"):
        assert sequence in output
    assert "tetra-ARMS" in output
    assert "deliberate -2" in output
    assert "deliberate -3" in output


def test_order_sheet_uses_classical_ref_and_alt_as_primers():
    pair = {
        "name": "SNP2",
        "forward": "PARENTF",
        "reverse": "PARENTR",
        "tm_f": 60.0,
        "tm_r": 60.0,
        "preferred_genotyping_mode": "AS-PCR",
        "aspcr": {
            "best_classical_ref": {
                "role": "F",
                "primer": "REFAS",
                "tm": 58.0,
                "common_role": "R",
                "common_primer": "COMMONR",
                "deliberate_mismatch": None,
            },
            "best_classical_alt": {
                "role": "F",
                "primer": "ALTAS",
                "tm": 58.5,
                "common_role": "R",
                "common_primer": "COMMONR",
                "deliberate_mismatch": {
                    "position_from_3prime": 2,
                    "from": "A",
                    "to": "G",
                },
            },
        },
    }
    output = order_table([pair])
    assert "REFAS" in output and "ALTAS" in output
    assert "COMMONR" in output
    assert "PARENTF" not in output
    assert "AS-PCR" in output
