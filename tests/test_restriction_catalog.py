from primerblast_oss.caps import (
    ENZYME_METADATA, caps_scan, cut_events, enzyme_records, result_to_dict,
)
from primerblast_oss.restriction_catalog import catalog, enzyme_info


def test_snapshot_has_provenance_and_all_records():
    data = catalog()
    assert len(data["enzymes"]) == 1088
    assert len(data["source_sha256"]) == 64
    assert data["rebase_version"] == "404 (2024)"
    assert len({row["name"] for row in data["enzymes"]}) == 1088


def test_equivalent_enzymes_are_not_confused_with_different_cleavage():
    assert "HpaII" in enzyme_info("MspI")["same_cut_enzymes"]
    assert "Acc65I" in enzyme_info("KpnI")["different_cut_enzymes"]
    assert "Acc65I" not in enzyme_info("KpnI")["same_cut_enzymes"]
    assert enzyme_info("KpnI")["cuts"] != enzyme_info("Acc65I")["cuts"]


def test_product_names_full_width_and_case_are_resolved():
    assert enzyme_info(" ＥｃｏＲＩ－ＨＦ ")["name"] == "EcoRI"
    assert enzyme_info("bsai-hfv2")["name"] == "BsaI"
    assert enzyme_info("EcoRI-made-up") is None
    assert cut_events("GGGAATTCCC", "EcoRI-HF") == cut_events("GGGAATTCCC", ENZYME_METADATA["EcoRI"])
    assert enzyme_records({"ＥｃｏＲＩ－ＨＦ": "GAATTC"}, recommended_only=True)[0].name == "EcoRI"


def test_catalog_geometry_handles_upstream_and_multiple_cuts():
    eco = enzyme_info("EcoRI")["pattern"]
    assert eco["top"].replace(" ", "") == "5′G|AATTC3′"
    assert eco["bottom"].replace(" ", "") == "3′CTTAA|G5′"
    bsai = enzyme_info("BsaI")
    assert bsai["cuts"] == [[7, 11]]
    assert bsai["pattern"]["flank_right"] == 5
    aju = enzyme_info("AjuI")
    assert len(aju["cuts"]) == 2
    assert aju["pattern"]["flank_left"] == 12
    assert aju["pattern"]["top"].count("|") == 2
    assert aju["prediction_supported"] is False


def test_unknown_geometry_and_special_substrates_not_automatically_recommended():
    records = enzyme_records(recommended_only=True)
    assert len(records) > 200
    assert all(record.name not in {"DpnI", "MspJI", "AbaSI", "EcoP15I", "AjuI"} for record in records)
    assert enzyme_info("Aba13301I")["pattern"] is None
    assert not enzyme_info("MspJI")["prediction_supported"]
    geometries = [(record.recognition, record.top_cut, record.bottom_cut) for record in records]
    assert len(geometries) == len(set(geometries))


def test_extended_caps_candidates_include_actual_patterns_and_relationships():
    # BsmFI was absent from the original 40-enzyme panel.
    enzyme = ENZYME_METADATA["BsmFI"]
    sequence = "AAAAGGGAC" + "A" * 40
    events = cut_events(sequence, enzyme)
    assert [(event.top_cut, event.bottom_cut) for event in events] == [(19, 23)]
    results = caps_scan(sequence, sequence.replace("GGGAC", "GGGAT"), enzymes={"BsmFI": enzyme})
    assert len(results) == 1
    data = result_to_dict(results[0])
    assert data["pattern"]
    assert sum(data["allele_a_fragments"]) == len(sequence)
    assert data["reference_url"].startswith("https://identifiers.org/rebase:")


def test_custom_geometry_does_not_inherit_catalog_cut_diagram():
    from primerblast_oss.caps import RestrictionEnzyme
    enzyme = RestrictionEnzyme("EcoRI", "GAATTC", 2, 4)
    result = caps_scan("GGGAATTCCC", "GGGACTTCCC", enzymes={"EcoRI": enzyme})[0]
    assert "pattern" not in result_to_dict(result)
