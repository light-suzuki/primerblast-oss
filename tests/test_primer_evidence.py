"""Literal oligos, label swaps and changed-sequence hypotheses stay distinct."""
from primerblast_oss.primer_evidence import annotate_input_evidence


PRIMERS = {"F": "ACGA", "R": "TTGC"}


def amp(left="F", right="R", start=1):
    return {"subject": "chr1", "start": start, "end": start + 99,
            "fwd_primer": left, "rev_primer": right, "fwd_end3": start + 3,
            "rev_end3": start + 96, "fwd_mismatch": 0, "rev_mismatch": 0}


def row(products, changed=(), complete=True, db="ref"):
    return {"db": db, "oligos": {"F": "TCGT" if "F" in changed else "ACGA",
                                  "R": "GCAA" if "R" in changed else "TTGC"},
            "reverse_complemented_inputs": list(changed), "products": products,
            "search_complete": complete}


def test_swapped_labels_predict_amplification_without_sequence_change():
    product = amp("R", "F")
    assessments = annotate_input_evidence([row([product])], PRIMERS)
    assert assessments[0]["distinct_primer_products"] == 1
    assert assessments[0]["input_error"] == "undetermined"
    assert product["input_evidence"]["sequence_status"] == "as_supplied"
    assert product["input_evidence"]["label_orientation"] == "swapped"
    left, right = product["primer_bindings"]
    assert (left["primer"], left["reference_strand"], left["extends"], left["end3"]) == ("R", "+", "right", 4)
    assert (right["primer"], right["reference_strand"], right["extends"], right["end3"]) == ("F", "-", "left", 97)


def test_reverse_complement_candidate_is_not_literal_amplification():
    candidate = amp()
    summaries = annotate_input_evidence([row([]), row([candidate], ["R"])], PRIMERS)
    assert summaries[0]["as_supplied_products"] == 0
    assert summaries[0]["reverse_complement_candidates"] == 1
    assert candidate["input_evidence"]["changed_primers"] == ["R"]
    assert candidate["input_evidence"]["as_supplied_locus_observed"] is False
    assert candidate["primer_bindings"][1]["input_sequence"] == "TTGC"
    assert candidate["primer_bindings"][1]["oligo"] == "GCAA"


def test_unused_changed_oligo_does_not_change_self_product_prediction():
    self_product = amp("F", "F")
    summaries = annotate_input_evidence([row([amp("F", "F")]), row([self_product], ["R"])], PRIMERS)
    assert summaries[0]["same_primer_products"] == 1
    assert summaries[0]["distinct_primer_products"] == 0
    assert summaries[0]["reverse_complement_candidates"] == 0
    assert self_product["input_evidence"]["sequence_status"] == "as_supplied"
    assert self_product["input_evidence"]["label_orientation"] == "same_primer"
    assert len(self_product["primer_bindings"]) == 2


def test_incomplete_original_search_remains_unresolved_per_database():
    candidate = amp()
    summaries = annotate_input_evidence([row([], complete=False), row([amp()], db="other"),
                                        row([candidate], ["R"])], PRIMERS)
    assert summaries[0]["original_search_complete"] is False
    assert summaries[1]["original_search_complete"] is True
    assert candidate["input_evidence"]["original_search_complete"] is False
    assert candidate["input_evidence"]["as_supplied_locus_observed"] is False


def test_changed_candidate_does_not_hide_literal_product_at_same_locus():
    candidate = amp()
    annotate_input_evidence([row([amp()]), row([candidate], ["R"])], PRIMERS)
    assert candidate["input_evidence"]["as_supplied_locus_observed"] is True
    assert candidate["input_evidence"]["sequence_status"] == "reverse_complement_candidate"


def test_binding_candidates_include_outward_sites_even_without_products(monkeypatch):
    from primerblast_oss import specificity as S, report
    sites = [S.PrimingSite("F", "chr1", "-", 1, 0, 0, plen=4),
             S.PrimingSite("R", "chr1", "+", 100, 0, 0, plen=4)]
    monkeypatch.setattr(S, "_detect_blastn", lambda _: "blastn")
    monkeypatch.setattr(S, "screen_primers_with_stats", lambda *args: (sites, {}))
    monkeypatch.setattr(S, "annotate_thermo", lambda *args: (sites, {}, {"evaluated_per_primer": {}}))
    raw = S.in_silico_pcr(PRIMERS, "ref")
    result = report.insilico_to_dict([raw], PRIMERS)["results"][0]
    assert result["products"] == []
    assert result["binding_site_counts"] == {"F": {"+": 0, "-": 1}, "R": {"+": 1, "-": 0}}
    assert result["binding_sites"][0]["extends"] == "left"
    assert result["binding_sites"][1]["end5"] == 97
    assert result["binding_sites"][0]["thermo_viable"] is None


def test_binding_candidate_limit_does_not_truncate_counts(monkeypatch):
    from primerblast_oss import specificity as S
    sites = [S.PrimingSite("F", "chr1", "+", i + 20, 0, 0) for i in range(250)]
    monkeypatch.setattr(S, "_detect_blastn", lambda _: "blastn")
    monkeypatch.setattr(S, "screen_primers_with_stats", lambda *args: (sites, {}))
    monkeypatch.setattr(S, "annotate_thermo", lambda *args: (sites, {}, {"evaluated_per_primer": {}}))
    raw = S.in_silico_pcr(PRIMERS, "ref")
    assert len(raw["binding_sites"]) == 200 and raw["binding_sites_truncated"] == 50
    assert raw["binding_site_counts"]["F"]["+"] == 250
