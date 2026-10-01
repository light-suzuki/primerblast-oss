"""Shared standalone tools and coordinate-based product exports."""
from types import SimpleNamespace

import pytest

from primerblast_oss import workflows, sequence_tools
from primerblast_oss.genome import Genome, revcomp
from primerblast_oss.specificity import Amplicon


def fixture_genome(tmp_path):
    sequence = "ACGTTGCAAGTCCGATCGTA" + "GATTACA" * 9 + "CGTAAGCTAGCATACGGTCA"
    path = tmp_path / "synthetic.fa"
    path.write_text(">chr1\n" + sequence + "\n", encoding="ascii")
    # Write bytes so this synthetic index is identical on Windows and Linux.
    path.write_bytes((">chr1\n" + sequence + "\n").encode("ascii"))
    (tmp_path / "synthetic.fa.fai").write_text("chr1\t%d\t6\t%d\t%d\n" % (len(sequence), len(sequence), len(sequence) + 1))
    return path, sequence


def mock_pcr(monkeypatch):
    calls = []
    def run(primers, db, **kwargs):
        calls.append(dict(primers))
        amp = Amplicon("chr1", 1, 103, 103, "F", "R", 0, 0)
        return {"db": db, "sites_per_primer": {name: 1 for name in primers},
                "n_products": 1, "products": [amp], "search_completeness": "incomplete", "search_complete": False}
    monkeypatch.setattr(workflows, "in_silico_pcr", run)
    return calls


def test_auto_orientation_keeps_hypotheses_separate_and_original_input(monkeypatch, tmp_path):
    path, sequence = fixture_genome(tmp_path)
    calls = mock_pcr(monkeypatch)
    forward, reverse = sequence[:20], sequence[-20:]
    output = workflows.execute("check", {"db": ["synthetic"], "genome": str(path),
                                        "forward": forward, "reverse": reverse, "input_orientation": "auto"})
    assert len(calls) == 4 and all(set(call) == {"F", "R"} for call in calls)
    assert {call["R"] for call in calls} == {reverse, revcomp(reverse)}
    assert output["primers"] == {"F": forward, "R": reverse}
    assert all(row["search_complete"] is False for row in output["results"])
    assert all(row["products"][0]["sequence"] == sequence for row in output["results"])
    assert output["results"][-1]["reverse_complemented_inputs"] == ["F", "R"]


def test_ordered_oligos_are_not_transformed(monkeypatch, tmp_path):
    path, sequence = fixture_genome(tmp_path)
    calls = mock_pcr(monkeypatch)
    workflows.execute("check", {"db": ["synthetic"], "genome": str(path),
                                "forward": "  " + sequence[:20].lower() + "\n", "reverse": revcomp(sequence[-20:])})
    assert calls == [{"F": sequence[:20], "R": revcomp(sequence[-20:])}]


def test_auto_orientation_rejects_large_pool_before_search(monkeypatch):
    monkeypatch.setattr(workflows, "in_silico_pcr", lambda *a, **k: pytest.fail("must not search"))
    with pytest.raises(ValueError, match="up to two"):
        workflows.execute("check", {"db": ["synthetic"], "primers": ["A=ACGT", "B=CCGT", "C=GCAT"], "input_orientation": "auto"})


def test_duplicate_names_and_invalid_bases_are_rejected():
    for params in ({"forward": "ACGT", "primers": ["F=CCGT"]}, {"forward": "ACG$"}):
        with pytest.raises(ValueError):
            workflows.execute("check", dict(params, db=["synthetic"]))


def test_product_export_and_gene_overlap_do_not_assume_other_assembly(monkeypatch, tmp_path):
    path, sequence = fixture_genome(tmp_path)
    gff = tmp_path / "synthetic.gff3"
    gff.write_text("##gff-version 3\nchr1\ttest\tgene\t10\t80\t.\t-\t.\tID=g1;Name=ExampleGene\n")
    mock_pcr(monkeypatch)
    monkeypatch.setattr(sequence_tools.shutil, "which", lambda _: None)
    output = workflows.execute("check", {"db": ["a", "b"], "genome": str(path), "gff3": str(gff), "forward": sequence[:20]})
    first, second = [row["products"][0] for row in output["results"]]
    assert first["annotations"]["genes"][0]["name"] == "ExampleGene"
    assert first["fasta"].splitlines()[0] == ">chr1:1-103_reference_plus"
    assert "".join(first["fasta"].splitlines()[1:]) == sequence
    assert second["sequence_status"] == "unavailable"
    assert second["annotations"]["status"] == "not_provided"


def test_extraction_length_failure_is_distinct_from_no_sequence(tmp_path):
    path, sequence = fixture_genome(tmp_path)
    product = {"subject": "chr1", "start": 1, "end": 105, "size": 105}
    sequence_tools.product_sequence(product, "synthetic", Genome(str(path)))
    assert product["sequence_status"] == "unavailable" and "sequence" not in product
    assert "length" in product["sequence_error"]


def test_blastdb_extraction_retains_reference_plus_coordinates(monkeypatch):
    monkeypatch.setattr(sequence_tools.shutil, "which", lambda _: "blastdbcmd")
    def run(command, **kwargs):
        assert command[command.index("-range") + 1] == "10-13"
        assert command[command.index("-strand") + 1] == "plus"
        return SimpleNamespace(returncode=0, stdout=b"ACGT\n", stderr=b"")
    monkeypatch.setattr(sequence_tools.subprocess, "run", run)
    product = {"subject": "chr1", "start": 10, "end": 13, "size": 4}
    sequence_tools.product_sequence(product, "synthetic")
    assert product["sequence"] == "ACGT" and product["sequence_source"] == "blast_database"


def test_regular_blast_preserves_reverse_alignment_and_multiple_queries(monkeypatch):
    monkeypatch.setattr(sequence_tools.shutil, "which", lambda _: "blastn")
    def run(command, **kwargs):
        assert command[command.index("-task") + 1] == "blastn"
        assert command[command.index("-strand") + 1] == "both"
        return SimpleNamespace(returncode=0, stdout=b"query2\tchr1\t100\t4\t0\t0\t1\t4\t13\t10\t1e-6\t20\tACGT\tACGT\n", stderr=b"")
    monkeypatch.setattr(sequence_tools.subprocess, "run", run)
    output = workflows.execute("blast", {"db": ["synthetic"], "template": ">first\nACGT\n>second\nCGTA\n", "max_target_seqs": 1})
    result = output["results"][0]
    assert [q["name"] for q in output["queries"]] == ["first", "second"]
    assert result["hits"][0]["strand"] == "-" and result["hits"][0]["sstart"] == 13
    assert result["at_target_limit"] == ["query2"]


@pytest.mark.parametrize("params", [{"task": "bad"}, {"evalue": "nan"}, {"max_target_seqs": 0}, {"template": "ACG!"}])
def test_regular_blast_rejects_invalid_inputs_before_execution(monkeypatch, params):
    monkeypatch.setattr(sequence_tools.shutil, "which", lambda _: "blastn")
    monkeypatch.setattr(sequence_tools.subprocess, "run", lambda *a, **k: pytest.fail("must not execute"))
    with pytest.raises(ValueError):
        workflows.execute("blast", dict({"db": ["synthetic"], "template": "ACGT"}, **params))


def test_primer3_runs_without_database_and_never_claims_specificity(monkeypatch):
    from primerblast_oss import design
    def run(name, sequence, params):
        pair = design.PrimerPair(0, name, sequence[:20], revcomp(sequence[-20:]), 0, 20, 99, 20, 100, 60, 60, 50, 50)
        return [pair], "ok"
    monkeypatch.setattr(design, "design_primers", run)
    output = workflows.execute("primer3", {"template": "ACGT" * 25})
    pair = output["templates"][0]["pairs"][0]
    assert output["specificity_status"] == "not_evaluated"
    assert pair["sequence"] == "ACGT" * 25 and pair["specificity"] == {}
