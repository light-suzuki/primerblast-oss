"""Native commands support shell pipelines without requiring a GUI server."""
import csv
import io
import json
from types import SimpleNamespace

import pytest

from primerblast_oss import cli, design, sequence_tools, cli_sequence_tools
from primerblast_oss.genome import revcomp
from primerblast_oss.specificity import Amplicon


FWD = "ACGTTGCAAGTCCGATCGTA"
REV = "TGACCGTATGCTAGCTTACG"
REFERENCE = FWD + "GATTACA" * 9 + revcomp(REV)


def install_fake_blast(monkeypatch):
    calls = []
    monkeypatch.setattr(sequence_tools.shutil, "which", lambda name: name)
    def run(command, **kwargs):
        calls.append(command)
        query = command[command.index("-query") + 1]
        assert open(query).read().startswith(">query1\nACGT\n>query2\nCGTA\n")
        return SimpleNamespace(returncode=0, stderr=b"", stdout=
            b"query2\tchr1\t100\t4\t0\t0\t1\t4\t13\t10\t1e-6\t20\tCGTA\tCGTA\n")
    monkeypatch.setattr(sequence_tools.subprocess, "run", run)
    return calls


def install_fake_primer3(monkeypatch):
    calls = []
    def run(name, sequence, params, **kwargs):
        calls.append((name, sequence, params, kwargs))
        return [design.PrimerPair(0, name, sequence[:20], revcomp(sequence[-20:]),
                                 0, 20, 99, 20, 100, 60, 60, 50, 50)], "synthetic"
    monkeypatch.setattr(design, "design_primers", run)
    return calls


def install_fake_pcr(monkeypatch):
    calls = []
    def run(primers, db, **kwargs):
        calls.append((primers, db, kwargs))
        amps = [Amplicon("chr1", 1, 103, 103, "F", "R", 0, 0,
                         fwd_end3=20, rev_end3=84)] if primers == {"F": FWD, "R": REV} else []
        return {"db": db, "sites_per_primer": {name: 1 for name in primers},
                "n_products": len(amps), "products": amps, "search_completeness": "incomplete",
                "search_complete": False}
    monkeypatch.setattr(cli, "in_silico_pcr", run)
    return calls


def indexed_fasta(tmp_path):
    path = tmp_path / "synthetic.fa"
    path.write_bytes((">chr1\n" + REFERENCE + "\n").encode("ascii"))
    path.with_suffix(".fa.fai").write_text("chr1\t103\t6\t103\t104\n")
    return path


def test_blast_stdin_tsv_distinguishes_queries_and_references(monkeypatch, capsys):
    calls = install_fake_blast(monkeypatch)
    monkeypatch.setattr(cli.sys, "stdin", io.StringIO(">first\nACGT\n>second\nCGTA\n"))
    assert cli.main(["blast", "--query-fasta", "-", "--db", "a", "--db", "b",
                     "--task", "blastn-short", "--blastn-bin", "custom-blast",
                     "--num-threads", "3", "--evalue", "0.01"]) == 0
    rows = list(csv.DictReader(io.StringIO(capsys.readouterr().out), delimiter="\t"))
    assert [(r["db"], r["query_name"], r["sstart"], r["send"]) for r in rows] == [
        ("a", "second", "13", "10"), ("b", "second", "13", "10")]
    assert len(calls) == 2 and calls[0][0] == "custom-blast"
    assert all(r["strand"] == "-" for r in rows)
    assert calls[0][calls[0].index("-strand") + 1] == "both"
    assert calls[0][calls[0].index("-num_threads") + 1] == "3"
    assert calls[0][calls[0].index("-evalue") + 1] == "0.01"


def test_blast_single_json_file_and_limit_warning(monkeypatch, tmp_path, capsys):
    install_fake_blast(monkeypatch)
    path = tmp_path / "queries.fa"
    path.write_text(">first\nACGT\n>second\nCGTA\n")
    output = tmp_path / "blast.json"
    assert cli.main(["blast", "--query-fasta", str(path), "--db", "a", "--db", "b",
                     "--format", "json", "--max-target-seqs", "1", "--out", str(output)]) == 0
    data = json.loads(output.read_text())
    assert len(data["results"]) == 2 and data["results"][0]["hits"][0]["strand"] == "-"
    captured = capsys.readouterr()
    assert captured.out == "" and "target limit reached" in captured.err


@pytest.mark.parametrize("content", ["", "ACGT", ">empty\n", "> \nACGT\n", ">first\nACGT\n>bad\nACG!\n"])
def test_fasta_batches_rejected_before_any_tool_runs(monkeypatch, content, capsys):
    monkeypatch.setattr(cli.sys, "stdin", io.StringIO(content))
    monkeypatch.setattr(cli_sequence_tools, "nucleotide_search", lambda *a, **k: pytest.fail("must not run"))
    with pytest.raises(SystemExit) as error:
        cli.main(["blast", "--query-fasta", "-", "--db", "synthetic"])
    assert error.value.code == 2 and "Traceback" not in capsys.readouterr().err


def test_primer3_stdin_without_db_and_unique_fasta_exports(monkeypatch, tmp_path, capsys):
    calls = install_fake_primer3(monkeypatch)
    sequence = "ACGT" * 25
    monkeypatch.setattr(cli.sys, "stdin", io.StringIO(">same\n" + sequence + "\n>same\n" + sequence + "\n"))
    primers, products = tmp_path / "oligos.fa", tmp_path / "products.fa"
    assert cli.main(["primer3", "--template-fasta", "-", "--product-size", "80-110",
                     "--num-return", "2", "--target", "30,10", "--min-tm", "50",
                     "--primer3-bin", "custom-primer3", "--primers-out", str(primers),
                     "--products-fasta", str(products)]) == 0
    rows = list(csv.DictReader(io.StringIO(capsys.readouterr().out), delimiter="\t"))
    assert len(rows) == 2 and all(r["specificity_status"] == "not_evaluated" for r in rows)
    assert calls[0][2].product_size_ranges == [(80, 110)]
    assert calls[0][2].target == (30, 10) and calls[0][2].min_tm == 50
    assert calls[0][3] == {"primer3_bin": "custom-primer3"}
    records = design.read_fasta(str(primers))
    assert len(records) == len({name for name, _ in records}) == 4
    assert all(seq == sequence for _, seq in design.read_fasta(str(products)))


def test_primer3_multi_template_stdout_is_one_json_value(monkeypatch, capsys):
    install_fake_primer3(monkeypatch)
    assert cli.main(["primer3", "--template", ">a\n" + "ACGT" * 25 + "\n>b\n" + "TGCA" * 25,
                     "--format", "json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert len(data["templates"]) == 2 and data["specificity_status"] == "not_evaluated"
    assert "results" not in data and "databases" not in data


def test_output_path_collisions_fail_before_design(monkeypatch, tmp_path):
    monkeypatch.setattr(design, "design_primers", lambda *a, **k: pytest.fail("must not design"))
    path = str(tmp_path / "same")
    with pytest.raises(SystemExit) as error:
        cli.main(["primer3", "--template", "ACGT" * 25, "--out", path, "--products-fasta", path])
    assert error.value.code == 2 and not (tmp_path / "same").exists()


def test_check_cli_preserves_profiles_thermo_associations_and_literal_input(monkeypatch, tmp_path, capsys):
    calls = install_fake_pcr(monkeypatch)
    path = indexed_fasta(tmp_path)
    gff = tmp_path / "synthetic.gff3"
    gff.write_text("##gff-version 3\nchr1\tfixture\tgene\t10\t80\t.\t-\t.\tID=g1;Name=ExampleGene\n")
    fasta_out = tmp_path / "products.fa"
    monkeypatch.setattr(sequence_tools.shutil, "which", lambda _: None)
    assert cli.main(["check", "--forward", FWD.lower(), "--reverse", revcomp(REV),
                     "--input-orientation", "auto", "--db", "a", "--db", "b",
                     "--genome-fasta", str(path), "--gff3", str(gff),
                     "--max-total-mismatch", "2", "--specificity-profile", "ncbi",
                     "--exhaustive", "--evalue", "5", "--blastn-bin", "custom-blast",
                     "--no-thermo", "--format", "json", "--products-fasta", str(fasta_out)]) == 0
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert len(calls) == 8 and all(set(p) == {"F", "R"} for p, _, _ in calls)
    assert calls[0][2]["sp"].max_total_mismatch == 2 and calls[0][2]["sp"].max_target_seqs == 50000
    assert calls[0][2]["sp"].evalue == 5 and calls[0][2]["blastn_bin"] == "custom-blast"
    assert all(k["thermo_params"] is False for _, _, k in calls)
    assert calls[0][2]["genome"].fasta == str(path) and calls[1][2]["genome"] is None
    assert data["primers"] == {"F": FWD, "R": revcomp(REV)}
    assert data["input_assessments"][0]["original_search_complete"] is False
    corrected = [r for r in data["results"] if r["reverse_complemented_inputs"] == ["R"]]
    a, b = [r["products"][0] for r in corrected]
    assert a["input_evidence"]["sequence_status"] == "reverse_complement_candidate"
    assert a["sequence"] == REFERENCE and a["annotations"]["genes"][0]["name"] == "ExampleGene"
    assert b["sequence_status"] == "unavailable" and b["annotations"]["status"] == "not_provided"
    assert "partial" in captured.err
    records = design.read_fasta(str(fasta_out))
    assert len(records) == 1 and records[0][1] == REFERENCE
    assert "reverse_complement_candidate" in records[0][0]


def test_check_gene_only_does_not_extract_sequence(monkeypatch, tmp_path, capsys):
    install_fake_pcr(monkeypatch)
    monkeypatch.setattr(cli, "_thermo_setup", lambda _: ({}, False, True))
    monkeypatch.setattr(sequence_tools, "product_sequence", lambda *a: pytest.fail("must not extract"))
    gff = tmp_path / "synthetic.gff3"
    gff.write_text("##gff-version 3\nchr1\ttest\tgene\t10\t80\t.\t+\t.\tID=g1;Name=ExampleGene\n")
    assert cli.main(["check", "--forward", FWD, "--reverse", REV, "--db", "a", "--db", "b",
                     "--db-gff3", "b=" + str(gff), "--format", "tsv"]) == 0
    captured = capsys.readouterr()
    a, b = list(csv.DictReader(io.StringIO(captured.out), delimiter="\t"))
    assert a["genes"] == "" and a["annotation_status"] == "not_provided"
    assert b["genes"] == "ExampleGene" and b["annotation_status"] == "loaded"
    assert a["hypothesis"] == b["hypothesis"] == "1"
    assert a["original_search_complete"] == "False" and a["as_supplied_locus_observed"] == "True"
    assert a["left_direction"] == "+:right" and a["sequence_status"] == "as_supplied"
    assert "search=incomplete" in captured.err


def test_check_stdin_normalizes_sequences_and_text_shows_input_evidence(monkeypatch, capsys):
    calls = install_fake_pcr(monkeypatch)
    monkeypatch.setattr(cli.sys, "stdin", io.StringIO(">F\n" + FWD.lower() + "\n>R\n" + REV + "\n"))
    assert cli.main(["check", "--primers-fasta", "-", "--db", "a", "--no-thermo",
                     "--show-sequence-forms"]) == 0
    text = capsys.readouterr().out
    assert len(calls) == 1 and calls[0][0] == {"F": FWD, "R": REV}
    assert "AS SUPPLIED" in text and "original search complete=False" in text
    assert "as_supplied; labels=as_labeled" in text and "complement 3'-5'" in text


@pytest.mark.parametrize("extra", [["--primer", "F=ACGT"], ["--reverse", "ACG!"],
                                  ["--input-orientation", "auto", "--primer", "B=CGTA", "--primer", "C=GGAA"],
                                  ["--db-gff3", "other=synthetic.gff3"]])
def test_check_input_errors_fail_before_any_pcr(monkeypatch, extra, capsys):
    monkeypatch.setattr(cli, "in_silico_pcr", lambda *a, **k: pytest.fail("must not search"))
    with pytest.raises(SystemExit) as error:
        cli.main(["check", "--forward", FWD, "--db", "a"] + extra)
    assert error.value.code == 2 and "Traceback" not in capsys.readouterr().err


def test_missing_file_has_cli_error_without_traceback(tmp_path, capsys):
    with pytest.raises(SystemExit) as error:
        cli.main(["primer3", "--template-fasta", str(tmp_path / "missing.fa")])
    assert error.value.code == 2 and "Traceback" not in capsys.readouterr().err


def test_no_product_tsv_keeps_completeness_diagnostics(monkeypatch, capsys):
    install_fake_pcr(monkeypatch)
    assert cli.main(["check", "--forward", FWD, "--reverse", revcomp(REV), "--db", "a",
                     "--no-thermo", "--format", "tsv"]) == 0
    captured = capsys.readouterr()
    assert captured.out.startswith("db\thypothesis\t") and len(captured.out.splitlines()) == 1
    assert "0 product(s)" in captured.err and "search=incomplete" in captured.err
