"""GUI adapters preserve CLI evidence and reject malformed input."""
import json
from pathlib import Path

import pytest

from primerblast_oss.webapp import server


def test_associations_reject_unselected_database():
    with pytest.raises(ValueError, match="selected databases"):
        server._associated_genomes({"db": ["selected"], "db_genomes": "other=x.fa"})


def test_associations_require_explicit_mapping_syntax():
    with pytest.raises(ValueError, match="DB=FASTA"):
        server._associated_genomes({"db": ["selected"], "db_genomes": "x.fa"})


def test_associations_keep_other_assemblies_unmapped(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGTACGT\n")
    Path(str(fasta) + ".fai").write_text("chr1\t8\t6\t8\t9\n")
    mappings = server._associated_genomes({"db": ["a", "b"], "genome": str(fasta)})
    assert set(mappings) == {"a"}


def test_sequence_adapter_exports_cli_plan(monkeypatch):
    from primerblast_oss import cli

    def sequence(args):
        assert args.template == "ACGT" * 50
        assert args.amplicon_size == "100-150"
        assert args.m13_tails is True
        Path(args.out).write_text(json.dumps({"mode": "sequence", "plans": [{
            "template_id": "test", "amplicons": [], "coverage": {
                "full_coverage": False, "gaps": [[0, 199]]}}]}))

    monkeypatch.setattr(cli, "_cmd_sequence", sequence)
    result = server._run_sequence({"db": ["fixture"], "template": "ACGT" * 50,
                                   "amplicon_size": "100-150", "m13_tails": True})
    assert result["plans"][0]["coverage"]["full_coverage"] is False
    assert "order" in result["exports"]


def test_sequence_invalid_number_is_regular_job_error():
    with pytest.raises(ValueError, match="Invalid sequencing parameters"):
        server._run_sequence({"db": ["fixture"], "template": "ACGT",
                              "overlap": "invalid"})
