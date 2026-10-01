"""Regressions for gap coordinates, rejected sites and native batch JSON."""
import json
from types import SimpleNamespace

import pytest

from primerblast_oss import cli, report, specificity as S, workflows
from primerblast_oss.pipeline import PipelineResult


@pytest.mark.parametrize("strand,start,end,expected", [
    ("plus", 100, 122, (100, 122)), ("minus", 122, 100, (122, 100)),
])
def test_hsp_indels_preserve_reference_footprint(strand, start, end, expected):
    fields = ["primer", "chr1", "87", "23", "0", "1", "1", "20",
              str(start), str(end), "1e-3", "20", strand,
              "AAAAAAAAAA---AAAAAAAAAA", "AAAAAAAAAACCCAAAAAAAAAA", "20"]
    site = S._hit_to_site(fields, "F", S.SpecParams(max_total_mismatch=3))
    assert site is not None
    assert (site.end5, site.end3) == expected
    mate = S.PrimingSite("R", "chr1", "-", 181, 0, 0, plen=20)
    if strand == "plus":
        products = S.enumerate_amplicons([site, mate], S.SpecParams())
        assert [(a.start, a.end, a.size) for a in products] == [(100, 200, 101)]


def test_thermo_rejected_sites_remain_visible_without_products(monkeypatch):
    from primerblast_oss import thermo
    sites = [S.PrimingSite("F", "chr1", "+", 20, 0, 0),
             S.PrimingSite("R", "chr1", "-", 81, 0, 0)]
    monkeypatch.setattr(S, "_detect_blastn", lambda _: "blastn")
    monkeypatch.setattr(S, "screen_primers_with_stats", lambda *a: (sites, {}))
    monkeypatch.setattr(thermo, "available", lambda: True)
    monkeypatch.setattr(thermo, "evaluate", lambda *a: SimpleNamespace(tm=10, end3_dg=0, viable=False))
    genome = SimpleNamespace(fetch=lambda *a: "T" * 20)
    raw = S.in_silico_pcr({"F": "A" * 20, "R": "A" * 20}, "synthetic", genome=genome)
    data = report.insilico_to_dict([raw], raw["primers"])["results"][0]
    assert data["products"] == []
    assert data["binding_site_counts"] == {"F": {"+": 1, "-": 0}, "R": {"+": 0, "-": 1}}
    assert len(data["binding_sites"]) == 2
    assert all(s["thermo_viable"] is False for s in data["binding_sites"])


@pytest.mark.parametrize("mode", ["design", "tile"])
@pytest.mark.parametrize("names", [["a"], ["a", "b"]])
def test_native_multiple_templates_emit_one_json_document(mode, names, tmp_path, monkeypatch, capsys):
    fasta = tmp_path / "templates.fa"
    fasta.write_text("".join(">" + name + "\nACGT\n" for name in names), encoding="utf-8")
    monkeypatch.setattr(cli, "_thermo_setup", lambda _: ({}, False, True))
    monkeypatch.setattr(cli, "run_pipeline", lambda name, seq, db, **k:
                        PipelineResult(name, len(seq), [], "synthetic", db, {}))
    monkeypatch.setattr(cli, "design_tiling", lambda *a, **k: [])
    assert cli.main([mode, "--template-fasta", str(fasta), "--db", "synthetic", "--format", "json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["mode"] == mode
    if len(names) == 1:
        assert data["template_id"] == "a" and "templates" not in data
    else:
        assert [t["template_id"] for t in data["templates"]] == names


def test_gui_check_reports_thermo_coverage(monkeypatch):
    monkeypatch.setattr(workflows, "in_silico_pcr", lambda primers, db, **k: {
        "db": db, "sites_per_primer": {"F": 0}, "products": [], "n_products": 0,
        "search_complete": True, "search_completeness": "complete"})
    data = workflows.execute("check", {"forward": "ACGT", "db": ["synthetic"]})
    assert data["results"][0]["thermo_status"] == "skipped_no_associated_genome"


def test_failed_thermo_calculation_is_unresolved_not_absent(monkeypatch):
    from primerblast_oss import thermo
    from primerblast_oss.pipeline import thermo_metadata
    site = S.PrimingSite("F", "chr1", "+", 20, 0, 0)
    monkeypatch.setattr(thermo, "available", lambda: True)
    def fail(*args):
        raise ValueError("unsupported oligo")
    monkeypatch.setattr(thermo, "evaluate", fail)
    genome = SimpleNamespace(fetch=lambda *a: "T" * 20)
    kept, viable, stats = S.annotate_thermo([site], {"F": "A" * 20}, genome)
    assert kept == [site] and viable == {}
    assert site.thermo_viable is None and stats["unresolved_per_primer"] == {"F": 1}
    assert thermo_metadata(genome, None, True, "explicit_db_mapping", stats)["thermo_status"] == "failed_no_resolvable_sites"
