import json

import pytest

from primerblast_oss.gff3 import parse_gff3
from primerblast_oss.regions import resolve_gene
from primerblast_oss.webapp.references import reference_catalog, infer_gene_id_format
from primerblast_oss.webapp.server import _find_gene_seqid


def test_gene_typography_and_registered_aliases(tmp_path):
    path = tmp_path / "genes.gff3"
    path.write_text(
        "chr1\tx\tgene\t10\t80\t.\t+\t.\tID=gene:AT1G01010;Name=NAC001;Alias=ANAC001;gene_id=old_id\n",
        encoding="utf-8",
    )
    gff = parse_gff3(str(path))
    for query in ("　ＡＴ１Ｇ０１０１０　", "at1g01010", "Ｇｅｎｅ：ＡＴ１Ｇ０１０１０", "gene : AT1G01010", "anac001", "old_id"):
        assert gff.gene(query).id == "gene:AT1G01010"
        assert _find_gene_seqid(str(path), query) == "chr1"
        assert resolve_gene(str(path), query, feature="gene").start == 10
    assert gff.gene("AT1G01010.1") is None
    assert gff.gene("AT1G0101") is None


def test_ambiguous_alias_is_not_hidden_by_chromosome_prescan(tmp_path):
    path = tmp_path / "genes.gff3"
    path.write_text(
        "chr1\tx\tgene\t10\t80\t.\t+\t.\tID=g1;Name=Same\n"
        "chr2\tx\tgene\t10\t80\t.\t+\t.\tID=g2;Name=SAME;Alias=g1\n",
        encoding="utf-8",
    )
    assert _find_gene_seqid(str(path), "same") is None
    with pytest.raises(ValueError, match="Ambiguous.*g1, g2"):
        resolve_gene(str(path), "same", feature="gene")
    assert parse_gff3(str(path)).gene("g1").seqid == "chr1"
    assert parse_gff3(str(path)).gene("　Ｇ１　").seqid == "chr1"


def test_reference_catalog_only_pairs_unambiguous_siblings(tmp_path):
    prefix = tmp_path / "ref"
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\nACGT\n")
    (tmp_path / "ref.fa.fai").write_text("chr1\t4\t6\t4\t5\n")
    (tmp_path / "ref.nin").touch()
    (tmp_path / "ref.gff3").touch()
    databases = [{"name": "ref", "path": str(prefix)}]
    empty_config = tmp_path / "not-present.json"
    profiles = reference_catalog(databases, empty_config)["references"]
    assert profiles[0]["genome"] == str(fasta)
    assert profiles[0]["gff3"] == str(tmp_path / "ref.gff3")
    assert profiles[0]["available"]
    (tmp_path / "ref.gff3.gz").touch()
    assert reference_catalog(databases, empty_config)["references"][0]["gff3"] == ""
    (tmp_path / "ref.fna").touch()
    (tmp_path / "ref.fna.fai").touch()
    assert reference_catalog(databases, empty_config)["references"] == []


def test_explicit_reference_availability_and_bad_config(tmp_path):
    config = tmp_path / "references.json"
    config.write_text(json.dumps([{"name": "missing", "genome": str(tmp_path / "missing.fa")}]))
    profile = reference_catalog([], config)["references"][0]
    assert not profile["available"]
    assert profile["missing"] == ["FASTA", "FASTA index (.fai)"]
    config.write_text("{}")
    assert reference_catalog([], config)["warnings"]


def test_gene_format_follows_annotation_and_index(tmp_path):
    genome = tmp_path / "public-example.fa"
    (tmp_path / "public-example.fa.fai").write_text("chr1\t100\t0\t60\t61\n7LG7\t100\t0\t60\t61\n")
    annotation = tmp_path / "example.gff3"
    annotation.write_text("chr1\tx\tgene\t1\t30\t.\t+\t.\tID=gene:Ref.1g00050\n")
    format = infer_gene_id_format(genome, annotation)
    assert format == {"prefix": "Ref.", "separator": "g", "digits": 5,
                      "chromosomes": ["1", "7"], "source": "annotation_sample"}
    annotation.write_text("chr1\tx\tgene\t1\t30\t.\t+\t.\tID=Ref.1g00050\n"
                          "chr1\tx\tgene\t31\t60\t.\t+\t.\tID=Different.1g00060\n")
    assert infer_gene_id_format(genome, annotation) is None
