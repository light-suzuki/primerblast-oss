from primerblast_oss.annotations import template_annotations
from primerblast_oss.gff3 import parse_gff3
from primerblast_oss.webapp.server import _find_gene_seqid


def test_annotations_preserve_isoforms_and_genomic_coordinates(tmp_path):
    path = tmp_path / "genes.gff3"
    path.write_text(
        "chr1\ttest\tgene\t21\t90\t.\t-\t.\tID=g1;Name=Gene1\n"
        "chr1\ttest\tmRNA\t21\t90\t.\t-\t.\tID=t1;Parent=g1\n"
        "chr1\ttest\texon\t21\t30\t.\t-\t.\tParent=t1\n"
        "chr1\ttest\texon\t60\t90\t.\t-\t.\tParent=t1\n"
        "chr1\ttest\tCDS\t65\t85\t.\t-\t0\tParent=t1\n"
        "chr1\ttest\tmRNA\t21\t90\t.\t-\t.\tID=t2;Parent=g1\n"
        "chr1\ttest\texon\t21\t45\t.\t-\t.\tParent=t2\n",
        encoding="utf-8",
    )
    context = {"chrom": "chr1", "start": 15, "end": 70, "anchor": 70, "strand": "-"}
    result = template_annotations(str(path), context)
    gene = result["genes"][0]
    assert gene["start"] == 21 and gene["end"] == 90
    assert gene["strand"] == "-"
    assert [tx["id"] for tx in gene["transcripts"]] == ["t1", "t2"]
    assert gene["transcripts"][0]["segments"][-1] == {"type": "cds", "start": 65, "end": 85}
    assert template_annotations(str(path), {**context, "start": 91, "end": 100})["genes"] == []
    assert template_annotations(str(path), {**context, "chrom": "other"})["status"] == "seqid_not_found"
    assert template_annotations(None, context)["status"] == "not_provided"


def test_ensembl_gene_alias_and_prescan_do_not_match_prefixes(tmp_path):
    path = tmp_path / "ensembl.gff3"
    path.write_text(
        "wrong\ttest\tgene\t1\t10\t.\t+\t.\tID=gene:AT1G010100\n"
        "1\ttest\tgene\t21\t90\t.\t+\t.\tID=gene:AT1G01010;Name=NAC001\n",
        encoding="utf-8",
    )
    assert _find_gene_seqid(str(path), "AT1G01010") == "1"
    assert parse_gff3(str(path)).gene("AT1G01010").id == "gene:AT1G01010"
    assert _find_gene_seqid(str(path), "AT1G0101") is None
