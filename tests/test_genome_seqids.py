"""Synthetic regressions for BLAST accession wrappers and FASTA identity.

The runtime bug was observed with a public NCBI accession, but all sequences
and coordinates below are synthetic and require no tools or downloaded data.
"""
from pathlib import Path
from types import SimpleNamespace

import pytest

from primerblast_oss.genome import Genome, revcomp
from primerblast_oss.assay import reclassify_by_anchor
from primerblast_oss.regions import GenomicRegion, extract_template
from primerblast_oss.sequence_tools import product_sequence
from primerblast_oss.specificity import (
    SEARCH_COMPLETE, SEARCH_POSSIBLY_TRUNCATED, SpecParams,
    _priming_sites_from_output, enumerate_amplicons,
)


ACCESSION = "NC_123456.1"


def make_genome(tmp_path, records):
    fasta = tmp_path / "reference.fa"
    payload = bytearray()
    index = []
    for name, sequence in records:
        payload.extend((">" + name + "\n").encode())
        index.append("%s\t%d\t%d\t%d\t%d\n" % (
            name, len(sequence), len(payload), len(sequence), len(sequence) + 1))
        payload.extend((sequence + "\n").encode())
    fasta.write_bytes(payload)
    Path(str(fasta) + ".fai").write_text("".join(index))
    return Genome(str(fasta))


@pytest.mark.parametrize("prefix", ["ref", "gb", "emb", "dbj"])
@pytest.mark.parametrize("wrapped_fasta", [False, True])
def test_standard_accession_wrappers_resolve_uniquely(tmp_path, prefix, wrapped_fasta):
    wrapped = "%s|%s|locus" % (prefix, ACCESSION)
    key, query = (wrapped, ACCESSION) if wrapped_fasta else (ACCESSION, wrapped)
    genome = make_genome(tmp_path, [(key, "ACGTTGCA")])
    assert genome.resolve_name(query) == key
    assert query in genome
    assert genome.length(query) == 8
    assert genome.fetch(query, 2, 6) == "CGTTG"
    assert genome.fetch(query, 2, 6, "-") == "CAACG"
    assert genome.chroms() == [key]


def test_exact_key_wins_without_hiding_ambiguous_aliases(tmp_path):
    wrapped = "ref|%s|" % ACCESSION
    genome = make_genome(tmp_path, [(ACCESSION, "AAAA"), (wrapped, "CCCC")])
    assert genome.fetch(ACCESSION, 1, 4) == "AAAA"
    assert genome.fetch(wrapped, 1, 4) == "CCCC"
    assert genome.resolve_name(wrapped) == wrapped
    alternate = "gb|%s|" % ACCESSION
    assert alternate not in genome
    with pytest.raises(KeyError, match="[Aa]mbiguous"):
        genome.fetch(alternate, 1, 4)


def test_colliding_wrapped_keys_do_not_resolve_bare_accession(tmp_path):
    genome = make_genome(tmp_path, [
        ("ref|%s|" % ACCESSION, "AAAA"),
        ("gb|%s|locus" % ACCESSION, "CCCC"),
    ])
    assert ACCESSION not in genome
    with pytest.raises(KeyError, match="[Aa]mbiguous"):
        genome.length(ACCESSION)


@pytest.mark.parametrize("query", [
    "NC_123456", "ref|NC_123456.2|", "ref|nc_123456.1|",
    "lcl|NC_123456.1", "custom|NC_123456.1|", "gi|42|ref|NC_123456.1|",
    "ref|NC_123456.1", "ref|NC_123456.1|locus|extra", "chrNC_123456.1",
])
def test_versions_and_unsupported_aliases_are_not_guessed(tmp_path, query):
    genome = make_genome(tmp_path, [(ACCESSION, "ACGT")])
    assert query not in genome
    with pytest.raises(KeyError):
        genome.fetch(query, 1, 4)


def hsp(primer, subject, start, end):
    return "\t".join(map(str, [
        "query", subject, 100, len(primer), 0, 0, 1, len(primer),
        start, end, "1e-8", 40, "plus" if start < end else "minus",
        primer, primer, len(primer),
    ]))


def test_realign_product_coordinates_and_sequence_use_fasta_key(tmp_path):
    forward = "GCTAGCTACGATCGTACGTA"
    reverse = "TGCACGTAGGCTTACAGGTC"
    sequence = "N" * 100 + forward + "A" * 60 + revcomp(reverse) + "N" * 100
    genome = make_genome(tmp_path, [(ACCESSION, sequence)])
    sp = SpecParams(min_product=40, max_product=200)
    sites = []
    for name, primer, start, end in (
            ("F", forward, 101, 120), ("R", reverse, 200, 181)):
        found, stats = _priming_sites_from_output(
            primer, name, hsp(primer, "ref|%s|" % ACCESSION, start, end), sp, genome)
        assert stats.completeness == SEARCH_COMPLETE
        assert stats.realignment_failures == 0
        assert len(found) == 1
        assert (found[0].subject, found[0].end5, found[0].end3) == (ACCESSION, start, end)
        sites.extend(found)
    products = enumerate_amplicons(sites, sp)
    assert len(products) == 1
    anchored = reclassify_by_anchor(
        {"on_target": [], "off_target": products, "search_completeness": SEARCH_COMPLETE},
        template=extract_template(genome, GenomicRegion(ACCESSION, 1, len(sequence))),
        pair=SimpleNamespace(left_start=100, right_start=199),
    )
    assert anchored["intended_status"] == "unique"
    assert anchored["specific"] is True
    product = products[0].__dict__.copy()
    assert (product["subject"], product["start"], product["end"], product["size"]) == (
        ACCESSION, 101, 200, 100)
    product_sequence(product, "unused-db", genome)
    assert product["sequence_status"] == "available"
    assert product["sequence"] == sequence[100:200]


def test_without_genome_original_blast_subject_is_preserved():
    primer = "GCTAGCTACGATCGTACGTA"
    subject = "ref|%s|" % ACCESSION
    sites, stats = _priming_sites_from_output(
        primer, "F", hsp(primer, subject, 101, 120), SpecParams())
    assert stats.completeness == SEARCH_COMPLETE
    assert len(sites) == 1 and sites[0].subject == subject


def test_ambiguous_subject_remains_unresolved_not_negative_evidence(tmp_path):
    primer = "GCTAGCTACGATCGTACGTA"
    genome = make_genome(tmp_path, [
        ("ref|%s|" % ACCESSION, primer + "A" * 20),
        ("gb|%s|" % ACCESSION, primer + "C" * 20),
    ])
    sites, stats = _priming_sites_from_output(
        primer, "F", hsp(primer, ACCESSION, 1, 20), SpecParams(), genome)
    assert sites == []
    assert stats.realignment_failures == 1
    assert stats.completeness == SEARCH_POSSIBLY_TRUNCATED
