"""Synthetic regressions for nominated loci and duplicate HSP evidence.

The repeat oracle is exhaustive substring matching, independent of both BLAST
and the fitting aligner. All sequences and coordinates in this file are synthetic.
"""
from dataclasses import replace

import pytest

from primerblast_oss.genome import revcomp
from primerblast_oss.specificity import (
    SEARCH_COMPLETE, SEARCH_REPEAT_LIMITED, PrimingSite, SpecParams,
    _priming_sites_from_output, _realign_hit_to_site, enumerate_amplicons,
)


class MemoryGenome:
    def __init__(self, sequence):
        self.sequence = sequence

    def length(self, _subject):
        return len(self.sequence)

    def fetch(self, _subject, start, end, strand="+"):
        sequence = self.sequence[start - 1:end]
        return revcomp(sequence) if strand == "-" else sequence


def hsp(primer, start, end, strand="+", qstart=1, qend=None):
    qend = len(primer) if qend is None else qend
    aligned = primer[qstart - 1:qend]
    return list(map(str, [
        "query", "synthetic", 100, len(aligned), 0, 0, qstart, qend,
        start, end, "1e-8", 40, "plus" if strand == "+" else "minus",
        aligned, aligned, len(primer),
    ]))


def parse(primer, rows, genome=None, sp=None, name="F"):
    return _priming_sites_from_output(
        primer, name, "\n".join("\t".join(row) for row in rows),
        sp or SpecParams(), genome,
    )


@pytest.mark.parametrize("reverse_reference", [False, True])
def test_every_exact_repeat_locus_and_product_is_preserved(reverse_reference):
    primers = {"F": "AC" * 10, "R": "TGCACGTAGGCTTACAGGTC"}
    reference = "N" * 100 + "AC" * 15 + "N" * 90 + revcomp(primers["R"]) + "N" * 100
    if reverse_reference:
        reference = revcomp(reference)
    genome = MemoryGenome(reference)
    sp = SpecParams(max_total_mismatch=0, max_3prime_mismatch=0,
                    min_product=50, max_product=300)
    expected_sites, sites = [], []
    for name, primer in primers.items():
        rows = []
        for strand in ("+", "-"):
            query = primer if strand == "+" else revcomp(primer)
            for index in range(len(reference) - len(query) + 1):
                if reference[index:index + len(query)] != query:
                    continue
                low, high = index + 1, index + len(query)
                end5, end3 = (low, high) if strand == "+" else (high, low)
                expected_sites.append((name, strand, end5, end3))
                rows.append(hsp(primer, end5, end3, strand))
        found, stats = parse(primer, rows, genome, sp, name)
        assert stats.completeness == SEARCH_COMPLETE
        assert stats.raw_blast_hits == stats.priming_sites == stats.realigned_sites == len(rows)
        sites.extend(found)
    assert sorted((s.primer, s.strand, s.end5, s.end3) for s in sites) == sorted(expected_sites)
    expected_products = sorted(
        (left[2], right[2], right[2] - left[2] + 1, left[0], right[0])
        for left in expected_sites for right in expected_sites
        if left[1] == "+" and right[1] == "-" and left[3] <= right[3]
        and sp.min_product <= right[2] - left[2] + 1 <= sp.max_product
    )
    products = enumerate_amplicons(sites, sp)
    assert sorted((a.start, a.end, a.size, a.fwd_primer, a.rev_primer) for a in products) == expected_products
    assert sorted(a.size for a in products) == [130, 132, 134, 136, 138, 140]


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("qstart,qend", [(1, 20), (4, 20), (1, 17), (4, 17)])
def test_exact_projected_repeat_site_survives_partial_hsp(strand, qstart, qend):
    primer = "AC" * 10
    reference = "N" * 100 + ("AC" * 15 if strand == "+" else "GT" * 15) + "N" * 100
    end5, end3 = (105, 124) if strand == "+" else (126, 107)
    step = 1 if strand == "+" else -1
    fields = hsp(primer, end5 + step * (qstart - 1),
                 end3 - step * (20 - qend), strand, qstart, qend)
    site, resolved = _realign_hit_to_site(fields, "F", primer, SpecParams(), MemoryGenome(reference))
    assert resolved and site is not None
    assert (site.end5, site.end3) == (end5, end3)
    assert site.aligned_query == site.aligned_target == primer
    assert site.total_mismatch == 0


@pytest.mark.parametrize("strand", ["+", "-"])
def test_exact_hsp_claim_must_be_verified_against_genome(strand):
    primer = "AC" * 10
    reference = "N" * 100 + ("AC" * 11 if strand == "+" else "GT" * 11) + "N" * 100
    # The claimed perfect site starts one base out of phase. Preserve the
    # existing fitting-alignment fallback rather than trusting HSP metadata.
    end5, end3 = (102, 121) if strand == "+" else (121, 102)
    site, resolved = _realign_hit_to_site(
        hsp(primer, end5, end3, strand), "F", primer, SpecParams(), MemoryGenome(reference))
    assert resolved and site is not None
    assert (site.end5, site.end3) == ((101, 120) if strand == "+" else (122, 103))
    assert site.total_mismatch == 0


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("insertion", [False, True])
def test_indel_realign_keeps_genomic_footprint(strand, insertion):
    primer = "GCTAGCTACGATCGTACGTA"
    target = primer[:10] + ("A" + primer[10:] if insertion else primer[11:])
    genome = MemoryGenome("N" * 100 + (target if strand == "+" else revcomp(target)) + "N" * 100)
    low, high = 101, 100 + len(target)
    end5, end3 = (low, high) if strand == "+" else (high, low)
    site, resolved = _realign_hit_to_site(
        hsp(primer, end5, end3, strand), "F", primer, SpecParams(), genome)
    assert resolved and site is not None
    assert (site.end5, site.end3) == (end5, end3)
    assert site.total_mismatch == 1
    assert "-" in (site.aligned_query if insertion else site.aligned_target)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_different_hsps_converging_on_one_site_are_counted_once(strand):
    primer = "GCTAGCTACGATCGTACGTA"
    genome = MemoryGenome("N" * 100 + (primer if strand == "+" else revcomp(primer)) + "N" * 100)
    end5, end3 = (101, 120) if strand == "+" else (120, 101)
    step = 1 if strand == "+" else -1
    rows = [hsp(primer, end5, end3, strand),
            hsp(primer, end5 + step * 3, end3, strand, qstart=4)]
    sites, stats = parse(primer, rows, genome, SpecParams(high_copy_hit_threshold=2))
    assert len(sites) == stats.priming_sites == stats.realigned_sites == 1
    assert stats.raw_blast_hits == 2
    assert stats.high_copy and stats.completeness == SEARCH_REPEAT_LIMITED


def test_hsp_only_duplicates_collapse_but_alignment_alternatives_remain():
    primer = "ACGTAACGTC"
    row = hsp(primer, 101, 111)
    row[13:15] = ["ACGTA-ACGTC", "ACGTAAACGTC"]
    alternative = row[:]
    alternative[13] = "ACGT-AACGTC"
    sites, stats = parse(primer, [row, row, alternative])
    assert stats.raw_blast_hits == 3
    assert len(sites) == stats.priming_sites == 2
    assert sites[0].aligned_query != sites[1].aligned_query
    assert (sites[0].end5, sites[0].end3) == (sites[1].end5, sites[1].end3)


def test_pairing_deduplicates_sites_without_collapsing_distinct_evidence():
    forward = PrimingSite("R", "synthetic", "+", 120, 1, 0,
                          mapped_end5=101, aligned_query="AC-GT", aligned_target="ACAGT")
    reverse = PrimingSite("R", "synthetic", "-", 221, 0, 0, mapped_end5=240)
    alternative_end3 = replace(forward, end3=121)
    alternative_alignment = replace(forward, aligned_query="A-CGT")
    sites = [forward, replace(forward), reverse, replace(reverse),
             alternative_end3, alternative_alignment]
    products = enumerate_amplicons(sites, SpecParams())
    assert len(products) == 3
    assert {(a.start, a.end, a.size, a.orientation) for a in products} == {(101, 240, 140, "R/R")}
    assert sorted(a.fwd_end3 for a in products) == [120, 120, 121]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_imperfect_repeat_keeps_existing_fitting_alignment_behavior(strand):
    primer = "AT" + "AC" * 9
    reference = "N" * 100 + ("AC" * 15 if strand == "+" else "GT" * 15) + "N" * 100
    end5, end3 = (105, 124) if strand == "+" else (126, 107)
    site, resolved = _realign_hit_to_site(
        hsp(primer, end5, end3, strand), "F", primer, SpecParams(), MemoryGenome(reference))
    assert resolved and site is not None
    # No perfect match at the nominated locus: keep the pre-existing fitting
    # score/tie policy, including its ability to choose an adjacent placement.
    assert (site.end5, site.end3) == ((103, 122) if strand == "+" else (128, 109))
    assert site.total_mismatch == 1


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("at_start", [False, True])
@pytest.mark.parametrize("partial", [False, True])
def test_exact_site_at_contig_edge(strand, at_start, partial):
    primer = "AC" * 10
    reference = "AC" * 20 if strand == "+" else "GT" * 20
    low, high = (1, 20) if at_start else (21, 40)
    end5, end3 = (low, high) if strand == "+" else (high, low)
    qstart, qend = (4, 17) if partial else (1, 20)
    step = 1 if strand == "+" else -1
    fields = hsp(primer, end5 + step * (qstart - 1),
                 end3 - step * (20 - qend), strand, qstart, qend)
    site, resolved = _realign_hit_to_site(fields, "F", primer, SpecParams(), MemoryGenome(reference))
    assert resolved and site is not None
    assert (site.end5, site.end3) == (end5, end3)
