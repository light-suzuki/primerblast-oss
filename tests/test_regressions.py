from primerblast_oss.assay import expected_amplicon_from_design, reclassify_by_anchor
from primerblast_oss.design import DesignParams, PrimerPair, _build_boulder, clean_sequence
from primerblast_oss.regions import GenomicRegion, Template
from primerblast_oss.risk import assess_risk
from primerblast_oss.specificity import Amplicon, SpecParams, _hit_to_site
from primerblast_oss.variants import footprints_from_amplicon, snps_under_primers


class Variant:
    def __init__(self, chrom, pos, ref, alt):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt

    @property
    def end(self):
        return self.pos + len(self.ref) - 1


def pair():
    return PrimerPair(
        index=0, template_id="target", forward="A" * 20, reverse="C" * 21,
        left_start=100, left_len=20, right_start=468, right_len=21,
        product_size=369, tm_f=60, tm_r=60, gc_f=50, gc_r=50,
    )


def template(strand="+"):
    region = GenomicRegion("chr1", 5000, 6500, strand, "target")
    return Template(
        id="target", seq="A" * 1501, region=region,
        ext_start=5000, ext_end=6500,
        anchor_coord=5000 if strand == "+" else 6500,
        anchor_strand=strand,
    )


def test_multibase_gap_counts_every_column():
    qseq = "AAAAAAAAAA---AAAAAAAAAA"
    sseq = "AAAAAAAAAACCCAAAAAAAAAA"
    fields = ["primer", "chr1", "87", "23", "0", "1", "1", "20",
              "100", "122", "1e-3", "20", "plus", qseq, sseq, "20"]
    assert _hit_to_site(fields, "F", SpecParams(max_total_mismatch=1)) is None
    site = _hit_to_site(fields, "F", SpecParams(max_total_mismatch=3))
    assert site is not None and site.total_mismatch == 3


def test_ambiguity_preserves_coordinates_and_is_excluded():
    raw = "AAA CCC\nNNRYS--GGG"
    seq = clean_sequence(raw)
    assert seq == "AAACCCNNNNNNNGGG"
    assert len(seq) == len("AAACCCNNRYS--GGG")
    boulder = _build_boulder("x", seq, DesignParams(target=(13, 1), excluded_regions=((1, 2),)))
    assert f"SEQUENCE_TEMPLATE={seq}" in boulder
    assert "SEQUENCE_TARGET=13,1" in boulder
    assert "SEQUENCE_EXCLUDED_REGION=1,2 6,7" in boulder


def test_rf_footprints_use_actual_primer_lengths_and_positions():
    amp = Amplicon("chr1", 100, 500, 401, "R", "F", 0, 0, on_target=True)
    fps = footprints_from_amplicon(amp, len_f=20, len_r=21)
    fwd, rev = fps
    assert (fwd.primer, fwd.start, fwd.end, fwd.three_prime_coord, fwd.strand) == (
        "F", 481, 500, 481, "-"
    )
    assert (rev.primer, rev.start, rev.end, rev.three_prime_coord, rev.strand) == (
        "R", 100, 120, 120, "+"
    )
    hits = snps_under_primers(fps, [Variant("chr1", 482, "A", ["G"])])
    assert len(hits) == 1 and hits[0].primer == "F" and hits[0].in_3prime_5bp


def test_indel_reference_span_overlaps_primer():
    amp = Amplicon("chr1", 100, 500, 401, "F", "R", 0, 0, on_target=True)
    fps = footprints_from_amplicon(amp, 20, 20)
    deletion = Variant("chr1", 95, "A" * 10, ["A"])
    terminal_mnp = Variant("chr1", 116, "AAAAA", ["CCCCC"])
    hits = snps_under_primers(fps, [deletion, terminal_mnp])
    by_pos = {hit.pos: hit for hit in hits}
    assert by_pos[95].end == 104 and by_pos[95].kind == "deletion"
    assert not by_pos[95].in_3prime_5bp
    assert by_pos[116].end == 120 and by_pos[116].kind == "mnp"
    assert by_pos[116].in_3prime_5bp and by_pos[116].distance_from_3prime == 0


def test_strict_anchor_requires_exact_designed_coordinates():
    expected = expected_amplicon_from_design(pair(), template("+"))
    assert expected == {
        "subject": "chr1", "start": 5100, "end": 5468, "size": 369,
        "orientation": "F/R", "fwd_primer": "F", "rev_primer": "R",
    }
    exact = Amplicon("chr1", 5100, 5468, 369, "F", "R", 0, 0)
    local_duplicate = Amplicon("chr1", 5200, 5568, 369, "F", "R", 0, 0)
    result = reclassify_by_anchor(
        {"db": "db", "on_target": [exact, local_duplicate], "off_target": []},
        template=template("+"), pair=pair(), gel_min_gap=50,
    )
    assert result["intended_status"] == "unique"
    assert result["on_target"] == [exact]
    assert result["off_target"] == [local_duplicate]
    assert result["specific"] is False


def test_minus_template_requires_rf_orientation():
    expected = expected_amplicon_from_design(pair(), template("-"))
    assert expected["start"] == 6032 and expected["end"] == 6400
    assert expected["orientation"] == "R/F"
    exact = Amplicon("chr1", 6032, 6400, 369, "R", "F", 0, 0)
    result = reclassify_by_anchor(
        {"db": "db", "on_target": [exact], "off_target": []},
        template=template("-"), pair=pair(),
    )
    assert result["intended_status"] == "unique" and result["specific"] is True


def test_missing_and_ambiguous_intended_are_never_specific_or_low_risk():
    shifted = Amplicon("chr1", 5101, 5469, 369, "F", "R", 0, 0)
    missing = reclassify_by_anchor(
        {"db": "db", "on_target": [shifted], "off_target": []},
        template=template("+"), pair=pair(),
    )
    assert missing["intended_status"] == "missing" and not missing["specific"]
    assert assess_risk(intended_status="missing").level == "high"

    a = Amplicon("chr1", 5100, 5468, 369, "F", "R", 0, 0)
    b = Amplicon("chr1", 5100, 5468, 369, "F", "R", 0, 0)
    ambiguous = reclassify_by_anchor(
        {"db": "db", "on_target": [a, b], "off_target": []},
        template=template("+"), pair=pair(),
    )
    assert ambiguous["intended_status"] == "ambiguous"
    assert ambiguous["n_on_target"] == 0 and not ambiguous["specific"]
    assert assess_risk(intended_status="ambiguous").level == "high"
