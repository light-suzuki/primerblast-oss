"""Unit tests for the amplicon-pairing logic (no external tools needed)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from primerblast_oss.specificity import (  # noqa: E402
    PrimingSite, SpecParams, Amplicon, enumerate_amplicons,
    _conservative_intended_products, _count_3prime_mismatch, _hit_to_site,
    nearest_size_gap, priming_sites_with_stats, spec_params_for_profile,
)


def test_named_primer_pool_labels_products():
    sp = SpecParams()
    sites = [
        PrimingSite("PrimerA", "chr1", "+", 100, 0, 0, plen=20),
        PrimingSite("PrimerB", "chr1", "-", 600, 0, 0, plen=20),
    ]
    amps = enumerate_amplicons(sites, sp)
    assert len(amps) == 1
    assert amps[0].orientation == "PrimerA/PrimerB"


def test_nearest_size_gap():
    assert nearest_size_gap(500, [300, 540, 900]) == 40
    assert nearest_size_gap(500, []) is None


def _classify(amp, designed, tol=10):
    """Mirror the on-target rule in pair_specificity for testing."""
    perfect = amp.fwd_mismatch == 0 and amp.rev_mismatch == 0
    proper_pair = {amp.fwd_primer, amp.rev_primer} == {"F", "R"}
    size_ok = abs(amp.size - designed) <= tol
    return perfect and proper_pair and size_ok


def test_ortholog_on_opposite_strand_is_on_target():
    # R/F orientation, perfect, right size -> still the intended amplicon
    a = Amplicon("chrX", 100, 467, 368, "R", "F", 0, 0)
    assert _classify(a, designed=368) is True


def test_ff_product_is_never_on_target():
    a = Amplicon("chrX", 100, 467, 368, "F", "F", 0, 0)
    assert _classify(a, designed=368) is False


def test_perfect_wrong_size_is_off_target():
    a = Amplicon("chrX", 100, 900, 801, "F", "R", 0, 0)
    assert _classify(a, designed=368) is False


def test_duplicate_perfect_products_are_conservative_offtargets():
    a1 = Amplicon("chr1", 100, 467, 368, "F", "R", 0, 0, on_target=True)
    a2 = Amplicon("chr2", 1000, 1367, 368, "F", "R", 0, 0, on_target=True)
    on = _conservative_intended_products([a1, a2])
    assert on == [a1]
    assert a1.on_target is True
    assert a2.on_target is False
    assert a2.__dict__["ambiguous_intended_duplicate"] is True


def test_single_product():
    sp = SpecParams(min_product=40, max_product=4000)
    # F 3' at 100 (5' at 81), R 3' at 600 (5' at 619) -> product 81..619 = 539 bp
    sites = [
        PrimingSite("F", "chr1", "+", 100, 0, 0, plen=20),
        PrimingSite("R", "chr1", "-", 600, 0, 0, plen=20),
    ]
    amps = enumerate_amplicons(sites, sp)
    assert len(amps) == 1
    assert amps[0].size == 539
    assert amps[0].start == 81 and amps[0].end == 619
    assert amps[0].fwd_primer == "F" and amps[0].rev_primer == "R"


def test_offtarget_from_same_primer_both_strands():
    # a forward primer that also hits the minus strand downstream -> F/F product
    sp = SpecParams()
    sites = [
        PrimingSite("F", "chr2", "+", 1000, 0, 0, plen=20),   # 5' at 981
        PrimingSite("F", "chr2", "-", 1300, 0, 0, plen=20),   # 5' at 1319
        PrimingSite("R", "chr2", "-", 1500, 0, 0, plen=20),   # 5' at 1519
    ]
    amps = enumerate_amplicons(sites, sp)
    sizes = sorted(a.size for a in amps)
    assert sizes == [339, 539]  # F(981)->F(1319)=339, F(981)->R(1519)=539


def test_size_window_excludes_far_pairs():
    sp = SpecParams(min_product=40, max_product=1000)
    sites = [
        PrimingSite("F", "chr3", "+", 100, 0, 0, plen=20),
        PrimingSite("R", "chr3", "-", 5000, 0, 0, plen=20),  # too far
    ]
    assert enumerate_amplicons(sites, sp) == []


def test_different_subjects_do_not_pair():
    sp = SpecParams()
    sites = [
        PrimingSite("F", "chrA", "+", 100, 0, 0, plen=20),
        PrimingSite("R", "chrB", "-", 400, 0, 0, plen=20),
    ]
    assert enumerate_amplicons(sites, sp) == []


def test_3prime_mismatch_count():
    # 3'-terminal base matches; one mismatch two bases in
    mm, terminal = _count_3prime_mismatch("ACGTACGT", "ACGTAGGT", window=5)
    assert terminal is True
    assert mm == 1
    # 3'-terminal mismatch
    mm2, terminal2 = _count_3prime_mismatch("ACGTACGT", "ACGTACGA", window=5)
    assert terminal2 is False
    assert mm2 == 1


def test_ncbi_profile_keeps_terminal_mismatch_as_candidate_hit():
    q = "ACGTACGTACGTACGTACGT"
    s = "ACGTACGTACGTACGTACGA"  # one mismatch at the 3'-terminal base
    fields = ["primer", "chr1", "95.0", "20", "1", "0", "1", "20",
              "500", "519", "1e-5", "35", "plus", q, s, "20"]
    assert _hit_to_site(fields, "F", spec_params_for_profile("local-strict")) is None
    site = _hit_to_site(fields, "F", spec_params_for_profile("ncbi"))
    assert site is not None
    assert site.total_mismatch == 1
    assert site.tp_mismatch == 1


def test_ncbi_profile_keeps_five_total_mismatches_but_not_six():
    five_mm = ["primer", "chr1", "75.0", "20", "5", "0", "1", "20",
               "500", "519", "1e-5", "20", "plus",
               "AAAAACCCCCGGGGGTTTTT", "TTTTACCCCCGGGGGTTTTT", "20"]
    six_mm = list(five_mm)
    six_mm[4] = "6"
    six_mm[14] = "TTTTTCCCGGGGGTTTTT"
    sp = spec_params_for_profile("ncbi")
    assert _hit_to_site(five_mm, "F", sp) is not None
    assert _hit_to_site(six_mm, "F", sp) is None


def test_hit_rejected_when_3prime_not_aligned():
    sp = SpecParams()
    # qend (7) != qlen (20): 3' end of primer not in alignment -> no priming
    fields = ["primer", "chr1", "100.0", "17", "0", "0", "1", "7",
              "500", "516", "0.1", "30", "plus", "ACGTACG", "ACGTACG", "20"]
    assert _hit_to_site(fields, "F", sp) is None


def test_hit_accepted_full_length_perfect():
    sp = SpecParams()
    q = "ACGTACGTACGTACGTACGT"
    fields = ["primer", "chr1", "100.0", "20", "0", "0", "1", "20",
              "500", "519", "1e-5", "40", "plus", q, q, "20"]
    site = _hit_to_site(fields, "F", sp)
    assert site is not None
    assert site.strand == "+" and site.end3 == 519 and site.total_mismatch == 0


def test_priming_site_stats_warn_high_copy_hits(monkeypatch=None):
    class _Fake:
        pass

    def fake_run_blast(_primer, _db, _sp, _blastn):
        q = "ACGTACGTACGTACGTACGT"
        line = "\t".join([
            "primer", "chr1", "100.0", "20", "0", "0", "1", "20",
            "500", "519", "1e-5", "40", "plus", q, q, "20"
        ])
        return "\n".join([line, line])

    import primerblast_oss.specificity as S
    old = S._run_blast
    S._run_blast = fake_run_blast
    try:
        sites, stats = priming_sites_with_stats(
            "ACGTACGTACGTACGTACGT", "F", "db",
            SpecParams(max_target_seqs=10, high_copy_hit_threshold=2), "blastn")
    finally:
        S._run_blast = old
    assert len(sites) == 2
    assert stats.raw_blast_hits == 2
    assert stats.priming_sites == 2
    assert stats.unique_subjects == 1
    assert stats.near_target_limit is False
    assert stats.high_copy is True


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS {fn.__name__}")
        except Exception:
            failed += 1
            print(f"FAIL {fn.__name__}")
            traceback.print_exc()
    print(f"\n{len(fns)-failed}/{len(fns)} passed")
    raise SystemExit(1 if failed else 0)
