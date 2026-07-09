"""Pure-logic tests for the breeding pipeline modules (no external tools)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from primerblast_oss.genome import revcomp                                # noqa: E402
from primerblast_oss.design import PrimerPair                             # noqa: E402
from primerblast_oss.regions import GenomicRegion, tile_interval          # noqa: E402
from primerblast_oss.specificity import Amplicon                          # noqa: E402
from primerblast_oss.report import tiling_to_dict                         # noqa: E402
from primerblast_oss.variants import (                                    # noqa: E402
    footprints_from_amplicon, snps_under_primers, amplicon_variants,
    conservation_from_per_db,
)
from primerblast_oss.risk import assess_risk                              # noqa: E402
from primerblast_oss.assay import _orientation_kind, reclassify_by_anchor  # noqa: E402
from primerblast_oss.caps import caps_scan, find_sites                    # noqa: E402


class _V:  # minimal Variant stand-in
    def __init__(self, chrom, pos, ref, alt):
        self.chrom, self.pos, self.ref, self.alt = chrom, pos, ref, alt


def test_revcomp():
    assert revcomp("AAAACGT") == "ACGTTTT"
    assert revcomp("ACGTN") == "NACGT"


def test_tile_interval_n_markers():
    r = GenomicRegion("chr1", 1000, 11000, "+", "QTL")
    pts = tile_interval(r, n_markers=4)
    assert len(pts) == 4
    assert all(r.start < p.start < r.end for p in pts)
    # evenly spaced
    gaps = [pts[i+1].start - pts[i].start for i in range(len(pts)-1)]
    assert len(set(gaps)) == 1


def test_footprints_and_3prime_snp():
    # intended amplicon: F 5' at 100 (len 20 -> 3' at 119), R 5' at 500 (3' at 481)
    amp = Amplicon("chr1", 100, 500, 401, "F", "R", 0, 0, on_target=True)
    fps = footprints_from_amplicon(amp, 20, 20)
    assert fps[0].three_prime_coord == 119 and fps[1].three_prime_coord == 481
    # SNP at 118 -> within 3'-5bp of the forward primer
    hits = snps_under_primers(fps, [_V("chr1", 118, "A", ["G"])])
    assert len(hits) == 1 and hits[0].in_3prime_5bp is True
    # SNP at 105 -> under primer but not near 3'
    hits2 = snps_under_primers(fps, [_V("chr1", 105, "A", ["G"])])
    assert hits2[0].in_3prime_5bp is False


def test_amplicon_variants():
    vs = [_V("chr1", 300, "A", ["T"]), _V("chr1", 9000, "C", ["G"]), _V("chr2", 300, "A", ["T"])]
    got = amplicon_variants("chr1", 100, 500, vs)
    assert len(got) == 1 and got[0].pos == 300


def test_conservation_across_refs():
    a1 = Amplicon("chr1", 100, 468, 369, "F", "R", 0, 0, on_target=True)
    a2 = Amplicon("chrX", 100, 468, 369, "F", "R", 1, 0)   # 1 mismatch but binds -> conserved
    per_db = [
        {"db": "/x/cameor", "on_target": [a1], "off_target": []},
        {"db": "/x/ZW6", "on_target": [], "off_target": [a2]},
        {"db": "/x/JI2694", "on_target": [], "off_target": []},   # not conserved
    ]
    c = conservation_from_per_db(per_db, designed_size=369)
    assert c["n_conserved"] == 2
    assert set(c["conserved_in"]) == {"cameor", "ZW6"}
    assert c["fully_conserved"] is False


def test_orientation_kind():
    assert _orientation_kind(Amplicon("c", 1, 2, 2, "F", "R", 0, 0)) == "FR"
    assert _orientation_kind(Amplicon("c", 1, 2, 2, "F", "F", 0, 0)) == "FF"
    assert _orientation_kind(Amplicon("c", 1, 2, 2, "R", "R", 0, 0)) == "RR"


def test_risk_levels():
    assert assess_risk().level == "low"
    # off-target that primes (perfect 3') -> not low
    assert assess_risk(n_ff=1, offtarget_min_tp5=0).level in ("medium", "high")
    # co-migrating off-target -> high
    assert assess_risk(n_comigrating_offtarget=1).level == "high"
    # SNP in 3' end -> never "low" (may fail, or be deliberately allele-specific)
    assert assess_risk(snp_in_primer_3prime=True).level in ("medium", "high")


def test_reclassify_by_anchor_paralogs():
    # three perfect same-size products; only the one at the template locus
    # (chr1:5000-6500) is intended, the paralogs on chr1 far away and chr2 are off
    intended = Amplicon("chr1", 5100, 5468, 369, "F", "R", 0, 0)
    para1 = Amplicon("chr1", 900000, 900368, 369, "F", "R", 0, 0)
    para2 = Amplicon("chr2", 100, 468, 369, "F", "R", 0, 0)
    res = {"db": "/x/cameor", "on_target": [intended, para1, para2], "off_target": []}
    out = reclassify_by_anchor(res, "chr1", 5000, 6500, 369)
    assert out["n_on_target"] == 1
    assert out["on_target"][0] is intended
    assert out["n_off_target"] == 2          # paralogs surfaced
    assert out["n_comigrating"] == 2          # same size -> co-migrating
    assert out["gel_distinguishable"] is False


def test_thermo_viability_optional():
    from primerblast_oss import thermo
    if not thermo.available():
        return  # primer3-py not installed -> feature is optional, skip
    rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    P = "GCACTCTAGAGGTTCAAGGCC"
    perfect = thermo.evaluate(P, rc(P))
    assert perfect is not None and perfect.viable and perfect.tm > 55
    # an unrelated target does not anneal -> very low Tm -> not viable
    weak = thermo.evaluate(P, "TTATTATGATCGATAGCTAGC")
    assert weak is not None and not weak.viable


def test_caps_ecori_distinguishable():
    left = "ACGT" * 30
    right = "TGCA" * 30
    a = left + "GAATTC" + right      # has EcoRI
    b = left + "GACTTC" + right      # SNP destroys EcoRI
    res = caps_scan(a, b)
    enz = [r for r in res if r.enzyme == "EcoRI"]
    assert enz and enz[0].distinguishable
    assert find_sites("GGGAATTCCC")  # sanity


def test_tiling_json_keeps_amplicons_structured():
    pair = PrimerPair(
        index=0, template_id="t", forward="ACGT", reverse="TGCA",
        left_start=0, left_len=4, right_start=99, right_len=4,
        product_size=100, tm_f=60.0, tm_r=60.0, gc_f=50.0, gc_r=50.0,
    )
    amp = Amplicon("chr1", 1, 100, 100, "F", "R", 0, 0, on_target=True)
    pair.specificity = {
        "per_db": [{
            "db": "/x/db", "on_target": [amp], "off_target": [],
            "n_products": 1, "n_on_target": 1, "n_off_target": 0,
            "n_comigrating": 0, "specific": True,
        }]
    }
    out = tiling_to_dict([{"index": 1, "covers": (0, 99), "pair": pair}],
                         "t", (0, 99), ["/x/db"])
    product = out["tiles"][0]["specificity"]["per_db"][0]["on_target"][0]
    assert product["subject"] == "chr1"
    assert product["orientation"] == "F/R"


if __name__ == "__main__":
    import traceback
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn(); print(f"PASS {fn.__name__}")
        except Exception:
            failed += 1; print(f"FAIL {fn.__name__}"); traceback.print_exc()
    print(f"\n{len(fns)-failed}/{len(fns)} passed")
    raise SystemExit(1 if failed else 0)
