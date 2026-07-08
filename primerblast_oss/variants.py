"""Variant-aware checks: SNPs under primers and amplicon conservation.

Two breeding-critical questions NCBI Primer-BLAST cannot answer locally:
  * does a known SNP/indel (from a VCF) sit under a primer -- especially its
    3' end, where it wrecks amplification or, deliberately, enables an
    allele-specific assay?
  * is the amplicon conserved across several reference / parent genomes, so the
    same primers work on all of them?
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence


@dataclass
class PrimerFootprint:
    primer: str             # "F" or "R"
    chrom: str
    start: int              # 1-based genomic, 5' <= 3' span on the top strand
    end: int
    three_prime_coord: int  # genomic coordinate of the primer's 3' base


def footprints_from_amplicon(amp, len_f: int, len_r: int) -> List[PrimerFootprint]:
    """Genomic footprints of the two primers of an on-target amplicon.

    `amp` is an Amplicon: amp.start = 5' of the plus-strand primer, amp.end =
    5' of the minus-strand primer (both 1-based genomic)."""
    chrom = amp.subject
    fwd = PrimerFootprint("F", chrom, amp.start, amp.start + len_f - 1,
                          three_prime_coord=amp.start + len_f - 1)
    rev = PrimerFootprint("R", chrom, amp.end - len_r + 1, amp.end,
                          three_prime_coord=amp.end - len_r + 1)
    return [fwd, rev]


@dataclass
class SiteVariant:
    primer: str
    chrom: str
    pos: int
    ref: str
    alt: List[str]
    in_3prime_5bp: bool     # variant within the 3'-terminal 5 bp of the primer
    distance_from_3prime: int


def snps_under_primers(footprints: Sequence[PrimerFootprint],
                       variants: Sequence) -> List[SiteVariant]:
    """Variants overlapping any primer footprint, flagged if they fall in the
    3'-terminal 5 bp (the amplification-critical zone)."""
    out: List[SiteVariant] = []
    for fp in footprints:
        for v in variants:
            if v.chrom != fp.chrom:
                continue
            if v.pos < fp.start or v.pos > fp.end:
                continue
            dist = abs(v.pos - fp.three_prime_coord)
            out.append(SiteVariant(
                primer=fp.primer, chrom=v.chrom, pos=v.pos, ref=v.ref,
                alt=list(v.alt), in_3prime_5bp=dist < 5, distance_from_3prime=dist))
    return out


def amplicon_variants(chrom: str, start: int, end: int,
                      variants: Sequence) -> List:
    """Variants inside the amplicon span (the polymorphisms a CAPS/dCAPS marker
    would genotype)."""
    lo, hi = min(start, end), max(start, end)
    return [v for v in variants if v.chrom == chrom and lo <= v.pos <= hi]


def conservation_from_per_db(per_db: Sequence[Dict], designed_size: int,
                             tol: int = 15) -> Dict:
    """Which reference databases conserve the amplicon.

    A reference conserves it if it yields a proper forward+reverse product of
    ~the designed size (mismatches allowed -- the primers still bind, which is
    what "conserved" means for a working assay)."""
    conserved: List[str] = []
    details: Dict[str, Dict] = {}
    for d in per_db:
        db = d["db"].split("/")[-1]
        products = list(d.get("on_target", [])) + list(d.get("off_target", []))
        hit = None
        for a in products:
            proper = {a.fwd_primer, a.rev_primer} == {"F", "R"}
            if proper and abs(a.size - designed_size) <= tol:
                hit = a
                break
        details[db] = {
            "conserved": hit is not None,
            "size": hit.size if hit else None,
            "fwd_mismatch": hit.fwd_mismatch if hit else None,
            "rev_mismatch": hit.rev_mismatch if hit else None,
        }
        if hit is not None:
            conserved.append(db)
    return {
        "conserved_in": conserved,
        "n_conserved": len(conserved),
        "n_refs": len(per_db),
        "fully_conserved": len(conserved) == len(per_db) and len(per_db) > 0,
        "per_ref": details,
    }
