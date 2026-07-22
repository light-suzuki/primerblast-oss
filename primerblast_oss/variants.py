"""Variant-aware checks: primer footprints and amplicon conservation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence


@dataclass
class PrimerFootprint:
    primer: str
    chrom: str
    start: int
    end: int
    three_prime_coord: int
    strand: str


def footprints_from_amplicon(amp, len_f: int, len_r: int) -> List[PrimerFootprint]:
    """Return genomic footprints for named primers F and R.

    ``Amplicon.fwd_primer`` names the primer on the genomic plus strand (left
    side), while ``rev_primer`` names the primer on the genomic minus strand
    (right side). This may be F/R or R/F depending on target orientation.
    Results remain ordered as F then R for API compatibility.
    """
    lengths = {"F": len_f, "R": len_r}
    if {amp.fwd_primer, amp.rev_primer} != {"F", "R"}:
        raise ValueError("primer footprints require one F and one R site")

    left_len = lengths[amp.fwd_primer]
    right_len = lengths[amp.rev_primer]
    left = PrimerFootprint(
        primer=amp.fwd_primer,
        chrom=amp.subject,
        start=amp.start,
        end=amp.start + left_len - 1,
        three_prime_coord=amp.start + left_len - 1,
        strand="+",
    )
    right = PrimerFootprint(
        primer=amp.rev_primer,
        chrom=amp.subject,
        start=amp.end - right_len + 1,
        end=amp.end,
        three_prime_coord=amp.end - right_len + 1,
        strand="-",
    )
    by_name = {left.primer: left, right.primer: right}
    return [by_name["F"], by_name["R"]]


@dataclass
class SiteVariant:
    primer: str
    chrom: str
    pos: int
    end: int
    ref: str
    alt: List[str]
    kind: str
    in_3prime_5bp: bool
    distance_from_3prime: int


def _variant_end(variant) -> int:
    end = getattr(variant, "end", None)
    if end is not None:
        return int(end)
    return int(variant.pos) + max(1, len(str(variant.ref))) - 1


def variant_kind(variant) -> str:
    """Classify a VCF-like record conservatively from REF and ALT lengths."""
    ref = str(variant.ref)
    alts = [str(a) for a in variant.alt]
    if not alts or any(a.startswith("<") or "[" in a or "]" in a for a in alts):
        return "complex"
    ref_len = len(ref)
    alt_lens = [len(a) for a in alts]
    if ref_len == 1 and all(length == 1 for length in alt_lens):
        return "snp"
    if all(length == ref_len for length in alt_lens):
        return "mnp"
    if all(length > ref_len for length in alt_lens):
        return "insertion"
    if all(length < ref_len for length in alt_lens):
        return "deletion"
    return "complex"


def _distance_to_interval(coord: int, start: int, end: int) -> int:
    if start <= coord <= end:
        return 0
    return min(abs(coord - start), abs(coord - end))


def variants_under_primers(footprints: Sequence[PrimerFootprint],
                           variants: Sequence) -> List[SiteVariant]:
    """Return VCF records whose reference span overlaps a primer footprint.

    Insertions are treated conservatively at their VCF anchor base. For records
    spanning several reference bases, the closest affected base determines
    distance from the primer 3' end.
    """
    out: List[SiteVariant] = []
    for footprint in footprints:
        for variant in variants:
            if variant.chrom != footprint.chrom:
                continue
            start = int(variant.pos)
            end = _variant_end(variant)
            if end < footprint.start or start > footprint.end:
                continue
            distance = _distance_to_interval(footprint.three_prime_coord, start, end)
            out.append(SiteVariant(
                primer=footprint.primer,
                chrom=variant.chrom,
                pos=start,
                end=end,
                ref=variant.ref,
                alt=list(variant.alt),
                kind=variant_kind(variant),
                in_3prime_5bp=distance < 5,
                distance_from_3prime=distance,
            ))
    return out


def snps_under_primers(footprints: Sequence[PrimerFootprint],
                       variants: Sequence) -> List[SiteVariant]:
    """Backward-compatible alias for :func:`variants_under_primers`."""
    return variants_under_primers(footprints, variants)


def amplicon_variants(chrom: str, start: int, end: int,
                      variants: Sequence) -> List:
    """Variants whose reference span overlaps the amplicon."""
    lo, hi = min(start, end), max(start, end)
    return [variant for variant in variants
            if variant.chrom == chrom
            and _variant_end(variant) >= lo
            and int(variant.pos) <= hi]


def conservation_from_per_db(per_db: Sequence[Dict], designed_size: int,
                             tol: int = 15) -> Dict:
    conserved: List[str] = []
    details: Dict[str, Dict] = {}
    for result in per_db:
        db = result["db"].split("/")[-1]
        products = list(result.get("on_target", [])) + list(result.get("off_target", []))
        hit = None
        for amplicon in products:
            proper = {amplicon.fwd_primer, amplicon.rev_primer} == {"F", "R"}
            if proper and abs(amplicon.size - designed_size) <= tol:
                hit = amplicon
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
