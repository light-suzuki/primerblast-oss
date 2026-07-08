"""Restriction-enzyme / CAPS-dCAPS analysis for primerblast-oss.

Self-contained module (Python 3.8+, standard library only). No third-party
dependencies and no imports from other ``primerblast_oss`` modules.

Purpose
-------
Design CAPS (Cleaved Amplified Polymorphic Sequence) and dCAPS markers used in
plant breeding. Given a PCR amplicon that carries a SNP differing between two
alleles/parents, we look for restriction enzymes whose digestion pattern differs
between the two alleles, so the SNP can be scored on an agarose gel.

Coordinate / cut conventions (documented assumptions)
-----------------------------------------------------
* Sites are found by matching an enzyme's recognition sequence (IUPAC aware) on
  BOTH strands of a plain-ACGT template.
* We use the *start of the recognition match on the plus strand* as a proxy cut
  coordinate. Real enzymes cut at a defined offset inside/around the site, but
  for CAPS/dCAPS marker discovery the question is only "does the cut pattern
  differ between the two alleles", and the recognition-start proxy answers that
  faithfully as long as it is applied consistently to both alleles. Fragment
  sizes are therefore approximate but internally consistent.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

# ---------------------------------------------------------------------------
# IUPAC nucleotide handling
# ---------------------------------------------------------------------------

# Map each IUPAC ambiguity code to the set of plain bases it matches.
IUPAC_CODES: Dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",   # puRine
    "Y": "CT",   # pYrimidine
    "S": "GC",   # Strong
    "W": "AT",   # Weak
    "K": "GT",   # Keto
    "M": "AC",   # aMino
    "B": "CGT",  # not A
    "D": "AGT",  # not C
    "H": "ACT",  # not G
    "V": "ACG",  # not T
    "N": "ACGT",  # aNy
}

# IUPAC-aware complement table (covers ambiguity codes too).
_COMPLEMENT: Dict[str, str] = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",
}


def iupac_match(recognition: str, seq_window: str) -> bool:
    """Return True if ``seq_window`` (plain ACGT) matches ``recognition``.

    ``recognition`` may contain IUPAC ambiguity codes. Comparison is
    case-insensitive. Lengths must be equal.
    """
    recognition = recognition.upper()
    seq_window = seq_window.upper()
    if len(recognition) != len(seq_window):
        return False
    for r, s in zip(recognition, seq_window):
        allowed = IUPAC_CODES.get(r)
        if allowed is None:
            # Unknown symbol in the recognition site -> cannot match.
            return False
        if s not in allowed:
            return False
    return True


def revcomp(seq: str) -> str:
    """Reverse complement of ``seq`` (IUPAC aware, case-insensitive input).

    Output is upper-case. Unknown characters are passed through as 'N'.
    """
    seq = seq.upper()
    return "".join(_COMPLEMENT.get(base, "N") for base in reversed(seq))


# ---------------------------------------------------------------------------
# Enzyme table (~40 common commercially available enzymes)
# NAME -> recognition site (5'->3', plus strand). Recognition-based finding.
# ---------------------------------------------------------------------------

ENZYMES: Dict[str, str] = {
    "EcoRI": "GAATTC",
    "HindIII": "AAGCTT",
    "BamHI": "GGATCC",
    "XbaI": "TCTAGA",
    "XhoI": "CTCGAG",
    "PstI": "CTGCAG",
    "SacI": "GAGCTC",
    "KpnI": "GGTACC",
    "SmaI": "CCCGGG",
    "NcoI": "CCATGG",
    "NdeI": "CATATG",
    "SalI": "GTCGAC",
    "NotI": "GCGGCCGC",
    "SpeI": "ACTAGT",
    "EcoRV": "GATATC",
    "DraI": "TTTAAA",
    "HpaI": "GTTAAC",
    "ScaI": "AGTACT",
    "StuI": "AGGCCT",
    "AluI": "AGCT",
    "HaeIII": "GGCC",
    "RsaI": "GTAC",
    "TaqI": "TCGA",
    "MseI": "TTAA",
    "MboI": "GATC",
    "HinfI": "GANTC",
    "DdeI": "CTNAG",
    "BsaI": "GGTCTC",
    "BstNI": "CCWGG",
    "MspI": "CCGG",
    "HhaI": "GCGC",
    "Sau3AI": "GATC",
    "NlaIII": "CATG",
    "DpnI": "GATC",
    "TseI": "GCWGC",
    "ApaI": "GGGCCC",
    "BglII": "AGATCT",
    "ClaI": "ATCGAT",
    "NheI": "GCTAGC",
    "MfeI": "CAATTG",
}


# ---------------------------------------------------------------------------
# Site finding
# ---------------------------------------------------------------------------

@dataclass
class Site:
    """A single restriction recognition match on a template.

    Attributes
    ----------
    enzyme : str
        Enzyme name.
    recognition : str
        Recognition sequence (may contain IUPAC codes).
    pos : int
        0-based start of the recognition match, expressed in plus-strand
        coordinates of the input sequence.
    strand : str
        '+' if the recognition site was found on the plus strand, '-' if the
        recognition site matched on the minus (reverse-complement) strand.
    """

    enzyme: str
    recognition: str
    pos: int
    strand: str


def _is_palindrome(recognition: str) -> bool:
    """True if a recognition site equals its own IUPAC reverse complement.

    Palindromic sites (e.g. GAATTC) match at the same locus on both strands and
    would otherwise be double-counted.
    """
    return revcomp(recognition) == recognition.upper()


def find_sites(seq: str, enzymes: Optional[Dict[str, str]] = None) -> List[Site]:
    """Scan ``seq`` on BOTH strands for every enzyme's recognition site.

    Returns a list of :class:`Site` objects. Palindromic recognition sites are
    reported once (plus-strand) to avoid double-counting the same locus.
    Results are sorted by (pos, enzyme, strand) for determinism.
    """
    if enzymes is None:
        enzymes = ENZYMES
    seq = seq.upper()
    n = len(seq)
    sites: List[Site] = []

    for name, recognition in enzymes.items():
        rec = recognition.upper()
        rlen = len(rec)
        if rlen == 0 or rlen > n:
            continue
        palindrome = _is_palindrome(rec)
        rec_rc = revcomp(rec)

        for i in range(n - rlen + 1):
            window = seq[i:i + rlen]
            # Plus-strand match.
            if iupac_match(rec, window):
                sites.append(Site(name, recognition, i, "+"))
            # Minus-strand match: the recognition site appears on the minus
            # strand when the plus-strand window matches the reverse complement
            # of the recognition sequence. Skip for palindromes (same locus).
            if not palindrome and iupac_match(rec_rc, window):
                sites.append(Site(name, recognition, i, "-"))

    sites.sort(key=lambda s: (s.pos, s.enzyme, s.strand))
    return sites


def cut_positions(seq: str, recognition: str) -> List[int]:
    """Sorted, de-duplicated cut coordinates for one recognition site.

    We use the recognition-site start (plus-strand coordinate) as the proxy cut
    point (see module docstring). Both strands are scanned; palindromic sites
    are counted once per locus.
    """
    seq = seq.upper()
    rec = recognition.upper()
    rlen = len(rec)
    n = len(seq)
    if rlen == 0 or rlen > n:
        return []
    palindrome = _is_palindrome(rec)
    rec_rc = revcomp(rec)

    positions = set()
    for i in range(n - rlen + 1):
        window = seq[i:i + rlen]
        if iupac_match(rec, window):
            positions.add(i)
        elif not palindrome and iupac_match(rec_rc, window):
            positions.add(i)
    return sorted(positions)


def digest_fragment_sizes(seq: str, recognition: str, circular: bool = False) -> List[int]:
    """Fragment lengths from cutting ``seq`` at all recognition sites.

    Returns sizes sorted descending. If there are no sites, returns
    ``[len(seq)]`` (linear) or ``[len(seq)]`` (circular, single uncut circle).

    Cutting model: we cut at each proxy cut coordinate (recognition start).
    * Linear: cuts partition the sequence into successive intervals.
    * Circular: the region after the last cut wraps to the first cut.
    """
    n = len(seq)
    cuts = cut_positions(seq, recognition)
    if not cuts:
        return [n]

    if not circular:
        fragments: List[int] = []
        prev = 0
        for c in cuts:
            fragments.append(c - prev)
            prev = c
        fragments.append(n - prev)
        # A cut at position 0 yields a leading zero-length fragment; drop
        # non-positive fragments as they are not physical bands.
        fragments = [f for f in fragments if f > 0]
        return sorted(fragments, reverse=True)

    # Circular: fragment between consecutive cuts, last wraps around.
    fragments = []
    for idx in range(len(cuts)):
        start = cuts[idx]
        end = cuts[(idx + 1) % len(cuts)]
        length = (end - start) % n
        if length == 0:
            length = n
        fragments.append(length)
    return sorted(fragments, reverse=True)


# ---------------------------------------------------------------------------
# CAPS scanning
# ---------------------------------------------------------------------------

@dataclass
class CapsResult:
    """Result of comparing one enzyme's digest between two alleles."""

    enzyme: str
    recognition: str
    allele_a_fragments: List[int]
    allele_b_fragments: List[int]
    distinguishable: bool
    min_gel_gap: int
    note: str = ""


def _multiset_diff_min_gap(frags_a: List[int], frags_b: List[int]) -> int:
    """Smallest fragment-size difference that lets you tell A from B.

    Heuristic (documented): sort both fragment lists descending and pad the
    shorter with zeros so they align by band rank. The per-rank absolute size
    differences are the candidate "gaps"; we return the *largest* such gap,
    because on a gel you only need ONE clearly resolvable band-size difference
    to distinguish the alleles. When the two patterns are identical the value is
    0. Returns 0 for identical patterns.
    """
    a = sorted(frags_a, reverse=True)
    b = sorted(frags_b, reverse=True)
    length = max(len(a), len(b))
    a += [0] * (length - len(a))
    b += [0] * (length - len(b))
    diffs = [abs(x - y) for x, y in zip(a, b)]
    return max(diffs) if diffs else 0


def caps_scan(
    amplicon_a: str,
    amplicon_b: str,
    enzymes: Optional[Dict[str, str]] = None,
    gel_min_gap: int = 25,
) -> List[CapsResult]:
    """Compare digests of two alleles for every enzyme.

    For each enzyme, digest both alleles and report the enzyme when the
    fragment-size multiset differs. ``distinguishable`` is True when there is a
    band-size difference >= ``gel_min_gap`` between the two patterns (i.e. a
    difference resolvable on a typical agarose gel). Results are sorted best
    (largest gap) first.
    """
    if enzymes is None:
        enzymes = ENZYMES

    results: List[CapsResult] = []
    for name, recognition in enzymes.items():
        frags_a = digest_fragment_sizes(amplicon_a, recognition)
        frags_b = digest_fragment_sizes(amplicon_b, recognition)

        # Only interesting when the fragment patterns actually differ.
        if sorted(frags_a) == sorted(frags_b):
            continue

        gap = _multiset_diff_min_gap(frags_a, frags_b)
        distinguishable = gap >= gel_min_gap

        n_a = len(frags_a)
        n_b = len(frags_b)
        if n_a != n_b:
            note = "cut-count differs (%d vs %d sites+1)" % (n_a, n_b)
        else:
            note = "same band count, sizes differ"

        results.append(
            CapsResult(
                enzyme=name,
                recognition=recognition,
                allele_a_fragments=frags_a,
                allele_b_fragments=frags_b,
                distinguishable=distinguishable,
                min_gel_gap=gap,
                note=note,
            )
        )

    # Best first: distinguishable enzymes first, then by largest gel gap.
    results.sort(key=lambda r: (r.distinguishable, r.min_gel_gap), reverse=True)
    return results


def enzymes_gained_lost(
    seq_ref: str,
    seq_alt: str,
    enzymes: Optional[Dict[str, str]] = None,
) -> Dict[str, List[str]]:
    """Enzymes whose site count changes between reference and alternate allele.

    Returns ``{'gained': [...], 'lost': [...]}`` where:
      * 'gained' = enzymes with MORE sites in ``seq_alt`` than ``seq_ref``
        (site created by the variant).
      * 'lost'   = enzymes with FEWER sites in ``seq_alt`` than ``seq_ref``
        (site destroyed by the variant).

    This is the core of CAPS marker discovery: an enzyme that gains or loses a
    site between the two alleles cleaves them differently.
    """
    if enzymes is None:
        enzymes = ENZYMES

    gained: List[str] = []
    lost: List[str] = []
    for name, recognition in enzymes.items():
        n_ref = len(cut_positions(seq_ref, recognition))
        n_alt = len(cut_positions(seq_alt, recognition))
        if n_alt > n_ref:
            gained.append(name)
        elif n_alt < n_ref:
            lost.append(name)
    return {"gained": sorted(gained), "lost": sorted(lost)}


# ---------------------------------------------------------------------------
# dCAPS (best-effort)
# ---------------------------------------------------------------------------

def dcaps_candidates(
    seq: str,
    snp_index: int,
    ref_base: str,
    alt_base: str,
    enzymes: Optional[Dict[str, str]] = None,
    max_primer_mismatch: int = 2,
) -> List[dict]:
    """Best-effort dCAPS candidate finder.

    dCAPS (derived CAPS) is used when a SNP does not by itself create or destroy
    a restriction site. A primer is designed with a few deliberate mismatches
    near the SNP so that, combined with one of the two SNP alleles, an enzyme
    site is created for exactly one allele.

    Implementation / assumptions
    ----------------------------
    We search, within a window spanning any recognition site that overlaps
    ``snp_index``, for enzymes where introducing up to ``max_primer_mismatch``
    base changes in the bases *around* the SNP (excluding the SNP position
    itself) can produce a recognition site that is present for one SNP allele
    but not the other. Only the enzyme's recognition length is considered
    (fixed-length sites); IUPAC codes in the recognition are honoured. The SNP
    base itself is the discriminating position and is not counted as a primer
    mismatch.

    Returns a list of dict candidates:
        {
          'enzyme': name,
          'recognition': site,
          'window_start': int,       # 0-based start of the recognition frame
          'ref_context': str,        # engineered window for the ref allele
          'alt_context': str,        # engineered window for the alt allele
          'mismatches': int,         # engineered primer mismatches used
          'present_in': 'ref'|'alt', # which allele carries the created site
        }
    """
    if enzymes is None:
        enzymes = ENZYMES

    seq = seq.upper()
    ref_base = ref_base.upper()
    alt_base = alt_base.upper()
    n = len(seq)
    if not (0 <= snp_index < n):
        return []

    # Two allele templates that differ only at the SNP position.
    ref_seq = seq[:snp_index] + ref_base + seq[snp_index + 1:]
    alt_seq = seq[:snp_index] + alt_base + seq[snp_index + 1:]

    candidates: List[dict] = []
    seen = set()

    for name, recognition in enzymes.items():
        rec = recognition.upper()
        rlen = len(rec)
        if rlen == 0 or rlen > n:
            continue

        # Consider every recognition frame that overlaps the SNP position.
        frame_start_min = max(0, snp_index - rlen + 1)
        frame_start_max = min(n - rlen, snp_index)
        for start in range(frame_start_min, frame_start_max + 1):
            rel = snp_index - start  # SNP offset within the frame
            if not (0 <= rel < rlen):
                continue

            for allele_seq, present_in, other_seq in (
                (alt_seq, "alt", ref_seq),
                (ref_seq, "ref", alt_seq),
            ):
                window = allele_seq[start:start + rlen]
                other_window = other_seq[start:start + rlen]

                # Count how many positions (excluding the SNP) must be changed
                # so that this window matches the recognition site. The SNP base
                # is free (it is the discriminator, not a primer mismatch).
                mism = 0
                feasible = True
                for j in range(rlen):
                    if j == rel:
                        # SNP position must be compatible with the site for the
                        # allele that is supposed to CARRY the created site.
                        if window[j] not in IUPAC_CODES.get(rec[j], ""):
                            feasible = False
                            break
                        continue
                    if window[j] not in IUPAC_CODES.get(rec[j], ""):
                        mism += 1
                if not feasible or mism > max_primer_mismatch:
                    continue

                # The engineered site must be ALLELE-SPECIFIC: with the same
                # engineered primer bases, the OTHER allele must NOT form the
                # site (its SNP base breaks it).
                other_forms = other_window[rel] in IUPAC_CODES.get(rec[rel], "")
                if other_forms:
                    continue

                # Build the engineered contexts (primer bases forced to match).
                eng_ref = list(ref_seq[start:start + rlen])
                eng_alt = list(alt_seq[start:start + rlen])
                for j in range(rlen):
                    if j == rel:
                        continue
                    if window[j] not in IUPAC_CODES.get(rec[j], ""):
                        # Force a concrete base that satisfies the site.
                        forced = IUPAC_CODES[rec[j]][0]
                        eng_ref[j] = forced
                        eng_alt[j] = forced
                ref_context = "".join(eng_ref)
                alt_context = "".join(eng_alt)

                key = (name, start, present_in, mism)
                if key in seen:
                    continue
                seen.add(key)

                candidates.append(
                    {
                        "enzyme": name,
                        "recognition": recognition,
                        "window_start": start,
                        "ref_context": ref_context,
                        "alt_context": alt_context,
                        "mismatches": mism,
                        "present_in": present_in,
                    }
                )

    # Fewest engineered mismatches first (easier primer to build).
    candidates.sort(key=lambda c: (c["mismatches"], c["enzyme"], c["window_start"]))
    return candidates


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # (a) EcoRI found at the right position in GGGAATTCCC.
    #     GAATTC starts at index 2 (G G G A A T T C C C -> GAATTC at 2..7).
    demo = "GGGAATTCCC"
    all_sites = find_sites(demo)
    ecori_sites = [s for s in all_sites if s.enzyme == "EcoRI"]
    print("=== (a) find_sites('GGGAATTCCC') EcoRI ===")
    for s in ecori_sites:
        print("  EcoRI at pos %d strand %s" % (s.pos, s.strand))
    assert any(s.pos == 2 for s in ecori_sites), "EcoRI should be found at pos 2"
    print("  OK: EcoRI found at position 2")
    print()

    # (b) Two 200 bp alleles differing by one SNP that creates/destroys EcoRI.
    #     Build a random-ish but deterministic 200 bp backbone with NO EcoRI
    #     site, then place a GAATTC in allele A and break it in allele B.
    import random

    rng = random.Random(42)
    bases = "ACGT"
    backbone = "".join(rng.choice(bases) for _ in range(200))

    # Ensure the backbone has no accidental EcoRI site by rebuilding until clean.
    while "GAATTC" in backbone or "GAATTC" in revcomp(backbone):
        backbone = "".join(rng.choice(bases) for _ in range(200))

    # Insert an EcoRI site at position 100 in allele A.
    site_pos = 100
    allele_a = backbone[:site_pos] + "GAATTC" + backbone[site_pos + 6:]
    # Allele B: single SNP breaks the EcoRI site (GAATTC -> GACTTC), i.e. the
    # SNP is at site_pos+2 (A->C).
    allele_b = backbone[:site_pos] + "GACTTC" + backbone[site_pos + 6:]

    assert len(allele_a) == 200 and len(allele_b) == 200
    assert allele_a[:site_pos] == allele_b[:site_pos]
    # Exactly one base differs.
    diffs = [i for i in range(200) if allele_a[i] != allele_b[i]]
    assert diffs == [site_pos + 2], "alleles must differ by exactly one SNP"

    print("=== (b) caps_scan on 200 bp alleles (EcoRI SNP) ===")
    results = caps_scan(allele_a, allele_b, gel_min_gap=25)
    ecori_result = next((r for r in results if r.enzyme == "EcoRI"), None)
    assert ecori_result is not None, "EcoRI must appear as a CAPS-distinguishing enzyme"
    print("  EcoRI recognition       : %s" % ecori_result.recognition)
    print("  Allele A fragment sizes : %s" % ecori_result.allele_a_fragments)
    print("  Allele B fragment sizes : %s" % ecori_result.allele_b_fragments)
    print("  Distinguishable         : %s" % ecori_result.distinguishable)
    print("  Min gel gap (bp)        : %s" % ecori_result.min_gel_gap)
    print("  Note                    : %s" % ecori_result.note)
    assert ecori_result.distinguishable, "EcoRI should distinguish the two alleles"
    print("  Top CAPS enzymes (best first):")
    for r in results[:5]:
        print("    %-8s dist=%-5s gap=%-4d A=%s B=%s"
              % (r.enzyme, r.distinguishable, r.min_gel_gap,
                 r.allele_a_fragments, r.allele_b_fragments))
    print()

    # (c) enzymes_gained_lost for the same pair (ref = A, alt = B).
    print("=== (c) enzymes_gained_lost(ref=alleleA, alt=alleleB) ===")
    gl = enzymes_gained_lost(allele_a, allele_b)
    print("  gained (site present only in alt): %s" % gl["gained"])
    print("  lost   (site present only in ref): %s" % gl["lost"])
    assert "EcoRI" in gl["lost"], "EcoRI site is lost going A -> B"
    print("  OK: EcoRI correctly reported as lost")
    print()

    print("All self-tests passed.")
