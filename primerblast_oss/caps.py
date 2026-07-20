"""Restriction-enzyme, CAPS and dCAPS analysis.

The module is standard-library-only and keeps the historical ``ENZYMES`` mapping
for compatibility. New code uses structured :class:`RestrictionEnzyme` records
with strand-aware top/bottom cleavage offsets and PCR-substrate suitability.

Cut offsets are boundary coordinates measured from the first recognition base as
written 5'->3'. For example EcoRI ``G/AATTC`` is ``top_cut=1`` and
``bottom_cut=5``. Type IIS offsets may lie outside the recognition sequence;
BsaI ``GGTCTC(1/5)`` is represented as 7/11.

The curated recognition and cleavage notation follows the New England Biolabs
recognition-specificity chart and product documentation, curated 2026-07-20.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union


IUPAC_CODES: Dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}
_COMPLEMENT: Dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "D": "H",
    "H": "D", "V": "B", "N": "N",
}


def iupac_match(recognition: str, sequence_window: str) -> bool:
    recognition = recognition.upper()
    sequence_window = sequence_window.upper()
    if len(recognition) != len(sequence_window):
        return False
    return all(
        base in IUPAC_CODES.get(code, "")
        for code, base in zip(recognition, sequence_window)
    )


def revcomp(sequence: str) -> str:
    return "".join(
        _COMPLEMENT.get(base, "N") for base in reversed(sequence.upper()))


_METADATA_SOURCE = (
    "NEB recognition-specificity chart and product documentation; "
    "curated 2026-07-20"
)


@dataclass(frozen=True)
class RestrictionEnzyme:
    name: str
    recognition: str
    top_cut: int
    bottom_cut: int
    pcr_compatible: bool = True
    recommendable: bool = True
    methylation_note: str = ""
    cut_model: str = "verified"
    source: str = _METADATA_SOURCE

    @property
    def palindrome(self) -> bool:
        return revcomp(self.recognition) == self.recognition.upper()


@dataclass(frozen=True)
class Site:
    enzyme: str
    recognition: str
    pos: int
    strand: str


@dataclass(frozen=True)
class CutEvent:
    enzyme: str
    recognition: str
    site_pos: int
    site_strand: str
    top_cut: int
    bottom_cut: int
    complete: bool


@dataclass
class CapsResult:
    enzyme: str
    recognition: str
    allele_a_fragments: List[int]
    allele_b_fragments: List[int]
    distinguishable: bool
    min_gel_gap: int
    note: str = ""
    allele_a_cuts: List[CutEvent] = None  # type: ignore[assignment]
    allele_b_cuts: List[CutEvent] = None  # type: ignore[assignment]
    top_cut_offset: Optional[int] = None
    bottom_cut_offset: Optional[int] = None
    pcr_compatible: bool = True
    recommendation_eligible: bool = True
    methylation_note: str = ""
    cut_model: str = "verified"
    metadata_source: str = _METADATA_SOURCE

    def __post_init__(self) -> None:
        if self.allele_a_cuts is None:
            self.allele_a_cuts = []
        if self.allele_b_cuts is None:
            self.allele_b_cuts = []


def _enzyme(name: str, recognition: str, top: int, bottom: int,
            **kwargs) -> RestrictionEnzyme:
    return RestrictionEnzyme(name, recognition, top, bottom, **kwargs)


# Cleavage notation is encoded explicitly. DpnI is retained for discovery and
# provenance but is not recommended for ordinary unmethylated PCR products.
ENZYME_METADATA: Dict[str, RestrictionEnzyme] = {
    "EcoRI": _enzyme("EcoRI", "GAATTC", 1, 5),
    "HindIII": _enzyme("HindIII", "AAGCTT", 1, 5),
    "BamHI": _enzyme("BamHI", "GGATCC", 1, 5),
    "XbaI": _enzyme("XbaI", "TCTAGA", 1, 5),
    "XhoI": _enzyme("XhoI", "CTCGAG", 1, 5),
    "PstI": _enzyme("PstI", "CTGCAG", 5, 1),
    "SacI": _enzyme("SacI", "GAGCTC", 5, 1),
    "KpnI": _enzyme("KpnI", "GGTACC", 5, 1),
    "SmaI": _enzyme("SmaI", "CCCGGG", 3, 3),
    "NcoI": _enzyme("NcoI", "CCATGG", 1, 5),
    "NdeI": _enzyme("NdeI", "CATATG", 2, 4),
    "SalI": _enzyme("SalI", "GTCGAC", 1, 5),
    "NotI": _enzyme("NotI", "GCGGCCGC", 2, 6),
    "SpeI": _enzyme("SpeI", "ACTAGT", 1, 5),
    "EcoRV": _enzyme("EcoRV", "GATATC", 3, 3),
    "DraI": _enzyme("DraI", "TTTAAA", 3, 3),
    "HpaI": _enzyme("HpaI", "GTTAAC", 3, 3),
    "ScaI": _enzyme("ScaI", "AGTACT", 3, 3),
    "StuI": _enzyme("StuI", "AGGCCT", 3, 3),
    "AluI": _enzyme("AluI", "AGCT", 2, 2),
    "HaeIII": _enzyme("HaeIII", "GGCC", 2, 2),
    "RsaI": _enzyme("RsaI", "GTAC", 2, 2),
    "TaqI": _enzyme("TaqI", "TCGA", 1, 3),
    "MseI": _enzyme("MseI", "TTAA", 1, 3),
    "MboI": _enzyme(
        "MboI", "GATC", 0, 4,
        methylation_note=(
            "blocked by dam methylation; ordinary PCR products are generally "
            "unmethylated and therefore digestible")),
    "HinfI": _enzyme("HinfI", "GANTC", 1, 4),
    "DdeI": _enzyme("DdeI", "CTNAG", 1, 4),
    "BsaI": _enzyme("BsaI", "GGTCTC", 7, 11),
    "BstNI": _enzyme("BstNI", "CCWGG", 2, 3),
    "MspI": _enzyme("MspI", "CCGG", 1, 3),
    "HhaI": _enzyme("HhaI", "GCGC", 3, 1),
    "Sau3AI": _enzyme(
        "Sau3AI", "GATC", 0, 4,
        methylation_note="not blocked by dam methylation"),
    "NlaIII": _enzyme("NlaIII", "CATG", 4, 0),
    "DpnI": _enzyme(
        "DpnI", "GATC", 2, 2,
        pcr_compatible=False,
        recommendable=False,
        methylation_note=(
            "requires adenine-methylated GATC; ordinary PCR products are "
            "unmethylated and are not a substrate")),
    "TseI": _enzyme("TseI", "GCWGC", 1, 4),
    "ApaI": _enzyme("ApaI", "GGGCCC", 5, 1),
    "BglII": _enzyme("BglII", "AGATCT", 1, 5),
    "ClaI": _enzyme("ClaI", "ATCGAT", 2, 4),
    "NheI": _enzyme("NheI", "GCTAGC", 1, 5),
    "MfeI": _enzyme("MfeI", "CAATTG", 1, 5),
}

# Historical public API.
ENZYMES: Dict[str, str] = {
    name: enzyme.recognition for name, enzyme in ENZYME_METADATA.items()
}


EnzymeInput = Optional[Mapping[str, Union[str, RestrictionEnzyme]]]


def _coerce_enzyme(name: str,
                   value: Union[str, RestrictionEnzyme]) -> RestrictionEnzyme:
    if isinstance(value, RestrictionEnzyme):
        return value
    recognition = str(value).upper()
    known = ENZYME_METADATA.get(name)
    if known is not None and known.recognition == recognition:
        return known
    # Unknown custom enzymes remain searchable, but are excluded from automatic
    # experimental recommendations until cleavage metadata is supplied.
    return RestrictionEnzyme(
        name=name,
        recognition=recognition,
        top_cut=0,
        bottom_cut=0,
        pcr_compatible=False,
        recommendable=False,
        cut_model="recognition_start_fallback",
        source="user-supplied recognition sequence; cleavage metadata absent",
    )


def enzyme_records(enzymes: EnzymeInput = None,
                   recommended_only: bool = False) -> List[RestrictionEnzyme]:
    mapping: Mapping[str, Union[str, RestrictionEnzyme]] = (
        ENZYME_METADATA if enzymes is None else enzymes)
    records = [_coerce_enzyme(name, value) for name, value in mapping.items()]
    if recommended_only:
        records = [
            record for record in records
            if record.recommendable and record.pcr_compatible
            and record.cut_model == "verified"
        ]
    return records


def _is_palindrome(recognition: str) -> bool:
    return revcomp(recognition) == recognition.upper()


def find_sites(sequence: str, enzymes: EnzymeInput = None) -> List[Site]:
    sequence = sequence.upper()
    sites: List[Site] = []
    for enzyme in enzyme_records(enzymes):
        recognition = enzyme.recognition.upper()
        length = len(recognition)
        if length == 0 or length > len(sequence):
            continue
        reverse_recognition = revcomp(recognition)
        for start in range(len(sequence) - length + 1):
            window = sequence[start:start + length]
            if iupac_match(recognition, window):
                sites.append(Site(enzyme.name, recognition, start, "+"))
            if (not enzyme.palindrome
                    and iupac_match(reverse_recognition, window)):
                sites.append(Site(enzyme.name, recognition, start, "-"))
    sites.sort(key=lambda site: (site.pos, site.enzyme, site.strand))
    return sites


def cut_event_for_site(site: Site, enzyme: RestrictionEnzyme,
                       sequence_length: int) -> CutEvent:
    recognition_length = len(enzyme.recognition)
    if site.strand == "+":
        top_cut = site.pos + enzyme.top_cut
        bottom_cut = site.pos + enzyme.bottom_cut
    else:
        # The recognition strand is the genomic minus strand. Convert offsets
        # back to plus-strand boundary coordinates and swap strand roles.
        top_cut = site.pos + recognition_length - enzyme.bottom_cut
        bottom_cut = site.pos + recognition_length - enzyme.top_cut
    complete = (
        0 <= top_cut <= sequence_length
        and 0 <= bottom_cut <= sequence_length
    )
    return CutEvent(
        enzyme=enzyme.name,
        recognition=enzyme.recognition,
        site_pos=site.pos,
        site_strand=site.strand,
        top_cut=top_cut,
        bottom_cut=bottom_cut,
        complete=complete,
    )


def cut_events(sequence: str,
               enzyme: Union[str, RestrictionEnzyme],
               name: Optional[str] = None) -> List[CutEvent]:
    if isinstance(enzyme, RestrictionEnzyme):
        record = enzyme
    else:
        recognition = enzyme.upper()
        record = next((
            candidate for candidate in ENZYME_METADATA.values()
            if candidate.recognition == recognition and candidate.recommendable
        ), None)
        if record is None:
            record = RestrictionEnzyme(
                name or "custom", recognition, 0, 0,
                pcr_compatible=False, recommendable=False,
                cut_model="recognition_start_fallback",
                source="recognition-only compatibility fallback")
    sites = find_sites(sequence, {record.name: record})
    return [cut_event_for_site(site, record, len(sequence)) for site in sites]


def cut_positions(sequence: str,
                  recognition: Union[str, RestrictionEnzyme]) -> List[int]:
    """Return complete top-strand cleavage boundaries.

    Passing a known recognition string resolves to curated metadata. Unknown
    recognition-only inputs retain the historical recognition-start fallback.
    """
    return sorted(set(
        event.top_cut for event in cut_events(sequence, recognition)
        if event.complete
    ))


def digest_fragment_sizes(sequence: str,
                          recognition: Union[str, RestrictionEnzyme],
                          circular: bool = False) -> List[int]:
    sequence_length = len(sequence)
    cuts = sorted(set(
        cut for cut in cut_positions(sequence, recognition)
        if 0 < cut < sequence_length
    ))
    if not cuts:
        return [sequence_length]
    if not circular:
        boundaries = [0] + cuts + [sequence_length]
        return sorted([
            high - low for low, high in zip(boundaries, boundaries[1:])
            if high > low
        ], reverse=True)
    fragments = []
    for index, start in enumerate(cuts):
        end = cuts[(index + 1) % len(cuts)]
        length = (end - start) % sequence_length
        fragments.append(sequence_length if length == 0 else length)
    return sorted(fragments, reverse=True)


def _multiset_diff_min_gap(fragments_a: List[int],
                           fragments_b: List[int]) -> int:
    first = sorted(fragments_a, reverse=True)
    second = sorted(fragments_b, reverse=True)
    length = max(len(first), len(second))
    first += [0] * (length - len(first))
    second += [0] * (length - len(second))
    differences = [abs(a - b) for a, b in zip(first, second)]
    return max(differences) if differences else 0


def caps_scan(amplicon_a: str, amplicon_b: str,
              enzymes: EnzymeInput = None,
              gel_min_gap: int = 25,
              include_nonrecommended: bool = False) -> List[CapsResult]:
    records = enzyme_records(enzymes, recommended_only=not include_nonrecommended)
    results: List[CapsResult] = []
    for enzyme in records:
        events_a = cut_events(amplicon_a, enzyme)
        events_b = cut_events(amplicon_b, enzyme)
        fragments_a = digest_fragment_sizes(amplicon_a, enzyme)
        fragments_b = digest_fragment_sizes(amplicon_b, enzyme)
        if sorted(fragments_a) == sorted(fragments_b):
            continue
        gap = _multiset_diff_min_gap(fragments_a, fragments_b)
        incomplete = any(not event.complete for event in events_a + events_b)
        eligible = (
            enzyme.pcr_compatible and enzyme.recommendable
            and enzyme.cut_model == "verified" and not incomplete
        )
        distinguishable = eligible and gap >= gel_min_gap
        note_parts = []
        if len(fragments_a) != len(fragments_b):
            note_parts.append("cut-count differs")
        else:
            note_parts.append("same band count, sizes differ")
        if incomplete:
            note_parts.append("one or more cleavage positions lie outside amplicon")
        if not enzyme.pcr_compatible:
            note_parts.append("not compatible with ordinary PCR substrate")
        if enzyme.methylation_note:
            note_parts.append(enzyme.methylation_note)
        results.append(CapsResult(
            enzyme=enzyme.name,
            recognition=enzyme.recognition,
            allele_a_fragments=fragments_a,
            allele_b_fragments=fragments_b,
            distinguishable=distinguishable,
            min_gel_gap=gap,
            note="; ".join(note_parts),
            allele_a_cuts=events_a,
            allele_b_cuts=events_b,
            top_cut_offset=enzyme.top_cut,
            bottom_cut_offset=enzyme.bottom_cut,
            pcr_compatible=enzyme.pcr_compatible,
            recommendation_eligible=eligible,
            methylation_note=enzyme.methylation_note,
            cut_model=enzyme.cut_model,
            metadata_source=enzyme.source,
        ))
    results.sort(key=lambda result: (
        result.distinguishable,
        result.recommendation_eligible,
        result.min_gel_gap,
    ), reverse=True)
    return results


def enzymes_gained_lost(seq_ref: str, seq_alt: str,
                        enzymes: EnzymeInput = None,
                        include_nonrecommended: bool = False) -> Dict[str, List[str]]:
    records = enzyme_records(enzymes, recommended_only=not include_nonrecommended)
    gained: List[str] = []
    lost: List[str] = []
    for enzyme in records:
        reference_count = len(cut_events(seq_ref, enzyme))
        alternate_count = len(cut_events(seq_alt, enzyme))
        if alternate_count > reference_count:
            gained.append(enzyme.name)
        elif alternate_count < reference_count:
            lost.append(enzyme.name)
    return {"gained": sorted(gained), "lost": sorted(lost)}


def _pattern_orientations(enzyme: RestrictionEnzyme) -> List[Tuple[str, str]]:
    patterns = [("+", enzyme.recognition.upper())]
    reverse = revcomp(enzyme.recognition)
    if reverse != enzyme.recognition.upper():
        patterns.append(("-", reverse))
    return patterns


def dcaps_candidates(sequence: str, snp_index: int, ref_base: str,
                     alt_base: str, enzymes: EnzymeInput = None,
                     max_primer_mismatch: int = 2,
                     include_nonrecommended: bool = False) -> List[dict]:
    """Find allele-specific recognition frames on both strands.

    Engineered changes are returned as absolute plus-strand coordinates. The SNP
    itself is never counted as a primer mismatch and is not overwritten.
    """
    sequence = sequence.upper()
    ref_base = ref_base.upper()
    alt_base = alt_base.upper()
    if not (0 <= snp_index < len(sequence)):
        return []
    reference_sequence = (
        sequence[:snp_index] + ref_base + sequence[snp_index + 1:])
    alternate_sequence = (
        sequence[:snp_index] + alt_base + sequence[snp_index + 1:])
    records = enzyme_records(enzymes, recommended_only=not include_nonrecommended)
    candidates: List[dict] = []
    seen = set()

    for enzyme in records:
        recognition_length = len(enzyme.recognition)
        if recognition_length == 0 or recognition_length > len(sequence):
            continue
        minimum_start = max(0, snp_index - recognition_length + 1)
        maximum_start = min(len(sequence) - recognition_length, snp_index)
        for orientation, pattern in _pattern_orientations(enzyme):
            for start in range(minimum_start, maximum_start + 1):
                relative_snp = snp_index - start
                for present_sequence, present_in, other_sequence in (
                    (alternate_sequence, "alt", reference_sequence),
                    (reference_sequence, "ref", alternate_sequence),
                ):
                    window = present_sequence[start:start + recognition_length]
                    other_window = other_sequence[start:start + recognition_length]
                    if window[relative_snp] not in IUPAC_CODES.get(
                            pattern[relative_snp], ""):
                        continue
                    if other_window[relative_snp] in IUPAC_CODES.get(
                            pattern[relative_snp], ""):
                        continue
                    changes = []
                    feasible = True
                    for offset, (base, code) in enumerate(zip(window, pattern)):
                        if offset == relative_snp:
                            continue
                        allowed = IUPAC_CODES.get(code, "")
                        if not allowed:
                            feasible = False
                            break
                        if base not in allowed:
                            changes.append({
                                "position": start + offset,
                                "from": base,
                                "to": allowed[0],
                            })
                    if not feasible or len(changes) > max_primer_mismatch:
                        continue

                    engineered_ref = list(
                        reference_sequence[start:start + recognition_length])
                    engineered_alt = list(
                        alternate_sequence[start:start + recognition_length])
                    for change in changes:
                        offset = change["position"] - start
                        engineered_ref[offset] = change["to"]
                        engineered_alt[offset] = change["to"]
                    ref_context = "".join(engineered_ref)
                    alt_context = "".join(engineered_alt)
                    present_context = (
                        ref_context if present_in == "ref" else alt_context)
                    other_context = (
                        alt_context if present_in == "ref" else ref_context)
                    if not iupac_match(pattern, present_context):
                        continue
                    if iupac_match(pattern, other_context):
                        continue
                    key = (
                        enzyme.name, start, orientation, present_in,
                        tuple((change["position"], change["to"])
                              for change in changes),
                    )
                    if key in seen:
                        continue
                    seen.add(key)
                    candidates.append({
                        "enzyme": enzyme.name,
                        "recognition": enzyme.recognition,
                        "recognition_oriented": pattern,
                        "orientation": orientation,
                        "window_start": start,
                        "window_end": start + recognition_length - 1,
                        "snp_index": snp_index,
                        "snp_offset": relative_snp,
                        "ref_context": ref_context,
                        "alt_context": alt_context,
                        "mismatches": len(changes),
                        "engineered_changes": changes,
                        "present_in": present_in,
                        "top_cut_offset": enzyme.top_cut,
                        "bottom_cut_offset": enzyme.bottom_cut,
                        "pcr_compatible": enzyme.pcr_compatible,
                        "methylation_note": enzyme.methylation_note,
                        "metadata_source": enzyme.source,
                    })
    candidates.sort(key=lambda candidate: (
        candidate["mismatches"], candidate["enzyme"],
        candidate["window_start"], candidate["orientation"],
    ))
    return candidates


def materialize_dcaps_primers(sequence: str, snp_index: int,
                              candidates: Sequence[dict],
                              min_size: int = 18, opt_size: int = 20,
                              max_size: int = 25) -> List[dict]:
    """Turn recognition-frame candidates into orderable one-sided primers.

    The SNP is deliberately left outside the primer footprint so its allele is
    retained in the PCR product. A forward dCAPS primer ends immediately before
    the SNP; a reverse dCAPS primer ends immediately after it. Candidates whose
    engineered changes fall on both sides of the SNP require two modified
    primers and are rejected.
    """
    sequence = sequence.upper()
    output: List[dict] = []
    seen = set()
    for candidate in candidates:
        changes = list(candidate.get("engineered_changes", []))
        if not changes:
            continue  # natural CAPS, not dCAPS
        positions = [int(change["position"]) for change in changes]
        roles = []
        if all(position < snp_index for position in positions):
            roles.append("F")
        if all(position > snp_index for position in positions):
            roles.append("R")
        for role in roles:
            if role == "F":
                end_exclusive = snp_index
                required_start = min(positions)
                length = max(min_size, opt_size, end_exclusive - required_start)
                if length > max_size:
                    continue
                start = end_exclusive - length
                end = end_exclusive - 1
                if start < 0:
                    continue
                plus_bases = list(sequence[start:end_exclusive])
                for change in changes:
                    plus_bases[change["position"] - start] = change["to"]
                primer_sequence = "".join(plus_bases)
            else:
                start = snp_index + 1
                required_end = max(positions) + 1
                length = max(min_size, opt_size, required_end - start)
                if length > max_size:
                    continue
                end_exclusive = start + length
                end = end_exclusive - 1
                if end_exclusive > len(sequence):
                    continue
                plus_bases = list(sequence[start:end_exclusive])
                for change in changes:
                    plus_bases[change["position"] - start] = change["to"]
                primer_sequence = revcomp("".join(plus_bases))
            if not set(primer_sequence) <= set("ACGT"):
                continue
            key = (candidate["enzyme"], role, primer_sequence,
                   candidate["present_in"])
            if key in seen:
                continue
            seen.add(key)
            materialized = dict(candidate)
            materialized.update({
                "primer_role": role,
                "primer_sequence": primer_sequence,
                "primer_start": start,
                "primer_end": end,
                "primer_length": len(primer_sequence),
                "three_prime_distance_to_snp": 1,
                "orderable": True,
            })
            output.append(materialized)
    output.sort(key=lambda candidate: (
        candidate["mismatches"], candidate["primer_length"],
        candidate["enzyme"], candidate["primer_role"],
    ))
    return output


def apply_engineered_changes(sequence: str, changes: Sequence[dict]) -> str:
    bases = list(sequence.upper())
    for change in changes:
        position = int(change["position"])
        if 0 <= position < len(bases):
            bases[position] = str(change["to"]).upper()
    return "".join(bases)


def result_to_dict(result: CapsResult) -> dict:
    return asdict(result)


if __name__ == "__main__":
    demo = "GGGAATTCCC"
    eco_sites = [site for site in find_sites(demo) if site.enzyme == "EcoRI"]
    assert any(site.pos == 2 for site in eco_sites)
    eco_events = cut_events(demo, ENZYME_METADATA["EcoRI"])
    assert eco_events[0].top_cut == 3
    assert eco_events[0].bottom_cut == 7

    allele_a = "A" * 100 + "GAATTC" + "A" * 94
    allele_b = "A" * 100 + "GACTTC" + "A" * 94
    eco_result = next(
        result for result in caps_scan(allele_a, allele_b)
        if result.enzyme == "EcoRI")
    assert eco_result.distinguishable
    assert sum(eco_result.allele_a_fragments) == len(allele_a)
    assert sum(eco_result.allele_b_fragments) == len(allele_b)

    reverse_candidates = dcaps_candidates(
        "GAGACT", 5, "C", "T",
        enzymes={"BsaI": ENZYME_METADATA["BsaI"]},
        max_primer_mismatch=0,
    )
    assert any(
        candidate["orientation"] == "-"
        and candidate["present_in"] == "ref"
        for candidate in reverse_candidates)

    assert all(
        result.enzyme != "DpnI"
        for result in caps_scan("GATC", "AATC"))
    print("All self-tests passed.")
