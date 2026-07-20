"""Specificity and in-silico PCR.

Primer BLAST HSPs are converted to 3'-anchored priming sites and paired into
convergent PCR products. The module also records whether BLAST returned enough
evidence to support a definitive specificity statement. A clean *observed* hit
list is not called exhaustive when target enumeration may have been capped.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

_OUTFMT = (
    "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send "
    "evalue bitscore sstrand qseq sseq qlen"
)

SEARCH_COMPLETE = "complete"
SEARCH_POSSIBLY_TRUNCATED = "possibly_truncated"
SEARCH_REPEAT_LIMITED = "truncated_or_repeat_limited"
_SEARCH_SEVERITY = {
    SEARCH_COMPLETE: 0,
    SEARCH_POSSIBLY_TRUNCATED: 1,
    SEARCH_REPEAT_LIMITED: 2,
}


def combine_search_completeness(states: Sequence[str]) -> str:
    """Return the most conservative completeness state in ``states``."""
    known = [state for state in states if state in _SEARCH_SEVERITY]
    if not known:
        return SEARCH_COMPLETE
    return max(known, key=lambda state: _SEARCH_SEVERITY[state])


@dataclass
class SpecParams:
    max_total_mismatch: int = 4
    max_3prime_mismatch: int = 1
    three_prime_window: int = 5
    require_3prime_terminal_match: bool = True
    min_product: int = 40
    max_product: int = 4000
    gel_min_gap_bp: int = 50
    word_size: int = 7
    evalue: float = 30000.0
    max_target_seqs: int = 5000
    high_copy_hit_threshold: int = 10000
    high_copy_site_threshold: int = 500
    dust: str = "no"
    reward: int = 1
    penalty: int = -1
    gapopen: int = 5
    gapextend: int = 2
    num_threads: int = 4


SPECIFICITY_PROFILES = {
    "local-strict": {
        "max_total_mismatch": 4,
        "max_3prime_mismatch": 1,
        "three_prime_window": 5,
        "require_3prime_terminal_match": True,
    },
    "ncbi": {
        "max_total_mismatch": 5,
        "max_3prime_mismatch": 1,
        "three_prime_window": 5,
        "require_3prime_terminal_match": False,
    },
}


def spec_params_for_profile(profile: str = "local-strict", **overrides) -> SpecParams:
    if profile not in SPECIFICITY_PROFILES:
        names = ", ".join(sorted(SPECIFICITY_PROFILES))
        raise ValueError("unknown specificity profile '%s'; choose one of: %s" % (profile, names))
    params = SpecParams()
    for key, value in SPECIFICITY_PROFILES[profile].items():
        setattr(params, key, value)
    for key, value in overrides.items():
        if value is not None:
            setattr(params, key, value)
    return params


@dataclass
class PrimingSite:
    primer: str
    subject: str
    strand: str
    end3: int
    total_mismatch: int
    tp_mismatch: int
    plen: int = 20
    tp5_mismatch: int = 0
    tp10_mismatch: int = 0
    tm: Optional[float] = None
    end3_dg: Optional[float] = None
    thermo_viable: Optional[bool] = None

    @property
    def end5(self) -> int:
        return self.end3 - (self.plen - 1) if self.strand == "+" else self.end3 + (self.plen - 1)

    @property
    def extends(self) -> str:
        return "right" if self.strand == "+" else "left"

    def describe(self) -> str:
        return ("%s binds %s %s strand, 5'=%s 3'=%s, extends %s, "
                "mm total=%s 3'5bp=%s" % (
                    self.primer, self.subject, self.strand, self.end5, self.end3,
                    self.extends, self.total_mismatch, self.tp5_mismatch))


@dataclass
class PrimerHitStats:
    primer: str
    raw_blast_hits: int
    priming_sites: int
    unique_subjects: int
    near_target_limit: bool
    high_copy: bool
    at_target_limit: bool = False
    completeness: str = SEARCH_COMPLETE


@dataclass
class Amplicon:
    subject: str
    start: int
    end: int
    size: int
    fwd_primer: str
    rev_primer: str
    fwd_mismatch: int
    rev_mismatch: int
    on_target: bool = False
    fwd_tp5: int = 0
    rev_tp5: int = 0
    fwd_tm: Optional[float] = None
    rev_tm: Optional[float] = None
    fwd_end3_dg: Optional[float] = None
    rev_end3_dg: Optional[float] = None

    @property
    def orientation(self) -> str:
        return "%s/%s" % (self.fwd_primer, self.rev_primer)


def _detect_blastn(explicit: Optional[str]) -> str:
    for candidate in (explicit, "blastn"):
        if candidate and shutil.which(candidate):
            return shutil.which(candidate)  # type: ignore[return-value]
    raise RuntimeError("blastn not found. Install BLAST+ or pass blastn_bin.")


def _run_blast(primer: str, db: str, sp: SpecParams, blastn: str) -> str:
    query = (">primer\n%s\n" % primer).encode()
    cmd = [
        blastn, "-task", "blastn-short", "-db", db, "-query", "-",
        "-outfmt", _OUTFMT,
        "-word_size", str(sp.word_size),
        "-evalue", str(sp.evalue),
        "-max_target_seqs", str(sp.max_target_seqs),
        "-dust", sp.dust,
        "-reward", str(sp.reward), "-penalty", str(sp.penalty),
        "-gapopen", str(sp.gapopen), "-gapextend", str(sp.gapextend),
        "-soft_masking", "false",
        "-num_threads", str(sp.num_threads),
    ]
    proc = subprocess.run(cmd, input=query, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError("blastn failed: %s" % proc.stderr.decode(errors="ignore"))
    return proc.stdout.decode(errors="ignore")


def _count_3prime_mismatch(qseq: str, sseq: str, window: int) -> Tuple[int, bool]:
    """Count mismatches in the final query bases; every gapped column counts."""
    mismatches = 0
    query_bases_seen = 0
    terminal_match = False
    first_query_base = True
    for query_base, subject_base in zip(reversed(qseq), reversed(sseq)):
        if query_base == "-":
            mismatches += 1
            continue
        if first_query_base:
            terminal_match = query_base == subject_base
            first_query_base = False
        if query_base != subject_base:
            mismatches += 1
        query_bases_seen += 1
        if query_bases_seen >= window:
            break
    return mismatches, terminal_match


def _alignment_edit_count(qseq: str, sseq: str) -> int:
    """Count substitutions and every gapped alignment column."""
    return sum(q != s for q, s in zip(qseq, sseq)) + abs(len(qseq) - len(sseq))


def _hit_to_site(fields: List[str], primer_id: str, sp: SpecParams) -> Optional[PrimingSite]:
    (_query_id, subject_id, _identity, _length, _mismatch, _gapopen,
     query_start, query_end, _subject_start, subject_end, _evalue, _bits,
     subject_strand, query_sequence, subject_sequence, query_length) = fields
    query_end_i = int(query_end)
    query_length_i = int(query_length)
    if query_end_i != query_length_i:
        return None

    tp_mismatch, terminal_match = _count_3prime_mismatch(
        query_sequence, subject_sequence, sp.three_prime_window)
    if sp.require_3prime_terminal_match and not terminal_match:
        return None

    total_mismatch = _alignment_edit_count(query_sequence, subject_sequence) + int(query_start) - 1
    if total_mismatch > sp.max_total_mismatch or tp_mismatch > sp.max_3prime_mismatch:
        return None

    strand = "+" if subject_strand == "plus" else "-"
    tp5, _ = _count_3prime_mismatch(query_sequence, subject_sequence, 5)
    tp10, _ = _count_3prime_mismatch(query_sequence, subject_sequence, 10)
    return PrimingSite(
        primer=primer_id,
        subject=subject_id,
        strand=strand,
        end3=int(subject_end),
        total_mismatch=total_mismatch,
        tp_mismatch=tp_mismatch,
        plen=query_length_i,
        tp5_mismatch=tp5,
        tp10_mismatch=tp10,
    )


def _classify_hit_list(raw_hits: int, priming_sites: int, unique_subjects: int,
                       sp: SpecParams) -> Tuple[bool, bool, str]:
    """Return ``(near_limit, high_copy, completeness)`` conservatively.

    BLAST does not expose a reliable truncation bit for ``max_target_seqs``.
    Reaching the subject cap, or producing a repeat-scale hit/site list, is
    therefore treated as repeat-limited. The final five percent below the cap is
    considered possibly truncated.
    """
    cap = max(1, int(sp.max_target_seqs))
    near_limit = unique_subjects >= max(1, int(cap * 0.95))
    at_limit = unique_subjects >= cap
    high_copy = (
        raw_hits >= sp.high_copy_hit_threshold
        or priming_sites >= sp.high_copy_site_threshold
    )
    if at_limit or high_copy:
        completeness = SEARCH_REPEAT_LIMITED
    elif near_limit:
        completeness = SEARCH_POSSIBLY_TRUNCATED
    else:
        completeness = SEARCH_COMPLETE
    return near_limit, high_copy, completeness


def priming_sites_with_stats(
    primer: str, primer_id: str, db: str, sp: SpecParams, blastn: str
) -> Tuple[List[PrimingSite], PrimerHitStats]:
    output = _run_blast(primer, db, sp, blastn)
    sites: List[PrimingSite] = []
    raw_hits = 0
    subjects = set()
    for line in output.splitlines():
        if not line.strip():
            continue
        raw_hits += 1
        fields = line.split("\t")
        if len(fields) < 16:
            continue
        subjects.add(fields[1])
        site = _hit_to_site(fields, primer_id, sp)
        if site is not None:
            sites.append(site)

    near_limit, high_copy, completeness = _classify_hit_list(
        raw_hits, len(sites), len(subjects), sp)
    stats = PrimerHitStats(
        primer=primer_id,
        raw_blast_hits=raw_hits,
        priming_sites=len(sites),
        unique_subjects=len(subjects),
        near_target_limit=near_limit,
        high_copy=high_copy,
        at_target_limit=len(subjects) >= max(1, int(sp.max_target_seqs)),
        completeness=completeness,
    )
    return sites, stats


def priming_sites(primer: str, primer_id: str, db: str, sp: SpecParams,
                  blastn: str) -> List[PrimingSite]:
    sites, _stats = priming_sites_with_stats(primer, primer_id, db, sp, blastn)
    return sites


def screen_primers(primers: Dict[str, str], db: str, sp: SpecParams,
                   blastn: str) -> List[PrimingSite]:
    sites: List[PrimingSite] = []
    for name, sequence in primers.items():
        sites.extend(priming_sites(sequence, name, db, sp, blastn))
    return sites


def screen_primers_with_stats(
    primers: Dict[str, str], db: str, sp: SpecParams, blastn: str
) -> Tuple[List[PrimingSite], Dict[str, PrimerHitStats]]:
    sites: List[PrimingSite] = []
    stats: Dict[str, PrimerHitStats] = {}
    for name, sequence in primers.items():
        primer_sites, primer_stats = priming_sites_with_stats(
            sequence, name, db, sp, blastn)
        sites.extend(primer_sites)
        stats[name] = primer_stats
    return sites, stats


def _search_metadata(hit_stats: Dict[str, PrimerHitStats], sp: SpecParams) -> Dict:
    per_primer = {name: stats.completeness for name, stats in hit_stats.items()}
    overall = combine_search_completeness(list(per_primer.values()))
    return {
        "search_completeness": overall,
        "search_complete": overall == SEARCH_COMPLETE,
        "primer_search_completeness": per_primer,
        "raw_hits_per_primer": {name: stats.raw_blast_hits for name, stats in hit_stats.items()},
        "unique_subjects_per_primer": {
            name: stats.unique_subjects for name, stats in hit_stats.items()
        },
        "near_blast_limit": [
            name for name, stats in hit_stats.items() if stats.near_target_limit
        ],
        "at_blast_limit": [
            name for name, stats in hit_stats.items() if stats.at_target_limit
        ],
        "high_copy_primers": [
            name for name, stats in hit_stats.items() if stats.high_copy
        ],
        "blast_limits": {
            "max_target_seqs": sp.max_target_seqs,
            "evalue": sp.evalue,
            "word_size": sp.word_size,
            "num_threads": sp.num_threads,
            "high_copy_hit_threshold": sp.high_copy_hit_threshold,
            "high_copy_site_threshold": sp.high_copy_site_threshold,
        },
        "completeness_recommendation": (
            None if overall == SEARCH_COMPLETE else
            "Rerun with --exhaustive or a larger --max-target-seqs; if the primer "
            "remains repeat-limited, redesign it or use an indexed alternative search."
        ),
    }


def _site_binding_strand(genome, site: PrimingSite) -> str:
    low, high = min(site.end5, site.end3), max(site.end5, site.end3)
    strand = "-" if site.strand == "+" else "+"
    try:
        return genome.fetch(site.subject, low, high, strand)
    except Exception:  # noqa: BLE001
        return ""


def annotate_thermo(sites: Sequence[PrimingSite], primers: Dict[str, str],
                    genome, tp=None, gate: bool = True):
    """Annotate sites and report exactly what could be evaluated.

    ``tp is False`` is an explicit disable sentinel used by ``--no-thermo``.
    Missing contigs, out-of-range coordinates, or failed calculations are
    counted as unresolved instead of being silently described as evaluated.
    """
    from . import thermo as thermo_module
    stats = {
        "attempted_per_primer": {},
        "evaluated_per_primer": {},
        "unresolved_per_primer": {},
        "gated_per_primer": {},
    }
    if tp is False or not thermo_module.available() or genome is None:
        return list(sites), {}, stats

    viable: Dict[str, int] = {}
    kept: List[PrimingSite] = []
    for site in sites:
        stats["attempted_per_primer"][site.primer] = (
            stats["attempted_per_primer"].get(site.primer, 0) + 1)
        sequence = primers.get(site.primer)
        target = _site_binding_strand(genome, site) if sequence else ""
        result = (
            thermo_module.evaluate(sequence, target, tp)
            if (sequence and target) else None
        )
        if result is None:
            stats["unresolved_per_primer"][site.primer] = (
                stats["unresolved_per_primer"].get(site.primer, 0) + 1)
            kept.append(site)
            continue

        stats["evaluated_per_primer"][site.primer] = (
            stats["evaluated_per_primer"].get(site.primer, 0) + 1)
        site.tm = result.tm
        site.end3_dg = result.end3_dg
        site.thermo_viable = result.viable
        if result.viable:
            viable[site.primer] = viable.get(site.primer, 0) + 1
        elif gate:
            stats["gated_per_primer"][site.primer] = (
                stats["gated_per_primer"].get(site.primer, 0) + 1)
            continue
        kept.append(site)
    return kept, viable, stats


def enumerate_amplicons(sites: Sequence[PrimingSite], sp: SpecParams) -> List[Amplicon]:
    by_subject: Dict[str, List[PrimingSite]] = {}
    for site in sites:
        by_subject.setdefault(site.subject, []).append(site)

    amplicons: List[Amplicon] = []
    for subject, group in by_subject.items():
        plus_sites = sorted(
            (site for site in group if site.strand == "+"), key=lambda site: site.end5)
        minus_sites = sorted(
            (site for site in group if site.strand == "-"), key=lambda site: site.end5)
        for forward_site in plus_sites:
            for reverse_site in minus_sites:
                if reverse_site.end3 < forward_site.end3:
                    continue
                start = forward_site.end5
                end = reverse_site.end5
                size = end - start + 1
                if size < sp.min_product:
                    continue
                if size > sp.max_product:
                    break
                amplicons.append(Amplicon(
                    subject=subject,
                    start=start,
                    end=end,
                    size=size,
                    fwd_primer=forward_site.primer,
                    rev_primer=reverse_site.primer,
                    fwd_mismatch=forward_site.total_mismatch,
                    rev_mismatch=reverse_site.total_mismatch,
                    fwd_tp5=forward_site.tp5_mismatch,
                    rev_tp5=reverse_site.tp5_mismatch,
                    fwd_tm=forward_site.tm,
                    rev_tm=reverse_site.tm,
                    fwd_end3_dg=forward_site.end3_dg,
                    rev_end3_dg=reverse_site.end3_dg,
                ))
    amplicons.sort(key=lambda amplicon: (amplicon.subject, amplicon.start))
    return amplicons


def nearest_size_gap(size: int, others: Sequence[int]) -> Optional[int]:
    gaps = [abs(size - other) for other in others]
    return min(gaps) if gaps else None


def _conservative_intended_products(amplicons: Sequence[Amplicon]) -> List[Amplicon]:
    candidates = [amplicon for amplicon in amplicons if amplicon.on_target]
    if len(candidates) <= 1:
        return candidates
    for candidate in candidates[1:]:
        candidate.on_target = False
        candidate.__dict__["ambiguous_intended_duplicate"] = True
    candidates[0].__dict__["generic_intended_candidate"] = True
    return candidates[:1]


def _specificity_verdict(observed_specific: bool, completeness: str):
    if not observed_specific:
        return False
    if completeness != SEARCH_COMPLETE:
        return None
    return True


def in_silico_pcr(
    primers: Dict[str, str],
    db: str,
    sp: Optional[SpecParams] = None,
    blastn_bin: Optional[str] = None,
    genome=None,
    thermo_params=None,
    thermo_gate: bool = True,
) -> Dict:
    sp = sp or SpecParams()
    blastn = _detect_blastn(blastn_bin)
    sites, hit_stats = screen_primers_with_stats(primers, db, sp, blastn)
    sites, viable_sites, thermo_site_stats = annotate_thermo(
        sites, primers, genome, thermo_params, thermo_gate)
    amplicons = enumerate_amplicons(sites, sp)

    sizes = [amplicon.size for amplicon in amplicons]
    for index, amplicon in enumerate(amplicons):
        amplicon.__dict__["nearest_gap"] = nearest_size_gap(
            amplicon.size,
            [size for other_index, size in enumerate(sizes) if other_index != index],
        )

    metadata = _search_metadata(hit_stats, sp)
    return {
        "db": db,
        "primers": dict(primers),
        "sites_per_primer": {
            name: sum(site.primer == name for site in sites) for name in primers
        },
        "thermo_evaluated": bool(sum(
            thermo_site_stats["evaluated_per_primer"].values())),
        "thermo_site_stats": thermo_site_stats,
        "viable_sites_per_primer": viable_sites,
        **metadata,
        "n_products": len(amplicons),
        "products": amplicons,
    }


def pair_specificity(
    forward: str,
    reverse: str,
    db: str,
    designed_size: Optional[int] = None,
    sp: Optional[SpecParams] = None,
    blastn_bin: Optional[str] = None,
    size_tolerance: int = 10,
    genome=None,
    thermo_params=None,
    thermo_gate: bool = True,
    allowed_primer_mismatches: Optional[Mapping[str, int]] = None,
) -> Dict:
    sp = sp or SpecParams()
    blastn = _detect_blastn(blastn_bin)
    primers = {"F": forward, "R": reverse}
    allowed_mismatches = {
        name: max(0, int(value))
        for name, value in (allowed_primer_mismatches or {}).items()
    }
    sites, hit_stats = screen_primers_with_stats(primers, db, sp, blastn)
    sites, viable_sites, thermo_site_stats = annotate_thermo(
        sites, primers, genome, thermo_params, thermo_gate)
    amplicons = enumerate_amplicons(sites, sp)

    for amplicon in amplicons:
        mismatch_ok = (
            amplicon.fwd_mismatch
            <= allowed_mismatches.get(amplicon.fwd_primer, 0)
            and amplicon.rev_mismatch
            <= allowed_mismatches.get(amplicon.rev_primer, 0)
        )
        proper_pair = {amplicon.fwd_primer, amplicon.rev_primer} == {"F", "R"}
        size_ok = designed_size is None or abs(amplicon.size - designed_size) <= size_tolerance
        amplicon.on_target = mismatch_ok and proper_pair and size_ok
    _conservative_intended_products(amplicons)

    on_target = [amplicon for amplicon in amplicons if amplicon.on_target]
    off_target = [amplicon for amplicon in amplicons if not amplicon.on_target]
    reference_size = designed_size if designed_size is not None else (
        on_target[0].size if on_target else None)
    comigrating = (
        [amplicon for amplicon in off_target
         if abs(amplicon.size - reference_size) < sp.gel_min_gap_bp]
        if reference_size is not None else []
    )
    nearest_off_gap = (
        nearest_size_gap(reference_size, [amplicon.size for amplicon in off_target])
        if reference_size is not None and off_target else None
    )

    metadata = _search_metadata(hit_stats, sp)
    observed_specific = len(amplicons) == 1 and len(on_target) == 1
    specific = _specificity_verdict(observed_specific, metadata["search_completeness"])
    if specific is True:
        specificity_status = "specific"
    elif specific is None:
        specificity_status = "indeterminate"
    else:
        specificity_status = "non_specific"

    return {
        "db": db,
        "n_forward_sites": sum(site.primer == "F" for site in sites),
        "n_reverse_sites": sum(site.primer == "R" for site in sites),
        "thermo_evaluated": bool(sum(
            thermo_site_stats["evaluated_per_primer"].values())),
        "thermo_site_stats": thermo_site_stats,
        "viable_sites_per_primer": viable_sites,
        **metadata,
        "n_products": len(amplicons),
        "n_on_target": len(on_target),
        "n_off_target": len(off_target),
        "n_comigrating": len(comigrating),
        "nearest_offtarget_gap": nearest_off_gap,
        "gel_distinguishable": len(comigrating) == 0,
        "on_target": on_target,
        "off_target": off_target,
        "specific_observed": observed_specific,
        "specific": specific,
        "specificity_status": specificity_status,
        "allowed_primer_mismatches": allowed_mismatches,
    }
