"""Specificity / in-silico PCR: the core Primer-BLAST-like step.

Each primer is searched with BLAST, accepted HSPs become 3'-anchored priming
sites, and convergent sites are paired into predicted PCR products. Pairing is
orientation-agnostic, so F/R, R/F, F/F and R/R products are all represented.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq qlen"


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
        raise ValueError(f"unknown specificity profile '{profile}'; choose one of: {names}")
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
        return (f"{self.primer} binds {self.subject} {self.strand} strand, "
                f"5'={self.end5} 3'={self.end3}, extends {self.extends}, "
                f"mm total={self.total_mismatch} 3'5bp={self.tp5_mismatch}")


@dataclass
class PrimerHitStats:
    primer: str
    raw_blast_hits: int
    priming_sites: int
    unique_subjects: int
    near_target_limit: bool
    high_copy: bool


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
        return f"{self.fwd_primer}/{self.rev_primer}"


def _detect_blastn(explicit: Optional[str]) -> str:
    for cand in (explicit, "blastn"):
        if cand and shutil.which(cand):
            return shutil.which(cand)  # type: ignore[return-value]
    raise RuntimeError("blastn not found. Install BLAST+ or pass blastn_bin.")


def _run_blast(primer: str, db: str, sp: SpecParams, blastn: str) -> str:
    query = f">primer\n{primer}\n".encode()
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
        raise RuntimeError(f"blastn failed: {proc.stderr.decode(errors='ignore')}")
    return proc.stdout.decode(errors="ignore")


def _count_3prime_mismatch(qseq: str, sseq: str, window: int) -> Tuple[int, bool]:
    """Count mismatches in the final query bases; every gapped column counts."""
    mm = 0
    seen = 0
    terminal_match = False
    first = True
    for q, s in zip(reversed(qseq), reversed(sseq)):
        if q == "-":
            mm += 1
            continue
        if first:
            terminal_match = q == s
            first = False
        if q != s:
            mm += 1
        seen += 1
        if seen >= window:
            break
    return mm, terminal_match


def _alignment_edit_count(qseq: str, sseq: str) -> int:
    """Count substitutions and every gapped alignment column.

    BLAST's ``gapopen`` field counts gap runs, not their length. Counting unequal
    aligned columns prevents a multi-base indel from being treated as one edit.
    """
    return sum(q != s for q, s in zip(qseq, sseq)) + abs(len(qseq) - len(sseq))


def _hit_to_site(fields: List[str], primer_id: str, sp: SpecParams) -> Optional[PrimingSite]:
    (_qid, sseqid, _pident, _length, _mismatch, _gapopen, qstart, qend,
     _sstart, send, _e, _bits, sstrand, qseq, sseq, qlen) = fields
    qend_i = int(qend)
    qlen_i = int(qlen)
    if qend_i != qlen_i:
        return None

    tp_mm, terminal_match = _count_3prime_mismatch(qseq, sseq, sp.three_prime_window)
    if sp.require_3prime_terminal_match and not terminal_match:
        return None

    # Unequal aligned columns include substitutions and each gap base. Unaligned
    # 5' query bases are additional edits against the full primer.
    total_mm = _alignment_edit_count(qseq, sseq) + int(qstart) - 1
    if total_mm > sp.max_total_mismatch or tp_mm > sp.max_3prime_mismatch:
        return None

    strand = "+" if sstrand == "plus" else "-"
    tp5, _ = _count_3prime_mismatch(qseq, sseq, 5)
    tp10, _ = _count_3prime_mismatch(qseq, sseq, 10)
    return PrimingSite(
        primer=primer_id, subject=sseqid, strand=strand, end3=int(send),
        total_mismatch=total_mm, tp_mismatch=tp_mm, plen=qlen_i,
        tp5_mismatch=tp5, tp10_mismatch=tp10,
    )


def priming_sites_with_stats(
    primer: str, primer_id: str, db: str, sp: SpecParams, blastn: str
) -> Tuple[List[PrimingSite], PrimerHitStats]:
    out = _run_blast(primer, db, sp, blastn)
    sites: List[PrimingSite] = []
    raw_hits = 0
    subjects = set()
    for line in out.splitlines():
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
    stats = PrimerHitStats(
        primer=primer_id,
        raw_blast_hits=raw_hits,
        priming_sites=len(sites),
        unique_subjects=len(subjects),
        near_target_limit=len(subjects) >= max(1, int(sp.max_target_seqs * 0.95)),
        high_copy=len(sites) >= sp.high_copy_site_threshold,
    )
    return sites, stats


def priming_sites(primer: str, primer_id: str, db: str, sp: SpecParams, blastn: str) -> List[PrimingSite]:
    sites, _stats = priming_sites_with_stats(primer, primer_id, db, sp, blastn)
    return sites


def screen_primers(primers: Dict[str, str], db: str, sp: SpecParams, blastn: str) -> List[PrimingSite]:
    sites: List[PrimingSite] = []
    for name, seq in primers.items():
        sites.extend(priming_sites(seq, name, db, sp, blastn))
    return sites


def screen_primers_with_stats(
    primers: Dict[str, str], db: str, sp: SpecParams, blastn: str
) -> Tuple[List[PrimingSite], Dict[str, PrimerHitStats]]:
    sites: List[PrimingSite] = []
    stats: Dict[str, PrimerHitStats] = {}
    for name, seq in primers.items():
        primer_sites, primer_stats = priming_sites_with_stats(seq, name, db, sp, blastn)
        sites.extend(primer_sites)
        stats[name] = primer_stats
    return sites, stats


def _site_binding_strand(genome, site: PrimingSite) -> str:
    lo, hi = min(site.end5, site.end3), max(site.end5, site.end3)
    strand = "-" if site.strand == "+" else "+"
    try:
        return genome.fetch(site.subject, lo, hi, strand)
    except Exception:  # noqa: BLE001
        return ""


def annotate_thermo(sites: Sequence[PrimingSite], primers: Dict[str, str],
                    genome, tp=None, gate: bool = True):
    from . import thermo as _thermo
    if not _thermo.available() or genome is None:
        return list(sites), {}
    viable: Dict[str, int] = {}
    kept: List[PrimingSite] = []
    for site in sites:
        seq = primers.get(site.primer)
        target = _site_binding_strand(genome, site) if seq else ""
        res = _thermo.evaluate(seq, target, tp) if (seq and target) else None
        if res is not None:
            site.tm, site.end3_dg, site.thermo_viable = res.tm, res.end3_dg, res.viable
            if res.viable:
                viable[site.primer] = viable.get(site.primer, 0) + 1
            elif gate:
                continue
        kept.append(site)
    return kept, viable


def enumerate_amplicons(sites: Sequence[PrimingSite], sp: SpecParams) -> List[Amplicon]:
    by_subject: Dict[str, List[PrimingSite]] = {}
    for site in sites:
        by_subject.setdefault(site.subject, []).append(site)

    amplicons: List[Amplicon] = []
    for subject, group in by_subject.items():
        plus = sorted((s for s in group if s.strand == "+"), key=lambda x: x.end5)
        minus = sorted((s for s in group if s.strand == "-"), key=lambda x: x.end5)
        for fwd in plus:
            for rev in minus:
                if rev.end3 < fwd.end3:
                    continue
                start = fwd.end5
                end = rev.end5
                size = end - start + 1
                if size < sp.min_product:
                    continue
                if size > sp.max_product:
                    break
                amplicons.append(Amplicon(
                    subject=subject, start=start, end=end, size=size,
                    fwd_primer=fwd.primer, rev_primer=rev.primer,
                    fwd_mismatch=fwd.total_mismatch, rev_mismatch=rev.total_mismatch,
                    fwd_tp5=fwd.tp5_mismatch, rev_tp5=rev.tp5_mismatch,
                    fwd_tm=fwd.tm, rev_tm=rev.tm,
                    fwd_end3_dg=fwd.end3_dg, rev_end3_dg=rev.end3_dg,
                ))
    amplicons.sort(key=lambda a: (a.subject, a.start))
    return amplicons


def nearest_size_gap(size: int, others: Sequence[int]) -> Optional[int]:
    gaps = [abs(size - other) for other in others]
    return min(gaps) if gaps else None


def _conservative_intended_products(amplicons: Sequence[Amplicon]) -> List[Amplicon]:
    candidates = [a for a in amplicons if a.on_target]
    if len(candidates) <= 1:
        return candidates
    for candidate in candidates[1:]:
        candidate.on_target = False
        candidate.__dict__["ambiguous_intended_duplicate"] = True
    candidates[0].__dict__["generic_intended_candidate"] = True
    return candidates[:1]


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
    sites, viable_sites = annotate_thermo(sites, primers, genome, thermo_params, thermo_gate)
    amplicons = enumerate_amplicons(sites, sp)

    sizes = [a.size for a in amplicons]
    for i, amplicon in enumerate(amplicons):
        amplicon.__dict__["nearest_gap"] = nearest_size_gap(
            amplicon.size, [size for j, size in enumerate(sizes) if j != i])

    return {
        "db": db,
        "primers": dict(primers),
        "sites_per_primer": {name: sum(s.primer == name for s in sites) for name in primers},
        "thermo_evaluated": bool(viable_sites) or (genome is not None and thermo_gate is False),
        "viable_sites_per_primer": viable_sites,
        "raw_hits_per_primer": {k: v.raw_blast_hits for k, v in hit_stats.items()},
        "unique_subjects_per_primer": {k: v.unique_subjects for k, v in hit_stats.items()},
        "near_blast_limit": [k for k, v in hit_stats.items() if v.near_target_limit],
        "high_copy_primers": [k for k, v in hit_stats.items() if v.high_copy],
        "blast_limits": {
            "max_target_seqs": sp.max_target_seqs,
            "evalue": sp.evalue,
            "word_size": sp.word_size,
            "num_threads": sp.num_threads,
            "high_copy_hit_threshold": sp.high_copy_hit_threshold,
        },
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
) -> Dict:
    sp = sp or SpecParams()
    blastn = _detect_blastn(blastn_bin)
    primers = {"F": forward, "R": reverse}
    sites, hit_stats = screen_primers_with_stats(primers, db, sp, blastn)
    sites, viable_sites = annotate_thermo(sites, primers, genome, thermo_params, thermo_gate)
    amplicons = enumerate_amplicons(sites, sp)

    for amplicon in amplicons:
        perfect = amplicon.fwd_mismatch == 0 and amplicon.rev_mismatch == 0
        proper_pair = {amplicon.fwd_primer, amplicon.rev_primer} == {"F", "R"}
        size_ok = designed_size is None or abs(amplicon.size - designed_size) <= size_tolerance
        amplicon.on_target = perfect and proper_pair and size_ok
    _conservative_intended_products(amplicons)

    on = [a for a in amplicons if a.on_target]
    off = [a for a in amplicons if not a.on_target]
    ref_size = designed_size if designed_size is not None else (on[0].size if on else None)
    comigrating = ([a for a in off if abs(a.size - ref_size) < sp.gel_min_gap_bp]
                   if ref_size is not None else [])
    nearest_off_gap = (nearest_size_gap(ref_size, [a.size for a in off])
                       if ref_size is not None and off else None)

    return {
        "db": db,
        "n_forward_sites": sum(s.primer == "F" for s in sites),
        "n_reverse_sites": sum(s.primer == "R" for s in sites),
        "thermo_evaluated": bool(viable_sites) or (genome is not None and thermo_gate is False),
        "viable_sites_per_primer": viable_sites,
        "raw_hits_per_primer": {k: v.raw_blast_hits for k, v in hit_stats.items()},
        "unique_subjects_per_primer": {k: v.unique_subjects for k, v in hit_stats.items()},
        "near_blast_limit": [k for k, v in hit_stats.items() if v.near_target_limit],
        "high_copy_primers": [k for k, v in hit_stats.items() if v.high_copy],
        "blast_limits": {
            "max_target_seqs": sp.max_target_seqs,
            "evalue": sp.evalue,
            "word_size": sp.word_size,
            "num_threads": sp.num_threads,
            "high_copy_hit_threshold": sp.high_copy_hit_threshold,
        },
        "n_products": len(amplicons),
        "n_on_target": len(on),
        "n_off_target": len(off),
        "n_comigrating": len(comigrating),
        "nearest_offtarget_gap": nearest_off_gap,
        "gel_distinguishable": len(comigrating) == 0,
        "on_target": on,
        "off_target": off,
        "specific": len(amplicons) == 1 and len(on) == 1,
    }
