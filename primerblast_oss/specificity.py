"""Specificity / in-silico PCR: the core Primer-BLAST-equivalent step.

We BLAST each primer, turn hits into 3'-anchored *priming sites*, then pair a
plus-strand site with a downstream minus-strand site on the same subject to
enumerate every predicted PCR product within the size window -- regardless of
which primer plays forward or reverse (F/R, R/F, F/F, R/R). This is exactly
what NCBI Primer-BLAST does and what a per-primer BLAST check omits.

The primer pool is generic: pass any number of named primers and get back all
products they would make together, each with its size -- so a caller can judge
tolerability from the size separation rather than a hard specific/non-specific
verdict.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

_OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq sseq qlen"


@dataclass
class SpecParams:
    """Thresholds controlling when a BLAST hit is treated as a priming site."""

    max_total_mismatch: int = 4      # over the full primer
    max_3prime_mismatch: int = 1     # within the 3'-terminal window
    three_prime_window: int = 5
    require_3prime_terminal_match: bool = True
    min_product: int = 40
    max_product: int = 4000          # off-target amplicons larger than this don't amplify
    gel_min_gap_bp: int = 50         # size gap needed to resolve two products on a gel
    # blastn-short tuning
    word_size: int = 7
    evalue: float = 30000.0
    max_target_seqs: int = 5000
    dust: str = "no"
    reward: int = 1
    penalty: int = -1
    gapopen: int = 5
    gapextend: int = 2
    num_threads: int = 4             # blastn threads (throughput, per PrimerServer2)


@dataclass
class PrimingSite:
    primer: str             # primer name/id that anneals here
    subject: str
    strand: str             # "+" primes toward increasing coord, "-" toward decreasing
    end3: int               # 1-based subject coordinate of the primer's 3' base
    total_mismatch: int
    tp_mismatch: int        # mismatches within the 3' window
    plen: int = 20          # primer length, needed to place its 5' end
    tp5_mismatch: int = 0   # mismatches within the 3'-terminal 5 bases
    tp10_mismatch: int = 0  # mismatches within the 3'-terminal 10 bases

    @property
    def end5(self) -> int:
        """1-based subject coordinate of the primer's 5' base."""
        return self.end3 - (self.plen - 1) if self.strand == "+" else self.end3 + (self.plen - 1)

    @property
    def extends(self) -> str:
        return "right" if self.strand == "+" else "left"

    def describe(self) -> str:
        return (f"{self.primer} binds {self.subject} {self.strand} strand, "
                f"5'={self.end5} 3'={self.end3}, extends {self.extends}, "
                f"mm total={self.total_mismatch} 3'5bp={self.tp5_mismatch}")


@dataclass
class Amplicon:
    subject: str
    start: int              # 1-based, 5' end of the plus-strand primer
    end: int                # 1-based, 5' end of the minus-strand primer
    size: int
    fwd_primer: str         # primer priming on the plus strand
    rev_primer: str         # primer priming on the minus strand
    fwd_mismatch: int
    rev_mismatch: int
    on_target: bool = False
    fwd_tp5: int = 0        # forward primer 3'-5bp mismatches at this site
    rev_tp5: int = 0        # reverse primer 3'-5bp mismatches at this site

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
    """Count mismatches within the last `window` aligned query bases and report
    whether the 3'-terminal query base matches. Gaps count as mismatch."""
    mm = 0
    seen = 0
    terminal_match = False
    first = True
    for q, s in zip(reversed(qseq), reversed(sseq)):
        if q == "-":
            mm += 1  # gap right at/near the 3' end kills priming
            continue
        if first:
            terminal_match = (q == s)
            first = False
        if q != s:
            mm += 1
        seen += 1
        if seen >= window:
            break
    return mm, terminal_match


def _hit_to_site(fields: List[str], primer_id: str, sp: SpecParams) -> Optional[PrimingSite]:
    (_qid, sseqid, _pident, length, mismatch, gapopen, qstart, qend,
     sstart, send, _e, _bits, sstrand, qseq, sseq, qlen) = fields
    qend_i = int(qend)
    qlen_i = int(qlen)
    total_mm = int(mismatch) + int(gapopen)
    # 3' end of the primer must be part of the alignment, or it cannot anneal.
    if qend_i != qlen_i:
        return None
    tp_mm, terminal_match = _count_3prime_mismatch(qseq, sseq, sp.three_prime_window)
    if sp.require_3prime_terminal_match and not terminal_match:
        return None
    # account for the 5' portion of the primer not covered by the HSP: those
    # unaligned bases are treated as mismatches against the full primer.
    total_mm += int(qstart) - 1
    if total_mm > sp.max_total_mismatch:
        return None
    if tp_mm > sp.max_3prime_mismatch:
        return None
    send_i = int(send)
    strand = "+" if sstrand == "plus" else "-"
    tp5, _ = _count_3prime_mismatch(qseq, sseq, 5)
    tp10, _ = _count_3prime_mismatch(qseq, sseq, 10)
    return PrimingSite(
        primer=primer_id, subject=sseqid, strand=strand, end3=send_i,
        total_mismatch=total_mm, tp_mismatch=tp_mm, plen=qlen_i,
        tp5_mismatch=tp5, tp10_mismatch=tp10,
    )


def priming_sites(primer: str, primer_id: str, db: str, sp: SpecParams, blastn: str) -> List[PrimingSite]:
    out = _run_blast(primer, db, sp, blastn)
    sites: List[PrimingSite] = []
    for line in out.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 16:
            continue
        site = _hit_to_site(fields, primer_id, sp)
        if site is not None:
            sites.append(site)
    return sites


def screen_primers(primers: Dict[str, str], db: str, sp: SpecParams, blastn: str) -> List[PrimingSite]:
    """BLAST every named primer and return the pooled priming sites."""
    sites: List[PrimingSite] = []
    for name, seq in primers.items():
        sites.extend(priming_sites(seq, name, db, sp, blastn))
    return sites


def enumerate_amplicons(sites: Sequence[PrimingSite], sp: SpecParams) -> List[Amplicon]:
    """Pair every plus-strand site with a downstream minus-strand site on the
    same subject within the product-size window. Orientation-agnostic: any
    primer may act as forward or reverse."""
    by_subject: Dict[str, List[PrimingSite]] = {}
    for s in sites:
        by_subject.setdefault(s.subject, []).append(s)

    amplicons: List[Amplicon] = []
    for subject, group in by_subject.items():
        plus = sorted([s for s in group if s.strand == "+"], key=lambda x: x.end5)
        minus = sorted([s for s in group if s.strand == "-"], key=lambda x: x.end5)
        for f in plus:
            for r in minus:
                # the 3' ends must converge (f upstream of r); the PCR product
                # spans from the plus primer's 5' end to the minus primer's 5' end.
                if r.end3 < f.end3:
                    continue
                start = f.end5           # 5' of plus-strand primer
                end = r.end5             # 5' of minus-strand primer
                size = end - start + 1
                if size < sp.min_product:
                    continue
                if size > sp.max_product:
                    break  # minus sorted by end5 ascending; further ones only larger
                amplicons.append(Amplicon(
                    subject=subject, start=start, end=end, size=size,
                    fwd_primer=f.primer, rev_primer=r.primer,
                    fwd_mismatch=f.total_mismatch, rev_mismatch=r.total_mismatch,
                    fwd_tp5=f.tp5_mismatch, rev_tp5=r.tp5_mismatch,
                ))
    amplicons.sort(key=lambda a: (a.subject, a.start))
    return amplicons


def nearest_size_gap(size: int, others: Sequence[int]) -> Optional[int]:
    """Smallest absolute size difference between `size` and any other product."""
    gaps = [abs(size - o) for o in others]
    return min(gaps) if gaps else None


def in_silico_pcr(
    primers: Dict[str, str],
    db: str,
    sp: Optional[SpecParams] = None,
    blastn_bin: Optional[str] = None,
) -> Dict:
    """Predict every PCR product a pool of primers makes against one database.

    Orientation is not constrained: the result lists all products with their
    sizes so the caller can judge whether extra products are resolvable.
    """
    sp = sp or SpecParams()
    blastn = _detect_blastn(blastn_bin)
    sites = screen_primers(primers, db, sp, blastn)
    amplicons = enumerate_amplicons(sites, sp)

    sizes = [a.size for a in amplicons]
    for a in amplicons:
        others = [s for i, s in enumerate(sizes) if amplicons[i] is not a]
        gap = nearest_size_gap(a.size, others)
        a.__dict__["nearest_gap"] = gap  # attach for reporting/GUI

    per_primer = {name: sum(1 for s in sites if s.primer == name) for name in primers}
    return {
        "db": db,
        "primers": {k: v for k, v in primers.items()},
        "sites_per_primer": per_primer,
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
) -> Dict:
    """Specificity analysis for one designed primer pair against one database.

    An amplicon is on_target when both primers anneal perfectly, form a proper
    forward+reverse pair, and match the designed size. Everything else is an
    off-target product. Off-targets are further judged by whether their size is
    resolvable from the intended product on a gel (gel_min_gap_bp).
    """
    sp = sp or SpecParams()
    blastn = _detect_blastn(blastn_bin)
    sites = screen_primers({"F": forward, "R": reverse}, db, sp, blastn)
    amplicons = enumerate_amplicons(sites, sp)

    for a in amplicons:
        perfect = a.fwd_mismatch == 0 and a.rev_mismatch == 0
        proper_pair = {a.fwd_primer, a.rev_primer} == {"F", "R"}
        size_ok = designed_size is None or abs(a.size - designed_size) <= size_tolerance
        a.on_target = perfect and proper_pair and size_ok

    on = [a for a in amplicons if a.on_target]
    off = [a for a in amplicons if not a.on_target]

    # gel resolvability: off-targets whose size is close to the intended product
    # would co-migrate and confound the result; far-sized ones are tolerable.
    ref_size = designed_size if designed_size is not None else (on[0].size if on else None)
    comigrating = []
    if ref_size is not None:
        comigrating = [a for a in off if abs(a.size - ref_size) < sp.gel_min_gap_bp]
    nearest_off_gap = nearest_size_gap(ref_size, [a.size for a in off]) if (ref_size is not None and off) else None

    return {
        "db": db,
        "n_forward_sites": sum(1 for s in sites if s.primer == "F"),
        "n_reverse_sites": sum(1 for s in sites if s.primer == "R"),
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
