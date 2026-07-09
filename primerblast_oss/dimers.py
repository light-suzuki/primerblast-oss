"""Primer-dimer, cross-dimer and hairpin analysis (optional, primer3-py).

Design-time filters (Primer3's per-oligo self/hairpin limits) are not enough for
two things NCBI Primer-BLAST does not address:

  * the **forward x reverse cross-dimer** of a chosen pair (especially a 3'-end
    dimer, which is what actually kills a PCR), and
  * **multiplex compatibility**: every primer against every other primer in a
    pool, to pick sets that can be run together.

Uses primer3-py (calc_hairpin / calc_homodimer / calc_heterodimer, the same
thermodynamic model as Primer3 itself). Optional: `available()` is False without
primer3-py and callers simply skip dimer analysis.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import primer3 as _primer3
    _HAVE = True
except Exception:  # pragma: no cover
    _primer3 = None
    _HAVE = False


def available() -> bool:
    return _HAVE


@dataclass
class DimerParams:
    # a structure is concerning if it is strong (very negative ΔG) OR it forms at
    # the reaction temperature (high Tm) with at least a modest ΔG.
    dg_warn: float = -8.0       # ΔG (kcal/mol) strong enough to matter regardless of Tm
    tm_warn: float = 45.0       # Tm (degC) at/above which a weaker structure still forms
    dg_at_tm_warn: float = -4.0  # ΔG floor applied together with tm_warn
    mv_conc: float = 50.0
    dv_conc: float = 1.5
    dntp_conc: float = 0.6
    dna_conc: float = 50.0


@dataclass
class Structure:
    kind: str                   # "hairpin" | "self-dimer" | "cross-dimer"
    a: str                      # primer name(s) involved
    b: str
    tm: float
    dg: float                   # kcal/mol
    concerning: bool


def _kw(dp: DimerParams) -> Dict:
    return dict(mv_conc=dp.mv_conc, dv_conc=dp.dv_conc,
                dntp_conc=dp.dntp_conc, dna_conc=dp.dna_conc)


def _res(kind, a, b, thermo, dp: DimerParams) -> Structure:
    tm = float(thermo.tm)
    dg = float(thermo.dg) / 1000.0
    concerning = dg <= dp.dg_warn or (tm >= dp.tm_warn and dg <= dp.dg_at_tm_warn)
    return Structure(kind=kind, a=a, b=b, tm=round(tm, 1), dg=round(dg, 2),
                     concerning=concerning)


def hairpin(name: str, seq: str, dp: Optional[DimerParams] = None) -> Optional[Structure]:
    if not _HAVE:
        return None
    dp = dp or DimerParams()
    return _res("hairpin", name, name, _primer3.calc_hairpin(seq, **_kw(dp)), dp)


def self_dimer(name: str, seq: str, dp: Optional[DimerParams] = None) -> Optional[Structure]:
    if not _HAVE:
        return None
    dp = dp or DimerParams()
    return _res("self-dimer", name, name, _primer3.calc_homodimer(seq, **_kw(dp)), dp)


def cross_dimer(name_a: str, seq_a: str, name_b: str, seq_b: str,
                dp: Optional[DimerParams] = None) -> Optional[Structure]:
    if not _HAVE:
        return None
    dp = dp or DimerParams()
    return _res("cross-dimer", name_a, name_b,
                _primer3.calc_heterodimer(seq_a, seq_b, **_kw(dp)), dp)


def analyze_pair(forward: str, reverse: str,
                 dp: Optional[DimerParams] = None) -> Optional[Dict]:
    """Hairpin + self-dimer for each primer and the F x R cross-dimer."""
    if not _HAVE:
        return None
    dp = dp or DimerParams()
    structs = [
        hairpin("F", forward, dp), hairpin("R", reverse, dp),
        self_dimer("F", forward, dp), self_dimer("R", reverse, dp),
        cross_dimer("F", forward, "R", reverse, dp),
    ]
    structs = [s for s in structs if s is not None]
    worst = min((s.dg for s in structs), default=None)
    concerning = [s for s in structs if s.concerning]
    return {
        "structures": structs,
        "worst_dg": worst,
        "n_concerning": len(concerning),
        "cross_dimer_dg": next((s.dg for s in structs if s.kind == "cross-dimer"), None),
        "ok": len(concerning) == 0,
    }


def select_multiplex_set(candidates: Sequence[Tuple[str, Sequence[Tuple[str, str]]]],
                         dp: Optional[DimerParams] = None,
                         max_steps: int = 200000) -> Optional[Dict]:
    """Pick one primer pair per target so the whole set is multiplex-compatible.

    `candidates` is a list of (target_name, [(fwd, rev), ...]) -- several
    candidate pairs per target (e.g. the top design hits). Returns a chosen pair
    per target such that no forward/reverse primer of any target forms a
    concerning cross-dimer with a primer of another target. Uses backtracking to
    find a fully compatible set; if none exists within the step budget, falls
    back to a greedy partial selection and lists the targets left unassigned.

    NCBI Primer-BLAST does not do this: it designs each amplicon independently.
    """
    if not _HAVE:
        return None
    dp = dp or DimerParams()
    targets = [(name, list(pairs)) for name, pairs in candidates]
    n = len(targets)
    assignment: List[Optional[int]] = [None] * n
    chosen: List[Tuple[str, str]] = []          # committed (label, seq)
    conflicts: List[Structure] = []
    steps = [0]

    def compatible(new_primers: List[Tuple[str, str]]) -> bool:
        for la, sa in new_primers:
            for lb, sb in chosen:
                s = cross_dimer(la, sa, lb, sb, dp)
                if s is not None and s.concerning:
                    conflicts.append(s)
                    return False
        return True

    def dfs(i: int) -> bool:
        if i == n:
            return True
        name, pairs = targets[i]
        for idx, (fwd, rev) in enumerate(pairs):
            steps[0] += 1
            if steps[0] > max_steps:
                return False
            new_primers = [(f"{name}_F", fwd), (f"{name}_R", rev)]
            if compatible(new_primers):
                assignment[i] = idx
                chosen.extend(new_primers)
                if dfs(i + 1):
                    return True
                del chosen[-2:]
                assignment[i] = None
        return False

    complete = dfs(0)

    if not complete:
        # greedy partial: keep whatever we can place, skip the rest
        assignment = [None] * n
        chosen = []
        for i, (name, pairs) in enumerate(targets):
            for idx, (fwd, rev) in enumerate(pairs):
                new_primers = [(f"{name}_F", fwd), (f"{name}_R", rev)]
                if compatible(new_primers):
                    assignment[i] = idx
                    chosen.extend(new_primers)
                    break

    selection = []
    unassigned = []
    for i, (name, pairs) in enumerate(targets):
        idx = assignment[i]
        if idx is None:
            unassigned.append(name)
            selection.append({"target": name, "pair": None, "candidate_index": None})
        else:
            fwd, rev = pairs[idx]
            selection.append({"target": name, "forward": fwd, "reverse": rev,
                              "candidate_index": idx})
    return {
        "n_targets": n,
        "complete": complete and not unassigned,
        "selection": selection,
        "unassigned": unassigned,
        "n_assigned": n - len(unassigned),
    }


def analyze_multiplex(primers: Sequence[Tuple[str, str]],
                      dp: Optional[DimerParams] = None) -> Optional[Dict]:
    """Every-primer-vs-every-primer cross-dimer scan for a multiplex pool.

    `primers` is a list of (name, sequence). Returns the concerning cross-dimers
    (ΔG <= dg_warn), sorted worst first -- the set NCBI Primer-BLAST won't check.
    """
    if not _HAVE:
        return None
    dp = dp or DimerParams()
    pairs: List[Structure] = []
    for (na, sa), (nb, sb) in combinations(primers, 2):
        s = cross_dimer(na, sa, nb, sb, dp)
        if s is not None:
            pairs.append(s)
    pairs.sort(key=lambda s: s.dg)
    concerning = [s for s in pairs if s.concerning]
    return {
        "n_primers": len(primers),
        "n_pairs_checked": len(pairs),
        "concerning": concerning,
        "n_concerning": len(concerning),
        "compatible": len(concerning) == 0,
        "worst_dg": pairs[0].dg if pairs else None,
    }
