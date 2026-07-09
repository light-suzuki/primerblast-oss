"""Tiling design: cover an entire region with overlapping amplicons.

NCBI Primer-BLAST designs one amplicon around a single target and tends to
cluster primers on one side of a long template. For sequencing or scanning a
whole gene you want a *series* of overlapping amplicons spanning the region.
This module walks left-to-right, forcing each successive amplicon into the
next window and preferring specific / gel-resolvable pairs.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Dict, List, Optional, Sequence, Tuple

from .design import DesignParams, PrimerPair, clean_sequence, design_primers
from .specificity import SpecParams, pair_specificity
from .pipeline import _score_pair


def _evaluate(pair: PrimerPair, databases: Sequence[str], sp: SpecParams,
              blastn_bin: Optional[str], size_tolerance: int,
              genome=None, thermo_params=None, thermo_gate: bool = True,
              dimer_params=None) -> None:
    per_db = [
        pair_specificity(pair.forward, pair.reverse, db,
                         designed_size=pair.product_size, sp=sp,
                         blastn_bin=blastn_bin, size_tolerance=size_tolerance,
                         genome=genome, thermo_params=thermo_params,
                         thermo_gate=thermo_gate)
        for db in databases
    ]
    from . import dimers as _dimers
    dimer = (_dimers.analyze_pair(pair.forward, pair.reverse, dimer_params)
             if _dimers.available() else None)
    _score_pair(pair, per_db, dimer)


def design_tiling(
    template_id: str,
    sequence: str,
    databases: Sequence[str],
    region: Optional[Tuple[int, int]] = None,   # 0-based inclusive [start, end]
    amplicon_min: int = 400,
    amplicon_max: int = 800,
    overlap: int = 40,
    design_params: Optional[DesignParams] = None,
    spec_params: Optional[SpecParams] = None,
    primer3_bin: Optional[str] = None,
    blastn_bin: Optional[str] = None,
    size_tolerance: int = 10,
    candidates_per_tile: int = 8,
    max_tiles: int = 200,
    genome=None,
    thermo_params=None,
    thermo_gate: bool = True,
    dimer_params=None,
) -> List[Dict]:
    seq = clean_sequence(sequence)
    L = len(seq)
    r0, r1 = region if region else (0, L - 1)
    r0 = max(0, r0)
    r1 = min(L - 1, r1)

    dp_base = design_params or DesignParams()
    sp = spec_params or SpecParams()

    tiles: List[Dict] = []
    cur = r0
    prev_right: Optional[int] = None
    covered_to = r0 - 1

    while covered_to < r1 and len(tiles) < max_tiles:
        win_start = cur
        # near the region end, pull the window back so the final amplicon can
        # still reach r1 instead of leaving an uncovered tail.
        final_window = (r1 - win_start + 1) < amplicon_max
        if final_window:
            win_start = max(r0, r1 - amplicon_max + 1)
        win_len = min(amplicon_max, r1 - win_start + 1)
        if win_len < amplicon_min:
            break  # not enough template left for a full amplicon

        # for a short final window, let the product shrink to fit it
        amp_lo = min(amplicon_min, win_len)
        dp = replace(
            dp_base,
            product_size_ranges=[(amp_lo, win_len)],
            included_region=(win_start, win_len),
            target=None,
            num_return=candidates_per_tile,
        )
        try:
            pairs, _ = design_primers(template_id, seq, dp, primer3_bin)
        except Exception:
            pairs = []
        if not pairs:
            # nothing fits here; jump forward by most of a window and retry
            if win_start >= r1 - amplicon_min:
                break
            cur = win_start + max(1, amplicon_min // 2)
            continue

        for pair in pairs:
            _evaluate(pair, databases, sp, blastn_bin, size_tolerance,
                      genome=genome, thermo_params=thermo_params,
                      thermo_gate=thermo_gate, dimer_params=dimer_params)

        # prefer specific, then gel-resolvable; then position: normally the
        # leftmost amplicon (walks coverage forward), but in the final
        # end-anchored window the rightmost one (reaches the region end).
        def key(p: PrimerPair):
            s = p.specificity
            pos = -p.right_start if final_window else p.left_start
            return (
                0 if s.get("specific_all_db") else 1,
                0 if s.get("gel_distinguishable") else 1,
                pos,
                -s.get("score", 0.0),
            )
        best = sorted(pairs, key=key)[0]

        gap_to_prev = None
        if prev_right is not None:
            # positive => overlap with previous amplicon; negative => a gap
            gap_to_prev = prev_right - best.left_start + 1

        tiles.append({
            "index": len(tiles) + 1,
            "pair": best,
            "covers": (best.left_start, best.right_start),
            "gap_to_prev": gap_to_prev,
        })

        prev_right = best.right_start
        covered_to = max(covered_to, best.right_start)
        if final_window or best.right_start >= r1 - overlap:
            break                          # reached (or did our best for) the end
        next_cur = best.right_start - overlap + 1
        if next_cur <= win_start:          # guarantee forward progress
            next_cur = win_start + max(1, amplicon_min // 2)
        cur = next_cur

    return tiles
