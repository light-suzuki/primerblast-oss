"""Tiling design: cover an entire region with overlapping amplicons."""
from __future__ import annotations

from dataclasses import replace
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .design import DesignParams, PrimerPair, clean_sequence, design_primers
from .specificity import SpecParams, pair_specificity
from .pipeline import (
    _score_pair,
    resolve_genome_for_database,
    thermo_metadata,
)


def _evaluate(pair: PrimerPair, databases: Sequence[str], sp: SpecParams,
              blastn_bin: Optional[str], size_tolerance: int,
              genome=None, genomes_by_db: Optional[Mapping[str, object]] = None,
              thermo_params=None, thermo_gate: bool = True,
              dimer_params=None) -> None:
    per_db = []
    for database in databases:
        database_genome, association = resolve_genome_for_database(
            database, databases, genome=genome, genomes_by_db=genomes_by_db)
        result = pair_specificity(
            pair.forward,
            pair.reverse,
            database,
            designed_size=pair.product_size,
            sp=sp,
            blastn_bin=blastn_bin,
            size_tolerance=size_tolerance,
            genome=database_genome,
            thermo_params=thermo_params,
            thermo_gate=thermo_gate,
        )
        result.update(thermo_metadata(
            database_genome, thermo_params, thermo_gate, association))
        per_db.append(result)

    from . import dimers as dimer_module
    dimer = (
        dimer_module.analyze_pair(pair.forward, pair.reverse, dimer_params)
        if dimer_module.available() else None
    )
    _score_pair(pair, per_db, dimer)


def design_tiling(
    template_id: str,
    sequence: str,
    databases: Sequence[str],
    region: Optional[Tuple[int, int]] = None,
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
    genomes_by_db: Optional[Mapping[str, object]] = None,
    thermo_params=None,
    thermo_gate: bool = True,
    dimer_params=None,
) -> List[Dict]:
    seq = clean_sequence(sequence)
    sequence_length = len(seq)
    region_start, region_end = region if region else (0, sequence_length - 1)
    region_start = max(0, region_start)
    region_end = min(sequence_length - 1, region_end)

    base_design = design_params or DesignParams()
    specificity = spec_params or SpecParams()

    tiles: List[Dict] = []
    current = region_start
    previous_right: Optional[int] = None
    covered_to = region_start - 1

    while covered_to < region_end and len(tiles) < max_tiles:
        window_start = current
        final_window = (region_end - window_start + 1) < amplicon_max
        if final_window:
            window_start = max(region_start, region_end - amplicon_max + 1)
        window_length = min(amplicon_max, region_end - window_start + 1)
        if window_length < amplicon_min:
            break

        amplicon_low = min(amplicon_min, window_length)
        design = replace(
            base_design,
            product_size_ranges=[(amplicon_low, window_length)],
            included_region=(window_start, window_length),
            target=None,
            num_return=candidates_per_tile,
        )
        try:
            pairs, _explain = design_primers(
                template_id, seq, design, primer3_bin)
        except Exception:
            pairs = []
        if not pairs:
            if window_start >= region_end - amplicon_min:
                break
            current = window_start + max(1, amplicon_min // 2)
            continue

        for pair in pairs:
            _evaluate(
                pair,
                databases,
                specificity,
                blastn_bin,
                size_tolerance,
                genome=genome,
                genomes_by_db=genomes_by_db,
                thermo_params=thermo_params,
                thermo_gate=thermo_gate,
                dimer_params=dimer_params,
            )

        def sort_key(pair: PrimerPair):
            result = pair.specificity
            position = -pair.right_start if final_window else pair.left_start
            return (
                0 if result.get("specific_all_db") else 1,
                0 if result.get("gel_distinguishable") else 1,
                result.get("rank") == "I",
                position,
                -result.get("score", 0.0),
            )

        best = sorted(pairs, key=sort_key)[0]
        gap_to_previous = None
        if previous_right is not None:
            gap_to_previous = previous_right - best.left_start + 1

        tiles.append({
            "index": len(tiles) + 1,
            "pair": best,
            "covers": (best.left_start, best.right_start),
            "gap_to_prev": gap_to_previous,
        })

        previous_right = best.right_start
        covered_to = max(covered_to, best.right_start)
        if final_window or best.right_start >= region_end - overlap:
            break
        next_current = best.right_start - overlap + 1
        if next_current <= window_start:
            next_current = window_start + max(1, amplicon_min // 2)
        current = next_current

    return tiles
