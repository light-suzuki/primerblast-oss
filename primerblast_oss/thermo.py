"""Optional thermodynamic scoring of primer binding, via primer3-py.

A BLAST hit that passes the mismatch / 3'-anchor rule may still be a weak
primer in practice — or a genuine mispriming site. primer3-py (the same
SantaLucia nearest-neighbour model used by Primer3, NCBI Primer-BLAST and
PrimerServer2) lets us score each candidate site by:

  * the duplex melting temperature (overall binding strength), and
  * the 3'-end stability ΔG (whether the 3' end anneals well enough to extend).

Both matter: a single 3'-terminal mismatch barely changes the duplex Tm but
sharply weakens the 3'-end ΔG, so Tm alone is not enough.

This module is OPTIONAL. If primer3-py is not installed, `available()` returns
False and callers fall back to the mismatch / 3'-anchor model.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

try:
    import primer3 as _primer3
    _HAVE = True
except Exception:  # pragma: no cover - environment without primer3-py
    _primer3 = None
    _HAVE = False


def available() -> bool:
    return _HAVE


@dataclass
class ThermoParams:
    """Reaction conditions and viability thresholds (defaults calibrated on a
    60 degC primer: perfect Tm~60/3'dG~-23; random Tm<0/3'dG~0)."""

    min_anneal_tm: float = 40.0     # duplex Tm floor for a site to prime (degC)
    max_3p_dg: float = -5.0         # 3'-end ΔG must be at least this stable (kcal/mol)
    mv_conc: float = 50.0           # monovalent cation (mM)
    dv_conc: float = 1.5            # divalent cation (mM)
    dntp_conc: float = 0.6          # dNTP (mM)
    dna_conc: float = 50.0          # primer/oligo (nM)


@dataclass
class SiteThermo:
    tm: float                       # duplex melting temperature (degC)
    duplex_dg: float                # duplex ΔG (kcal/mol)
    end3_dg: float                  # 3'-end stability ΔG (kcal/mol)
    viable: bool                    # would this site realistically prime?


def evaluate(primer: str, target_bind_strand: str,
             tp: Optional[ThermoParams] = None) -> Optional[SiteThermo]:
    """Score a primer against the strand it anneals to (given 5'->3').

    Returns None if primer3-py is unavailable or the target is empty.
    """
    if not _HAVE or not target_bind_strand:
        return None
    tp = tp or ThermoParams()
    kw = dict(mv_conc=tp.mv_conc, dv_conc=tp.dv_conc,
              dntp_conc=tp.dntp_conc, dna_conc=tp.dna_conc)
    h = _primer3.calc_heterodimer(primer, target_bind_strand, **kw)
    e = _primer3.calc_end_stability(primer, target_bind_strand, **kw)
    tm = float(h.tm)
    duplex_dg = float(h.dg) / 1000.0
    end3_dg = float(e.dg) / 1000.0
    viable = tm >= tp.min_anneal_tm and end3_dg <= tp.max_3p_dg
    return SiteThermo(tm=round(tm, 1), duplex_dg=round(duplex_dg, 2),
                      end3_dg=round(end3_dg, 2), viable=viable)
