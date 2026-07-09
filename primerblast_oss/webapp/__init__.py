"""Local, bilingual (EN/JA) web GUI for primerblast-oss.

Stdlib-only front end over the design / in-silico PCR / tiling / assay /
markers / makedb engine. Launch with ``python -m primerblast_oss.webapp``.
"""
from .server import serve, health, discover_databases

__all__ = ["serve", "health", "discover_databases"]
