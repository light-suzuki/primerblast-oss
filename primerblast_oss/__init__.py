"""primerblast-oss: a local, open-source reimplementation of NCBI Primer-BLAST.

Design PCR primers with Primer3 and evaluate their specificity by pairing
BLAST hits into predicted amplicons on each subject sequence -- exactly the
step that distinguishes Primer-BLAST from a plain "BLAST each primer" check.

Runs fully offline against local BLAST+ databases (including unpublished
genomes), and can screen against several databases at once.
"""

__version__ = "0.2.0"

from .design import PrimerPair, design_primers, read_fasta
from .specificity import (
    PrimingSite, Amplicon, pair_specificity, in_silico_pcr,
    enumerate_amplicons, screen_primers,
)
from .pipeline import run_pipeline, PipelineResult
from .tiling import design_tiling
from .genome import Genome, revcomp
from .regions import GenomicRegion, extract_template, resolve_gene, resolve_interval, resolve_snp, resolve_bed
from .assay import run_assay, run_batch, design_qtl_markers
from .risk import assess_risk

__all__ = [
    "__version__",
    "PrimerPair",
    "design_primers",
    "read_fasta",
    "PrimingSite",
    "Amplicon",
    "pair_specificity",
    "in_silico_pcr",
    "enumerate_amplicons",
    "screen_primers",
    "run_pipeline",
    "PipelineResult",
    "design_tiling",
    "Genome",
    "revcomp",
    "GenomicRegion",
    "extract_template",
    "resolve_gene",
    "resolve_interval",
    "resolve_snp",
    "resolve_bed",
    "run_assay",
    "run_batch",
    "design_qtl_markers",
    "assess_risk",
]
