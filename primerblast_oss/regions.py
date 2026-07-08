"""Resolve breeding/analysis targets into genomic regions and extract templates.

A target can be a gene (via GFF3), an explicit interval, a SNP to flank (for
CAPS/dCAPS), or a BED file. Each resolves to a `GenomicRegion` carrying the
chromosome, span, strand and name. Extracting a template also returns a genomic
*anchor* so predicted products can be tied back to real coordinates -- which is
what makes the intended-vs-off-target call precise instead of size-guessing.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .genome import Genome


@dataclass
class GenomicRegion:
    chrom: str
    start: int              # 1-based inclusive
    end: int
    strand: str = "+"
    name: str = "region"
    source: str = "interval"

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass
class Template:
    id: str
    seq: str
    region: GenomicRegion
    # extraction span in genome coords (includes flank)
    ext_start: int
    ext_end: int
    # anchor for mapping a 0-based template index -> genomic coordinate
    anchor_coord: int       # genomic coord of template index 0
    anchor_strand: str      # '+' -> genomic = anchor + i ; '-' -> anchor - i
    flank: int = 0
    extras: Dict = field(default_factory=dict)

    def to_genomic(self, local_index0: int) -> int:
        if self.anchor_strand == "-":
            return self.anchor_coord - local_index0
        return self.anchor_coord + local_index0


def extract_template(genome: Genome, region: GenomicRegion, flank: int = 0,
                     id_suffix: str = "") -> Template:
    """Extract the region ±flank from the genome, oriented 5'->3' on the
    region's strand."""
    ext_start = max(1, region.start - flank)
    ext_end = min(genome.length(region.chrom), region.end + flank)
    seq = genome.fetch(region.chrom, ext_start, ext_end, region.strand)
    if region.strand == "-":
        anchor_coord = ext_end          # template index 0 maps to the high coord
        anchor_strand = "-"
    else:
        anchor_coord = ext_start
        anchor_strand = "+"
    tid = f"{region.name}{id_suffix}"
    return Template(
        id=tid, seq=seq, region=region, ext_start=ext_start, ext_end=ext_end,
        anchor_coord=anchor_coord, anchor_strand=anchor_strand, flank=flank,
    )


# --------------------------------------------------------------------------- #
# resolvers
# --------------------------------------------------------------------------- #
def resolve_interval(chrom: str, start: int, end: int, strand: str = "+",
                     name: str = "region") -> GenomicRegion:
    return GenomicRegion(chrom, start, end, strand, name, source="interval")


def resolve_snp(chrom: str, pos: int, flank: int = 250, name: Optional[str] = None,
                strand: str = "+") -> GenomicRegion:
    """A region centered on a SNP position, for CAPS/dCAPS amplicon design."""
    return GenomicRegion(chrom, pos - flank, pos + flank, strand,
                         name or f"{chrom}_{pos}", source="snp")


def resolve_gene(gff3_path: str, gene_id: str, feature: str = "cds",
                 flank: int = 0, gff3_seqid: Optional[str] = None) -> GenomicRegion:
    """Resolve a gene id to a region using a GFF3. `feature` selects the span:
    'gene', 'mrna', 'exon' (union), or 'cds' (union)."""
    from .gff3 import parse_gff3   # provided by the gff3 module
    gff = parse_gff3(gff3_path) if gff3_seqid is None else parse_gff3(gff3_path, seqid=gff3_seqid)
    r = gff.gene_region(gene_id, feature=feature)
    if r is None:
        raise KeyError(f"gene '{gene_id}' (feature={feature}) not found in {gff3_path}")
    chrom, start, end, strand = r
    return GenomicRegion(chrom, start, end, strand, name=gene_id, source=f"gff3:{feature}")


def resolve_bed(bed_path: str) -> List[GenomicRegion]:
    from .bed import parse_bed
    out: List[GenomicRegion] = []
    for iv in parse_bed(bed_path):
        out.append(GenomicRegion(
            iv.chrom, iv.start + 1, iv.end,            # BED 0-based half-open -> 1-based incl
            iv.strand if iv.strand in ("+", "-") else "+",
            name=iv.name or f"{iv.chrom}:{iv.start}-{iv.end}", source="bed"))
    return out


def tile_interval(region: GenomicRegion, n_markers: int = 0,
                  spacing: int = 0) -> List[GenomicRegion]:
    """Split a large interval (e.g. a QTL) into evenly spaced marker anchor
    points. Give either n_markers (evenly spaced) or a fixed spacing (bp)."""
    L = region.length
    if n_markers > 0:
        step = max(1, L // (n_markers + 1))
    elif spacing > 0:
        step = spacing
    else:
        raise ValueError("give n_markers or spacing")
    points: List[GenomicRegion] = []
    pos = region.start + step
    idx = 1
    while pos < region.end and (n_markers == 0 or idx <= n_markers):
        points.append(GenomicRegion(region.chrom, pos, pos, region.strand,
                                    name=f"{region.name}_m{idx}", source="qtl-marker"))
        pos += step
        idx += 1
    return points
