"""GFF3 evidence for genomic template views; no inferred gene structures."""
from .annotation_index import region_annotation


def template_annotations(path, context):
    if not path:
        return {"status": "not_provided", "genes": []}
    chrom = context.get("chrom")
    if not chrom:
        return {"status": "no_genomic_coordinates", "genes": []}
    low, high = sorted((context["start"], context["end"]))
    annotation, exists = region_annotation(path, chrom, low, high)
    if not exists:
        return {"status": "seqid_not_found", "source": path, "genes": []}
    genes = []
    for gene in annotation.features_in(chrom, low, high, types=["gene"]):
        transcripts = []
        transcript_features = [child for child in gene.children if child.children]
        for transcript in transcript_features or [gene]:
            segments = [child for child in transcript.children
                        if child.type.lower() in ("exon", "cds")]
            if not segments:
                continue
            transcripts.append({
                "id": transcript.id or transcript.name or gene.id,
                "start": transcript.start, "end": transcript.end,
                "segments": [{"type": feature.type.lower(), "start": feature.start,
                              "end": feature.end} for feature in segments],
            })
        genes.append({"id": gene.id, "name": gene.name or gene.id,
                      "start": gene.start, "end": gene.end, "strand": gene.strand,
                      "transcripts": transcripts})
    return {"status": "loaded", "source": path, "genes": genes}
