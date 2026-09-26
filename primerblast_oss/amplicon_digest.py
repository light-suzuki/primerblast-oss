"""Reconstruct and restriction-digest predicted off-target PCR products."""
from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from .caps import ENZYME_METADATA, digest_fragment_sizes
from .genome import revcomp


def reconstruct_pcr_product(
    genome,
    product: Mapping,
    primers: Mapping[str, str],
) -> Tuple[Optional[str], Optional[str]]:
    """Return the plus-strand PCR product including primer-incorporated ends.

    Predicted products are stored from the plus-strand priming site's 5' end to
    the minus-strand priming site's 5' end. The left primer sequence is
    incorporated as written; the right primer contributes its reverse
    complement to the plus-strand product.
    """
    try:
        subject = str(product["subject"])
        start = int(product["start"])
        end = int(product["end"])
        orientation = str(product["orientation"])
    except (KeyError, TypeError, ValueError) as error:
        return None, "invalid_product_metadata:%s" % error

    roles = orientation.split("/")
    if len(roles) != 2:
        return None, "invalid_orientation:%s" % orientation
    left_name, right_name = roles
    left_primer = str(primers.get(left_name, "")).upper()
    right_primer = str(primers.get(right_name, "")).upper()
    if not left_primer or not right_primer:
        return None, "missing_primer_sequence:%s" % orientation
    if start > end:
        start, end = end, start
    try:
        template = str(genome.fetch(subject, start, end, "+")).upper()
    except Exception as error:  # noqa: BLE001
        return None, "genome_fetch_failed:%s" % error
    if not template:
        return None, "empty_genome_fetch"
    if len(left_primer) + len(right_primer) > len(template):
        return None, "primer_ends_overlap"

    # PCR products inherit primer-encoded bases. This matters for dCAPS, where
    # the engineered mismatch may itself create the restriction site.
    product_sequence = (
        left_primer
        + template[len(left_primer):len(template) - len(right_primer)]
        + revcomp(right_primer)
    )
    return product_sequence, None


def digest_offtarget_products(
    per_db_products: Sequence[Dict],
    genomes_by_db: Optional[Mapping[str, object]],
    primers: Mapping[str, str],
    enzyme_name: str,
) -> List[Dict]:
    """Attach exact off-target digest fragments to each database view.

    Databases without an associated FASTA remain explicit/unresolved. A clean
    database with zero off-targets is complete even without an associated FASTA.
    """
    enzyme = ENZYME_METADATA.get(enzyme_name)
    if enzyme is None:
        raise ValueError("unknown restriction enzyme %r" % enzyme_name)

    genomes = dict(genomes_by_db or {})
    annotated: List[Dict] = []
    for view in per_db_products:
        output = dict(view)
        database = str(view.get("db", ""))
        off_targets = [
            product for product in (view.get("products") or [])
            if not product.get("on_target")
        ]
        details: List[Dict] = []
        unresolved: List[Dict] = []
        background: List[int] = []
        genome = genomes.get(database)

        if off_targets and genome is None:
            for product in off_targets:
                unresolved.append({
                    "subject": product.get("subject"),
                    "start": product.get("start"),
                    "end": product.get("end"),
                    "size": product.get("size"),
                    "orientation": product.get("orientation"),
                    "reason": "no_associated_genome",
                })
        else:
            for product in off_targets:
                sequence, error = reconstruct_pcr_product(
                    genome, product, primers) if genome is not None else (None, None)
                if error or not sequence:
                    unresolved.append({
                        "subject": product.get("subject"),
                        "start": product.get("start"),
                        "end": product.get("end"),
                        "size": product.get("size"),
                        "orientation": product.get("orientation"),
                        "reason": error or "sequence_unavailable",
                    })
                    continue
                fragments = digest_fragment_sizes(sequence, enzyme)
                background.extend(fragments)
                details.append({
                    "subject": product.get("subject"),
                    "start": product.get("start"),
                    "end": product.get("end"),
                    "size": product.get("size"),
                    "orientation": product.get("orientation"),
                    "sequence_length": len(sequence),
                    "fragments": fragments,
                })

        complete = len(unresolved) == 0
        output["offtarget_digest"] = {
            "enzyme": enzyme_name,
            "complete": complete,
            "n_off_target_products": len(off_targets),
            "n_digested_products": len(details),
            "background_fragments": sorted(background, reverse=True),
            "products": details,
            "unresolved_products": unresolved,
            "model": "reconstructed_pcr_product_with_primer_incorporation",
        }
        annotated.append(output)
    return annotated
