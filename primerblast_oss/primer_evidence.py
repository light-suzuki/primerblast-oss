"""Separate literal oligo predictions from alternative input hypotheses."""


def annotate_input_evidence(results, primers):
    """Add per-product direction and per-database literal-input evidence.

    F/R are names, not strand constraints. Only changes to participating oligos
    affect whether a product is a prediction for the supplied sequences.
    """
    originals = {row["db"]: row for row in results
                 if not row["reverse_complemented_inputs"]}
    literal_loci = {
        db: {(p["subject"], p["start"], p["end"],
              tuple(sorted([p["fwd_primer"], p["rev_primer"]]))) for p in row["products"]}
        for db, row in originals.items()
    }
    assessments = []
    for row in results:
        baseline = originals[row["db"]]
        for product in row["products"]:
            names = [product["fwd_primer"], product["rev_primer"]]
            changed = list(dict.fromkeys(name for name in names
                                        if row["oligos"].get(name) != primers.get(name)))
            if names == ["R", "F"]:
                labels = "swapped"
            elif names[0] == names[1]:
                labels = "same_primer"
            elif names == ["F", "R"]:
                labels = "as_labeled"
            else:
                labels = "custom"
            product["input_evidence"] = {
                "sequence_status": "reverse_complement_candidate" if changed else "as_supplied",
                "changed_primers": changed,
                "label_orientation": labels,
                "original_search_complete": baseline["search_complete"],
                "as_supplied_locus_observed": (
                    product["subject"], product["start"], product["end"], tuple(sorted(names)))
                    in literal_loci[row["db"]],
            }
            product["primer_bindings"] = [
                {"primer": name, "input_sequence": primers.get(name),
                 "oligo": row["oligos"].get(name), "reference_strand": strand,
                 "extends": direction, "end5": product[end5],
                 "end3": product.get(end3), "mismatches": product.get(mm)}
                for name, strand, direction, end5, end3, mm in (
                    (names[0], "+", "right", "start", "fwd_end3", "fwd_mismatch"),
                    (names[1], "-", "left", "end", "rev_end3", "rev_mismatch"))
            ]
    for db, baseline in originals.items():
        products = baseline["products"]
        alternatives = {
            (p["subject"], p["start"], p["end"], p["fwd_primer"], p["rev_primer"],
             row["oligos"].get(p["fwd_primer"]), row["oligos"].get(p["rev_primer"]))
            for row in results if row["db"] == db
            for p in row["products"] if p["input_evidence"]["changed_primers"]
        }
        assessments.append({
            "db": db, "as_supplied_products": len(products),
            "distinct_primer_products": sum(p["fwd_primer"] != p["rev_primer"] for p in products),
            "same_primer_products": sum(p["fwd_primer"] == p["rev_primer"] for p in products),
            "reverse_complement_candidates": len(alternatives),
            "distinct_primer_reverse_complement_candidates": sum(p[3] != p[4] for p in alternatives),
            "original_search_complete": baseline["search_complete"],
            "input_error": "undetermined",
        })
    return assessments
