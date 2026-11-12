#!/usr/bin/env python3
"""
compare_ref_vs_glimmer.py - Score Glimmer3 predictions against reference CDS.

Usage: python3 compare_ref_vs_glimmer.py reference.json glimmer.json results_dir

Each reference CDS is paired with the Glimmer3 prediction on the same strand
that overlaps it by the most bases. The pair is then classified by whether the
5' (start codon) and 3' (stop codon) coordinates agree. Writes summary.json and
detailed.json into results_dir.
"""

import json
import os
import sys


def overlap_bases(segments_a, segments_b):
    """Count bases shared between two lists of [low, high] segments."""
    total = 0
    for a_low, a_high in segments_a:
        for b_low, b_high in segments_b:
            shared = min(a_high, b_high) - max(a_low, b_low) + 1
            if shared > 0:
                total += shared
    return total


def best_prediction(gene, predictions):
    """Return the same-strand prediction with the largest overlap, or None."""
    best, best_overlap = None, 0
    for prediction in predictions:
        if prediction["strand"] != gene["strand"]:
            continue
        shared = overlap_bases(gene["segments"], prediction["segments"])
        if shared > best_overlap:
            best, best_overlap = prediction, shared
    return best


def classify(gene, prediction):
    """Label one reference/prediction pair by 5' and 3' agreement."""
    if prediction is None:
        return "no_prediction"
    five = gene["five_prime"] == prediction["five_prime"]
    three = gene["three_prime"] == prediction["three_prime"]
    if five and three:
        return "exact_match"
    if five:
        return "5_agree_3_disagree"
    if three:
        return "5_disagree_3_agree"
    return "both_disagree"


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: python3 compare_ref_vs_glimmer.py reference.json glimmer.json results_dir")

    with open(sys.argv[1], encoding="utf-8") as handle:
        reference = json.load(handle)
    with open(sys.argv[2], encoding="utf-8") as handle:
        predictions = json.load(handle)

    detailed = []
    counts = {}
    for gene in reference:
        prediction = best_prediction(gene, predictions)
        category = classify(gene, prediction)
        counts[category] = counts.get(category, 0) + 1
        detailed.append({
            "locus_tag": gene["locus_tag"],
            "gene": gene["gene"],
            "product": gene["product"],
            "strand": gene["strand"],
            "ref_5_prime": gene["five_prime"],
            "ref_3_prime": gene["three_prime"],
            "glimmer_orf": prediction["orf_id"] if prediction else "",
            "glimmer_5_prime": prediction["five_prime"] if prediction else "",
            "glimmer_3_prime": prediction["three_prime"] if prediction else "",
            "category": category,
        })

    # Predictions that share no base with any reference CDS on the same strand.
    no_overlap = 0
    for prediction in predictions:
        hits = [g for g in reference
                if g["strand"] == prediction["strand"]
                and overlap_bases(g["segments"], prediction["segments"]) > 0]
        if not hits:
            no_overlap += 1

    summary = {
        "reference_cds_count": len(reference),
        "glimmer_predicted_count": len(predictions),
        "exact_match": counts.get("exact_match", 0),
        "5_agree_3_disagree": counts.get("5_agree_3_disagree", 0),
        "5_disagree_3_agree": counts.get("5_disagree_3_agree", 0),
        "both_disagree": counts.get("both_disagree", 0),
        "reference_with_no_prediction": counts.get("no_prediction", 0),
        "predictions_with_no_reference_overlap": no_overlap,
    }

    os.makedirs(sys.argv[3], exist_ok=True)
    with open(os.path.join(sys.argv[3], "summary.json"), "w", encoding="utf-8") as out:
        json.dump(summary, out, indent=2)
    with open(os.path.join(sys.argv[3], "detailed.json"), "w", encoding="utf-8") as out:
        json.dump(detailed, out, indent=2)

    for key, value in summary.items():
        print(f"{key:40s} {value}")


if __name__ == "__main__":
    main()
