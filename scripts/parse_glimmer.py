#!/usr/bin/env python3
"""
parse_glimmer.py - Parse a Glimmer3 .predict file.

Usage: python3 parse_glimmer.py input.predict output.json --length SEQ_LENGTH

Each prediction line holds: orf id, start, end, reading frame, score.
Glimmer3 writes the start codon position first, so on the - strand the first
coordinate is the larger one. A + strand gene whose start is larger than its
end (or a - strand gene whose start is smaller) wraps around the origin of a
circular sequence and is stored as two segments.
"""

import json
import sys


def parse_predict(predict_path, seq_length):
    """Return a list of predictions with strand, segments, 5' and 3' ends."""
    predictions = []

    with open(predict_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(">"):    # skip blanks and the header
                continue

            orf_id, start, end, frame, score = line.split()
            start, end = int(start), int(end)
            strand = "+" if frame.startswith("+") else "-"

            low, high = (start, end) if strand == "+" else (end, start)
            if low <= high:
                segments = [[low, high]]
            else:                                   # wraps around the origin
                segments = [[low, seq_length], [1, high]]

            predictions.append({
                "orf_id": orf_id,
                "strand": strand,
                "frame": frame,
                "score": float(score),
                "segments": segments,
                "five_prime": start,                # start codon side
                "three_prime": end,                 # stop codon side
            })

    return predictions


def main():
    if len(sys.argv) != 5 or sys.argv[3] != "--length":
        sys.exit("Usage: python3 parse_glimmer.py input.predict output.json --length SEQ_LENGTH")

    predictions = parse_predict(sys.argv[1], int(sys.argv[4]))
    with open(sys.argv[2], "w", encoding="utf-8") as out:
        json.dump(predictions, out, indent=2)

    print(f"Glimmer3 predictions: {len(predictions)}")


if __name__ == "__main__":
    main()
