#!/usr/bin/env python3
"""
parse_gbk.py - Extract CDS features from a GenBank flat file.

Usage: python3 parse_gbk.py input.gbk output.json [--include-partial]

Each CDS is reported with its strand and its biological 5' and 3' coordinates:
on the + strand the 5' end is the lowest coordinate, on the - strand it is the
highest. A join() location (a gene spanning the origin of a circular replicon)
is kept as a list of segments.

CDS features with partial coordinates (< or >) have no annotated start or stop,
so their boundaries cannot be scored. They are skipped unless --include-partial
is given.
"""

import json
import re
import sys


def read_cds_locations(gbk_path):
    """Return a list of (location_string, qualifiers) for every CDS feature."""
    features = []
    in_features = False
    current = None                      # [location, qualifier_text] of the CDS being read

    with open(gbk_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith("ORIGIN"):
                break
            if not in_features:
                continue

            key = line[5:21].strip()    # feature key column
            value = line[21:].rstrip("\n")

            if key:                     # a new feature starts on this line
                if current is not None:
                    features.append(current)
                current = [value, ""] if key == "CDS" else None
            elif current is not None:
                if value.startswith("/") or current[1]:
                    current[1] += value + "\n"      # qualifier lines
                else:
                    current[0] += value             # location continued on next line

    if current is not None:
        features.append(current)
    return features


def get_qualifier(qualifier_text, name):
    """Return the value of one qualifier (e.g. locus_tag), or an empty string."""
    match = re.search(r'/' + name + r'="([^"]*)"', qualifier_text)
    return match.group(1).replace("\n", " ") if match else ""


def parse_location(location):
    """Turn a GenBank location string into strand, segments, 5' and 3' ends."""
    strand = "-" if location.startswith("complement") else "+"
    partial = "<" in location or ">" in location
    pairs = re.findall(r"[<>]?(\d+)\.\.[<>]?(\d+)", location)
    segments = [[int(a), int(b)] for a, b in pairs]

    if strand == "+":
        five_prime = segments[0][0]     # first base of the first segment
        three_prime = segments[-1][1]   # last base of the last segment
    else:
        five_prime = segments[-1][1]
        three_prime = segments[0][0]

    return {
        "strand": strand,
        "segments": segments,
        "five_prime": five_prime,
        "three_prime": three_prime,
        "partial": partial,
    }


def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: python3 parse_gbk.py input.gbk output.json [--include-partial]")
    include_partial = "--include-partial" in sys.argv

    genes = []
    skipped = 0
    for location, qualifiers in read_cds_locations(sys.argv[1]):
        record = parse_location(location)
        if record["partial"] and not include_partial:
            skipped += 1
            continue
        record["location"] = location
        record["locus_tag"] = get_qualifier(qualifiers, "locus_tag")
        record["gene"] = get_qualifier(qualifiers, "gene")
        record["product"] = get_qualifier(qualifiers, "product")
        genes.append(record)

    with open(sys.argv[2], "w", encoding="utf-8") as out:
        json.dump(genes, out, indent=2)

    print(f"CDS features kept: {len(genes)}   partial CDS skipped: {skipped}")


if __name__ == "__main__":
    main()
