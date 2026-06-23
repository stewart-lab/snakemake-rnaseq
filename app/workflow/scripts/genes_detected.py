#!/usr/bin/env python3
"""Count detected genes per sample from RSEM .genes.results files.

A gene is "detected" if its RSEM expected_count is at least --min-count (default
10). Writes a sample,genes_detected CSV used by the QC report card to set a
dataset-relative target for the genes-detected metric.
"""
import argparse
import csv
import os


def count_detected(path, min_count):
    n = 0
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                if float(row["expected_count"]) >= min_count:
                    n += 1
            except (KeyError, ValueError):
                continue
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="output CSV: sample,genes_detected")
    ap.add_argument("--min-count", type=float, default=10.0,
                    help="minimum expected_count to call a gene detected [default %(default)s]")
    ap.add_argument("genes", nargs="+", help="RSEM .genes.results files")
    args = ap.parse_args()

    rows = []
    for path in args.genes:
        sample = os.path.basename(path).replace(".genes.results", "")
        rows.append((sample, count_detected(path, args.min_count)))

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "genes_detected"])
        writer.writerows(rows)

    print(f"[genes_detected] wrote {len(rows)} samples to {args.out}")


if __name__ == "__main__":
    main()
