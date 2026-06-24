#!/usr/bin/env python3
"""Low-count gene filter for the all-samples and samples-dropped count matrices.

A gene is kept if it has at least ``--min-count`` counts in at least
``--min-samples`` samples. The same threshold is applied independently to the
all-samples and samples-dropped matrices. The design tables and contrasts table
are passed through unchanged so every downstream step reads from one directory.
"""
import argparse
import shutil

import pandas as pd


def filter_genes(counts, min_count, min_samples):
    """Return the count matrix keeping only genes that clear the threshold."""
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts[keep]


def filter_file(counts_path, out_path, min_count, min_samples):
    counts = pd.read_csv(counts_path, index_col="feature_id")
    filter_genes(counts, min_count, min_samples).to_csv(out_path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--counts-all", required=True)
    ap.add_argument("--counts-dropped", required=True)
    ap.add_argument("--design-all", required=True)
    ap.add_argument("--design-dropped", required=True)
    ap.add_argument("--contrasts", required=True)
    ap.add_argument("--out-counts-all", required=True)
    ap.add_argument("--out-counts-dropped", required=True)
    ap.add_argument("--out-design-all", required=True)
    ap.add_argument("--out-design-dropped", required=True)
    ap.add_argument("--out-contrasts", required=True)
    ap.add_argument("--min-count", type=float, required=True)
    ap.add_argument("--min-samples", type=int, required=True)
    args = ap.parse_args()

    filter_file(args.counts_all, args.out_counts_all, args.min_count, args.min_samples)
    filter_file(args.counts_dropped, args.out_counts_dropped, args.min_count, args.min_samples)

    # Pass the design tables and contrasts through unchanged.
    pd.read_csv(args.design_all, index_col=0).to_csv(args.out_design_all)
    pd.read_csv(args.design_dropped, index_col=0).to_csv(args.out_design_dropped)
    shutil.copy(args.contrasts, args.out_contrasts)
    print(f"[filter_low_count_genes] kept genes with >= {args.min_count} counts "
          f"in >= {args.min_samples} samples")


if __name__ == "__main__":
    main()
