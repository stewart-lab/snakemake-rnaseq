#!/usr/bin/env python3
"""Drop excluded samples from the raw count matrix and the design table.

The samples to drop are passed as positional-style ``--drop`` arguments (zero or
more sample names). Both the count matrix columns and the design rows for those
samples are removed; a missing sample is an error so typos in the exclude table
fail loudly rather than silently dropping nothing.
"""
import argparse

import pandas as pd


def drop_samples(counts, design, drop):
    """Return (counts, design) with the named samples removed from both."""
    counts = counts.drop(columns=drop, errors="raise")
    design = design.drop(index=drop, errors="raise")
    return counts, design


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--counts", required=True, help="raw count matrix CSV (feature_id index)")
    ap.add_argument("--design", required=True, help="design CSV (sample index)")
    ap.add_argument("--out-counts", required=True, help="output count matrix CSV")
    ap.add_argument("--out-design", required=True, help="output design CSV")
    ap.add_argument("--drop", nargs="*", default=[],
                    help="sample names to drop (may be empty)")
    args = ap.parse_args()

    counts = pd.read_csv(args.counts, index_col="feature_id")
    design = pd.read_csv(args.design, index_col=0)

    counts, design = drop_samples(counts, design, args.drop)

    counts.to_csv(args.out_counts)
    design.to_csv(args.out_design)
    print(f"[drop_samples] dropped {len(args.drop)} sample(s): "
          f"{', '.join(args.drop) or '(none)'}")


if __name__ == "__main__":
    main()
