#!/usr/bin/env python3
"""Build a gene_id -> gene_name mapping from a GTF/GFF annotation.

Parses the annotation with gffutils, which handles both the GTF (key "value")
and GFF3 (key=value) attribute conventions, and writes a two-column CSV
(gene_id,gene_name). This replaces organism-specific annotation packages (e.g.
org.Hs.eg.db) so the pipeline works for any species whose annotation provides
gene names.

Different sources name the gene-name attribute differently, so several candidate
keys are tried in priority order (see NAME_KEYS). Gene IDs are taken verbatim
from the annotation, so they match the IDs RSEM emits (RSEM is built from the
same file). When a gene has no name attribute, its gene_name is set equal to its
gene_id so downstream steps always have a label.
"""
import argparse
import csv
import gzip

from gffutils.feature import feature_from_line

# Attribute keys that may hold a human-readable gene name, in priority order:
#   gene_name   - Ensembl, GENCODE
#   gene        - NCBI RefSeq
#   gene_symbol - some other annotation pipelines
#   Name        - GFF3-style exports
NAME_KEYS = ("gene_name", "gene", "gene_symbol", "Name")


def _open(path):
    """Open a plain or gzipped text file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def first_name(attributes):
    """First non-empty value among the candidate name keys, or None.

    gffutils stores each attribute as a list of values, so we take the first.
    """
    for key in NAME_KEYS:
        values = attributes.get(key)
        if values and values[0]:
            return values[0]
    return None


def build_map(gtf_path):
    """Return an ordered dict gene_id -> gene_name parsed from the annotation.

    A real gene name always wins over the gene_id fallback, regardless of which
    line it is first seen on.
    """
    mapping = {}
    with _open(gtf_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if line.count("\t") < 8:  # not a complete 9-column feature line
                continue

            feature = feature_from_line(line)
            gene_id = (feature.attributes.get("gene_id") or [None])[0]
            if not gene_id:
                continue

            name = first_name(feature.attributes)
            if gene_id not in mapping:
                mapping[gene_id] = name or gene_id
            elif name and mapping[gene_id] == gene_id:
                # We only had the fallback so far; upgrade to the real name.
                mapping[gene_id] = name
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True,
                    help="GTF/GFF annotation (optionally gzipped)")
    ap.add_argument("--out", required=True,
                    help="output CSV with columns gene_id,gene_name")
    args = ap.parse_args()

    mapping = build_map(args.gtf)

    with open(args.out, "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(["gene_id", "gene_name"])
        for gene_id, gene_name in mapping.items():
            writer.writerow([gene_id, gene_name])

    named = sum(1 for gid, name in mapping.items() if name != gid)
    print(f"[gtf_gene_map] wrote {len(mapping)} genes "
          f"({named} with a gene name) to {args.out}")


if __name__ == "__main__":
    main()
