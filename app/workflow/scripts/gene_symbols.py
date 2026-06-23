#!/usr/bin/env python3
"""Prepend gene_symbol/gene_id columns to a DESeq2 output table.

Maps the gene IDs (the first column of the input CSV) to gene names using a
gene_id,gene_name CSV built from the GTF (see gtf_gene_map.py). This replaces
the organism-specific org.Hs.eg.db lookup so the pipeline works for any species.
Genes with no name fall back to their gene ID so every row keeps a usable label.

Uses only the standard library so it can run in the base environment without a
dedicated conda env.
"""
import argparse
import csv


def load_gene_map(path):
    """Return {gene_id: gene_name} from a gene_id,gene_name CSV."""
    mapping = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            gene_id = (row.get("gene_id") or "").strip()
            gene_name = (row.get("gene_name") or "").strip()
            if gene_id:
                mapping[gene_id] = gene_name
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True,
                    help="DESeq2 CSV keyed by gene ID (first column)")
    ap.add_argument("--map", required=True,
                    help="gene_id,gene_name CSV built from the GTF")
    ap.add_argument("--output", required=True,
                    help="output CSV with gene_symbol + gene_id prepended")
    args = ap.parse_args()

    gene_map = load_gene_map(args.map)

    n_total = 0
    n_mapped = 0
    with open(args.input, newline="", encoding="utf-8") as fin, \
            open(args.output, "w", newline="", encoding="utf-8") as fout:
        reader = csv.reader(fin)
        writer = csv.writer(fout)

        header = next(reader)
        # The first column is the gene ID index; the rest are data columns.
        writer.writerow(["gene_symbol", "gene_id"] + header[1:])

        for row in reader:
            if not row:
                continue
            gene_id = row[0]
            symbol = gene_map.get(gene_id, "")
            if symbol:
                n_mapped += 1
            else:
                symbol = gene_id  # fall back to the gene ID when unnamed
            n_total += 1
            writer.writerow([symbol, gene_id] + row[1:])

    print(f"[gene_symbols] mapped {n_mapped}/{n_total} gene IDs to names "
          f"-> {args.output}")


if __name__ == "__main__":
    main()
