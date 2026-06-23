#!/usr/bin/env python3
"""Extract ribosomal-RNA features from a GTF into a BED.

`picard BedToIntervalList` turns this BED plus the sequence dictionary into the
rRNA interval_list that CollectRnaSeqMetrics needs for RIBOSOMAL_INTERVALS, so
this script only has to do the GTF parsing (the part no standard tool covers).
Only rRNA gene/transcript/exon features are emitted; if the GTF has none, an
empty BED is written and the resulting interval_list is header-only (Picard
then reports 0 ribosomal bases).
"""
import argparse
import os
import re

# Biotype attribute values that mark a feature as ribosomal RNA. Sources name
# the biotype attribute differently, so several keys are checked.
RRNA_BIOTYPES = {"rRNA", "Mt_rRNA", "rRNA_pseudogene"}
BIOTYPE_KEYS = ("gene_biotype", "gene_type", "transcript_biotype", "transcript_type")


def get_attr(attributes, key):
    """Value of a GTF attribute (key "value";) or None."""
    m = re.search(key + r'\s+"([^"]*)"', attributes)
    return m.group(1) if m else None


def is_rrna(attributes):
    for key in BIOTYPE_KEYS:
        value = get_attr(attributes, key)
        if value and (value in RRNA_BIOTYPES or value.endswith("rRNA")):
            return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True, help="output BED of rRNA features")
    args = ap.parse_args()

    seen = set()
    n = 0
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.gtf) as gtf, open(args.out, "w") as out:
        for line in gtf:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, _, feattype, start, end, _, strand, _, attributes = fields[:9]
            if feattype not in ("gene", "transcript", "exon"):
                continue
            if not is_rrna(attributes):
                continue
            start0 = int(start) - 1  # GTF is 1-based inclusive; BED is 0-based.
            key = (seqid, start0, end, strand)
            if key in seen:
                continue
            seen.add(key)
            name = get_attr(attributes, "gene_id") or get_attr(attributes, "gene_name") or "rRNA"
            strand = strand if strand in ("+", "-") else "+"
            out.write(f"{seqid}\t{start0}\t{end}\t{name}\t0\t{strand}\n")
            n += 1

    print(f"[rrna_bed] wrote {n} rRNA features to {args.out}")


if __name__ == "__main__":
    main()
