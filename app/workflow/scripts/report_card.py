#!/usr/bin/env python3
"""Render a text-based QC 'report card' for one sample.

Metrics are collected from the RSEM/STAR + Picard QC outputs:
  - uniquely_mapped_pct, total_reads_M   <- STAR log
  - genes_detected                       <- precomputed genes-detected table
  - rrna_pct, exonic_pct, intergenic_pct,
    median_5p_3p_bias                     <- Picard CollectRnaSeqMetrics

Any metric that cannot be read is left as None and shown as "n/a" (graded NA,
which is neutral: it does not trigger a drop and is not flagged).
"""
import argparse
import csv
import os
from math import log2

# --- four-tier grade scale ---
# Tiers, worst to best: FAIL < WARN < OK < GOOD.
# FAIL is the empty-bar state (value below the bottom edge); the three filled
# sections of the bar are WARN | OK | GOOD. A sample is dropped on any FAIL.
GOOD, OK, WARN, FAIL, NA = "GOOD", "OK", "WARN", "FAIL", "NA"
GRADE_MARK = {
    GOOD: "[good]",
    OK:   "[ok]  ",
    WARN: "[warn]",
    FAIL: "[fail]",
    NA:   "[n/a] ",
}

FULL = '█'
LIGHT = '░'
SECTION = 10                 # chars per section
BAR_WIDTH = SECTION * 3      # warn | ok | good


def _bias_distance(v):
    """Symmetric distance of a 5'/3' ratio from the ideal 1.0: abs(log2(v)).

    A ratio and its reciprocal map to the same distance, so over- and
    under-coverage of the 5' end are penalized equally. A non-positive or
    missing ratio is treated as maximally biased.
    """
    if v is None or v <= 0:
        return float("inf")
    return abs(log2(v))

# Each metric defines its zone edges as 4 points along the bar, ordered
# from the "worst" end to the "best" end:
#   edges = (fail_edge, warn|ok, ok|good, best)
# Direction is implied by ordering: "higher is better" is ascending
# (fail_edge < best); "lower is better" is descending (fail_edge > best).
#
# The fail_edge is the bottom of the bar. A value past it (worse than it) is
# FAIL and draws an empty bar. The three bands above it map 1:1 to the three
# bar sections, so the section where the fill ends is the metric's tier.
#
# grade_fn (last field) optionally transforms the raw value into the quantity
# that is actually graded and drawn. The displayed number column still shows
# the raw value; only grading/bar use the transform. It exists for metrics that
# are not monotonic in the raw value (see 5'/3' bias below). None = identity.
#
# Note: the genes_detected edges are placeholders; they are overridden per run
# from the dataset (see main()).
# key, label, unit, fmt, edges, critical_fn, grade_fn
METRICS = [
    ("uniquely_mapped_pct", "Uniquely mapped", "%", "{:>6.0f}", (70, 80, 90, 100),    lambda v: v < 70, None),
    ("total_reads_M",       "Total reads",     "M", "{:>6.1f}", (5, 10, 20, 30),       lambda v: v < 5,  None),
    ("genes_detected",      "Genes detected",  "k", "{:>6.1f}", (0, 0, 0, 0),          lambda v: False,  None),
    ("rrna_pct",            "rRNA",            "%", "{:>6.0f}", (15, 10, 5, 0),         lambda v: v > 30, None),
    ("exonic_pct",          "Exonic",          "%", "{:>6.0f}", (55, 70, 85, 100),     lambda v: v < 55, None),
    ("intergenic_pct",      "Intergenic",      "%", "{:>6.0f}", (15, 10, 5, 0),        lambda v: v > 30, None),
    # 5'/3' bias is a ratio: a value and its reciprocal are equally bad (a 5'
    # enrichment of 2.0 == a 3' enrichment of 0.5). It is graded on the symmetric
    # log2 distance from the ideal 1.0, so both over- and under-coverage of the 5'
    # end are penalized. Edges are in that distance space (lower = better):
    #   >1.32 fail (raw <0.40 or >2.50), >1.00 warn (<0.50 or >2.00),
    #   >0.30 ok   (<0.81 or >1.23),     else good (0.81-1.23).
    ("median_5p_3p_bias",   "5'/3' bias",      " ", "{:>6.2f}", (1.32, 1.0, 0.3, 0.0), lambda v: _bias_distance(v) > 1.32, _bias_distance),
]

# Tier -> bar section index (0 = leftmost filled section). FAIL has no section.
SECTION_BY_TIER = {WARN: 0, OK: 1, GOOD: 2}


def grade_from_edges(value, edges):
    """Return the tier for `value`. edges = (fail_edge, warn|ok, ok|good, best),
    ordered worst->best (ascending or descending)."""
    warn_threshold, ok_threshold, good_threshold, best = edges
    ascending = best > warn_threshold
    if ascending:
        if value < warn_threshold:   return FAIL
        if value < ok_threshold:     return WARN
        if value < good_threshold:   return OK
        return GOOD
    else:
        if value > warn_threshold:   return FAIL
        if value > ok_threshold:     return WARN
        if value > good_threshold:   return OK
        return GOOD


def _interp_section(value, lo, hi, section=SECTION):
    """Fractional fill (0..section) for value within band [lo, hi], any order."""
    if hi == lo:
        return section
    frac = (value - lo) / (hi - lo)
    frac = max(0.0, min(1.0, frac))
    return frac * section


def _empty_bar():
    return ' '.join((LIGHT * SECTION) for _ in range(3))


def make_zone_bar(value, edges, grade):
    """Bar of 3 sections (warn | ok | good). FAIL draws an empty bar; the other
    tiers fill left->right, ending inside their own section."""
    if grade == FAIL:
        return _empty_bar()

    fail_edge, warn_ok, ok_good, best = edges
    band = {WARN: (fail_edge, warn_ok), OK: (warn_ok, ok_good), GOOD: (ok_good, best)}[grade]
    section_index = SECTION_BY_TIER[grade]

    filled = section_index * SECTION + _interp_section(value, band[0], band[1])
    n = max(0, min(BAR_WIDTH, int(round(filled))))

    bar = FULL * n + LIGHT * (BAR_WIDTH - n)
    # space the sections apart: warn | ok | good
    return ' '.join(bar[i:i + SECTION] for i in range(0, BAR_WIDTH, SECTION))


def report_card(sample_id, m, edges_override=None):
    edges_override = edges_override or {}
    lines = []
    lines.append("=" * 78)
    lines.append("SAMPLE QC REPORT CARD: " + sample_id)
    lines.append("                                   [  warn   |   ok    |   good  ]")
    lines.append("=" * 78)

    fail_reasons = []

    for key, label, unit, fmt, edges, crit_fn, grade_fn in METRICS:
        edges = edges_override.get(key, edges)
        val = m.get(key)
        if val is None:
            grade, crit = NA, False
            shown = "   n/a"
            bar = _empty_bar()
        else:
            gval = grade_fn(val) if grade_fn else val
            grade = grade_from_edges(gval, edges)
            crit = bool(crit_fn(val))
            if crit:
                grade = FAIL
            shown = fmt.format(val)
            bar = make_zone_bar(gval, edges, grade)

        if grade == FAIL:
            tag = " [critical]" if crit else ""
            fail_reasons.append(f"{label} = {fmt.strip().format(val).strip()}{unit}{tag}")

        crit_mark = " !" if crit else "  "
        lines.append(
            f" {GRADE_MARK[grade]}{crit_mark} {label:<16}{shown}{unit:<2} {bar}"
        )

    # Recommendation: drop on any FAIL, otherwise keep.
    verdict = "drop" if fail_reasons else "keep"

    lines.append("-" * 78)
    lines.append(f"  recommendation: {verdict}")
    if fail_reasons:
        lines.append("  failed metrics:")
        for r in fail_reasons:
            lines.append("    - " + r)
    lines.append("=" * 78)
    lines.append("\n\n")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# metric collection
# --------------------------------------------------------------------------- #
def _float(value, ndigits=None):
    try:
        x = float(value)
        return round(x, ndigits) if ndigits is not None else x
    except (TypeError, ValueError):
        return None


def _pct(value):
    """Picard fraction (0-1, possibly blank) -> rounded percent, or None."""
    x = _float(value)
    return round(x * 100) if x is not None else None


def parse_star_log(path):
    """uniquely_mapped_pct and total_reads_M from a STAR Log"""
    values = {}
    with open(path) as fh:
        for line in fh:
            if "|" not in line:
                continue
            label, value = line.split("|", 1)
            values[label.strip()] = value.strip()

    out = {}
    if values.get("Number of input reads", "").isdigit():
        out["total_reads_M"] = int(values["Number of input reads"]) / 1e6
    pct = _float(values.get("Uniquely mapped reads %", "").rstrip("%"))
    if pct is not None:
        out["uniquely_mapped_pct"] = round(pct)
    return out


def _picard_metrics_row(path):
    """First METRICS data row of a Picard metrics file as a {column: value} dict."""
    with open(path) as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if line.startswith("## METRICS CLASS"):
            header = lines[i + 1].rstrip("\n").split("\t")
            data = lines[i + 2].rstrip("\n").split("\t")
            return dict(zip(header, data))
    return {}


def parse_rnaseq_metrics(path):
    """rRNA/exonic/intergenic % and 5'->3' bias from Picard CollectRnaSeqMetrics."""
    m = _picard_metrics_row(path)
    out = {
        "rrna_pct": _pct(m.get("PCT_RIBOSOMAL_BASES")),
        "exonic_pct": _pct(m.get("PCT_MRNA_BASES")),
        "intergenic_pct": _pct(m.get("PCT_INTERGENIC_BASES")),
    }
    bias = _float(m.get("MEDIAN_5PRIME_TO_3PRIME_BIAS"), ndigits=2)
    if bias is not None:
        out["median_5p_3p_bias"] = bias
    return out


def read_genes_table(path):
    """{sample: genes_detected} from the precomputed genes-detected table."""
    table = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            try:
                table[row["sample"]] = float(row["genes_detected"]) / 1000
            except (KeyError, ValueError):
                continue
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--star-log", required=True, help="STAR log")
    ap.add_argument("--rnaseq-metrics", required=True, help="Picard CollectRnaSeqMetrics output")
    ap.add_argument("--genes-table", required=True, help="genes_detected.csv (sample,genes_detected) for the whole dataset")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    table = read_genes_table(args.genes_table)

    metrics = {}
    metrics.update(parse_star_log(args.star_log))
    if args.sample in table:
        metrics["genes_detected"] = table[args.sample]
    metrics.update(parse_rnaseq_metrics(args.rnaseq_metrics))

    # Dynamic genes-detected target: the most any sample detects, or 10k, whichever
    # is higher. Each sample is then graded relative to that dataset-wide best.
    reference = max(max(table.values(), default=0), 10)
    edges_override = {
        "genes_detected": (
            round(0.75 * reference),   # fail_edge (bottom of bar)
            round(0.85 * reference),   # warn|ok
            round(0.95 * reference),   # ok|good
            reference,                 # best
        ),
    }

    card = report_card(args.sample, metrics, edges_override)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(card)
    print(card)


if __name__ == "__main__":
    main()