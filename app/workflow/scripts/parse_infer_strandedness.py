#!/usr/bin/env python3
"""Turn RSeQC `infer_experiment.py` output into a strandedness call.

`infer_experiment.py` reports, for a subsample of aligned reads, the fraction
explained by each read-orientation pattern. This script reads that report and
writes the inferred call -- one of forward / reverse / unstranded -- as a single
line. The RSEM and Picard rules read that word (or the explicit value from the
samples.csv strandedness column) and translate it to their own option.

Orientation patterns (from infer_experiment.py):
  paired-end  forward:  "1++,1--,2+-,2-+"   reverse:  "1+-,1-+,2++,2--"
  single-end  forward:  "++,--"             reverse:  "+-,-+"

"forward" means the read (or read 1) is on the same strand as the transcript
(fr-secondstrand / ISF); RSEM --strandedness forward, Picard
SECOND_READ_TRANSCRIPTION_STRAND. "reverse" is the dUTP/Illumina-stranded case
(fr-firststrand / ISR), the most common stranded protocol; RSEM --strandedness
reverse, Picard FIRST_READ_TRANSCRIPTION_STRAND. "unstranded" maps to RSEM
--strandedness none and Picard NONE.
"""
import argparse
import os
import re

# (label, regex capturing the trailing fraction) for the forward and reverse
# orientation lines, covering both paired-end and single-end report wording.
FORWARD_PATTERNS = (
    r'"1\+\+,1--,2\+-,2-\+":\s*([0-9.]+)',  # paired-end
    r'"\+\+,--":\s*([0-9.]+)',              # single-end
)
REVERSE_PATTERNS = (
    r'"1\+-,1-\+,2\+\+,2--":\s*([0-9.]+)',  # paired-end
    r'"\+-,-\+":\s*([0-9.]+)',              # single-end
)

def _first_fraction(text, patterns):
    """First trailing fraction matched by any of `patterns`, or None."""
    for pattern in patterns:
        m = re.search(pattern, text)
        if m:
            return float(m.group(1))
    return None


def call_strandedness(forward, reverse, threshold):
    """Classify the library from the forward/reverse fractions.

    The decision is made over the *explained* reads only (undetermined reads are
    ignored): if one orientation accounts for at least `threshold` of the
    explained reads, the library is stranded that way; otherwise unstranded.
    """
    forward = forward or 0.0
    reverse = reverse or 0.0
    explained = forward + reverse
    if explained <= 0:
        return "unstranded", 0.0
    forward_frac = forward / explained
    if forward_frac >= threshold:
        return "forward", forward_frac
    if forward_frac <= 1 - threshold:
        return "reverse", forward_frac
    return "unstranded", forward_frac


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--infer-experiment", required=True,
                    help="infer_experiment.py output text file")
    ap.add_argument("--threshold", type=float, default=0.8,
                    help="min fraction of explained reads in one orientation to "
                         "call the library stranded (default: 0.8)")
    ap.add_argument("--out", required=True,
                    help="output file; the inferred call (forward/reverse/"
                         "unstranded) is written as a single line")
    args = ap.parse_args()

    with open(args.infer_experiment) as fh:
        text = fh.read()

    forward = _first_fraction(text, FORWARD_PATTERNS)
    reverse = _first_fraction(text, REVERSE_PATTERNS)
    if forward is None and reverse is None:
        raise SystemExit(
            f"[strandedness] could not find orientation fractions in "
            f"{args.infer_experiment}; is this infer_experiment.py output?"
        )

    call, forward_frac = call_strandedness(forward, reverse, args.threshold)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        fh.write(call + "\n")

    print(f"[strandedness] call={call} "
          f"(forward={forward or 0.0:.4f}, reverse={reverse or 0.0:.4f}, "
          f"threshold={args.threshold})")


if __name__ == "__main__":
    main()
