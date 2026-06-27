import argparse
import csv
import os
import re
from collections import Counter

_tool_to_version = dict()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--versions-dir", required=True,
                   help="directory of per-env version files (key=value lines)")
    p.add_argument("--exclude", default="",
                   help="exclude CSV (sample, reason); may be empty")
    p.add_argument("--min-count", required=True)
    p.add_argument("--min-samples", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--strandedness-samplesheet", required=True,
                   help="samples.csv with the strandedness column populated "
                        "(003_strandedness/samples.csv)")
    p.add_argument("--samples", required=True,
                   help="the original input samples.csv, used to tell which "
                        "strandedness values were inferred vs. user-specified")
    p.add_argument("--out", required=True, help="output methods markdown file")
    args = p.parse_args()

    parse_tool_versions(args.versions_dir)
    dropped_samples_by_reason = read_exclusions_csv(args.exclude)
    final, inferred = read_strandedness(args.strandedness_samplesheet, args.samples)
    strand_clause = strandedness_clause(final, inferred)
    paragraph = render_paragraph(dropped_samples_by_reason, args.min_count, args.min_samples,
                                 args.fasta, args.gtf, strand_clause)
    references = render_references()

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write("## Methods\n\n")
        f.write(paragraph)
        f.write("\n\nReferences:\n\n")
        f.write(references + "\n")

    print(f"Wrote methods paragraph to {args.out}")


def version(tool: str):
    global _tool_to_version
    return _tool_to_version.get(tool, "NA")

def parse_tool_versions(versions_dir: str) -> None:
    """Read tool versions from versions_dir/*.txt into one dict."""
    if not versions_dir or not os.path.isdir(versions_dir):
        return
    
    global _tool_to_version
    for name in sorted(os.listdir(versions_dir)):
        if not name.endswith(".txt"):
            continue
        with open(os.path.join(versions_dir, name)) as f:
            for line in f.readlines():
                line = line.strip()

                # when we reach a line <50 chars long, we assume we have found
                # the line containing the version info. in most cases, this
                # is the first line. but sometimes there are errors/warnings/etc.
                # that we will ignore.
                if len(line) < 50:
                    break

        # Version files come in two shapes: "key=value" (from R packageVersion)
        # and a tool's own --version free text (e.g. "FastQC v0.12.1"). For
        # key=value, drop the key first since the key itself can contain digits
        # (e.g. "deseq2"); then take the first digit and everything after it.
        value = line.split("=", 1)[1] if "=" in line else line
        m = re.search(r"\d.*", value)
        version_number = m.group() if m else None

        if version_number and len(version_number) > 20:
            version_number = "[[** ERROR **]]"

        tool = os.path.basename(name).replace(".txt", "")
        _tool_to_version[tool] = version_number

    print(str(_tool_to_version))

def read_exclusions_csv(path):
    """Return {reason: [sample, ...]} from the exclude table."""
    dropped_samples_by_reason = dict()
    if not path or not os.path.exists(path):
        return dropped_samples_by_reason
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            sample = (row.get("sample") or "").strip()
            reason = (row.get("reason") or "").strip()

            if sample:
                dropped_samples_by_reason.setdefault(reason, []).append(sample)
    return dropped_samples_by_reason

def read_strandedness(populated_path, input_path):
    """Per-sample final strandedness and whether each value was RSeQC-inferred.

    Returns (final, inferred):
      final[sample]    -> 'forward' / 'reverse' / 'unstranded' (handed to RSEM)
      inferred[sample] -> True if RSeQC inferred it (i.e. it was blank in the
                          original input sample sheet)
    """
    final = {}
    with open(populated_path, newline="") as fh:
        for row in csv.DictReader(fh):
            final[row["sample"]] = (row.get("strandedness") or "").strip().lower()

    provided = {}
    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh)
        has_col = bool(reader.fieldnames) and "strandedness" in reader.fieldnames
        for row in reader:
            value = (row.get("strandedness") or "").strip().lower() if has_col else ""
            provided[row["sample"]] = value

    inferred = {s: provided.get(s, "") == "" for s in final}
    return final, inferred


def strandedness_clause(final: dict, inferred: dict):
    """A clause describing the strandedness handed to RSEM, e.g.
    'using reverse strandedness as determined by RSeQC infer_experiment.py'.

    When samples disagree, the per-sample breakdown is summarized instead of
    claiming a single library-wide strandedness.
    """
    values = [v for v in final.values() if v]
    if not values:
        return ""

    any_inferred = any(inferred.get(s) for s in final)
    all_inferred = all(inferred.get(s) for s in final)
    if all_inferred:
        provenance = " as determined by RSeQC infer_experiment.py"
    elif any_inferred:
        provenance = (" as specified in the sample sheet, or determined by RSeQC "
                      "infer_experiment.py where unspecified")
    else:
        provenance = " as specified in the sample sheet"

    distinct = sorted(set(values))
    if len(distinct) == 1:
        value = distinct[0]
        base = ("treating the libraries as unstranded" if value == "unstranded"
                else f"using {value} strandedness")
        return base + provenance

    # Mixed strandedness across samples: summarize the per-sample breakdown.
    counts = Counter(values)
    parts = [
        f"{counts[value]} {'unstranded' if value == 'unstranded' else value + '-stranded'}"
        for value in ("forward", "reverse", "unstranded") if counts.get(value)
    ]
    return f"using per-sample library strandedness ({', '.join(parts)}){provenance}"


def render_paragraph(dropped_samples_by_reason: dict, min_count: int, min_samples: int, fasta: str, gtf: str, strand_clause: str = ""):
    sentences = [
        _fastp(),
        _fastqc(),
        _alignment_reference(fasta, gtf),
        _rsem(strand_clause),
        _picard(),
        _multiqc(),
        _tximport(),
        _drop_samples(dropped_samples_by_reason),
        _low_count_filter(min_count, min_samples),
        _deseq2(),
        _gene_symbol_mapping(),
        _gsva_limma(),
        _gsea(),
        _fdr_etc(),
    ]
    paragraph = " ".join(sentences)
    paragraph = paragraph.replace("  ", " ")
    return paragraph

def _fastp():
    return f"Raw reads were preprocessed with fastp (v{version('fastp')}) to trim adapters, filter low-quality reads, and remove polyX tails, with paired-end adapter auto-detection enabled."

def _fastqc():
    return f"Read quality was assessed with FastQC (v{version('fastqc')})."

def _alignment_reference(fasta: str, gtf: str):
    fasta = os.path.basename(fasta)
    gtf = os.path.basename(gtf)
    return f"The alignment reference was prepared with RSEM using the [[FASTA goes here, e.g. GRCh38 DNA primary assembly; raw filename was {fasta}]] FASTA and the [[GTF goes here, e.g. Ensembl GRCh38.111; raw filename was {gtf}]] GTF annotation."

def _rsem(strand_clause: str = ""):
    sentence = (f"Trimmed reads were aligned to the reference and quantified at the "
                f"gene level using RSEM (v{version('rsem')}) with STAR "
                f"(v{version('star')}) in paired-end mode")
    if strand_clause:
        sentence += ", " + strand_clause
    return sentence + "."

def _picard():
    return f"Post-alignment RNA-seq metrics, including ribosomal, exonic, intronic and intergenic rates and 5'-to-3' transcript coverage bias, were collected with Picard CollectRnaSeqMetrics (v{version('picard')})."

def _multiqc():
    return f"Quality-control metrics from fastp, FastQC, RSeQC, RSEM, and Picard were aggregated into a single report with MultiQC (v{version('multiqc')})."

def _tximport():
    return f"Per-sample RSEM gene-level results were assembled into a count matrix with tximport (v{version('tximport')})."

def _drop_samples(dropped_samples_by_reason: dict):
    numbers = {
        1: "one",
        2: "two",
        3: "three",
        4: "four",
        5: "five",
        6: "six",
        7: "seven",
        8: "eight",
        9: "nine",
        10: "ten",
    }
    
    clauses = []
    for reason, samples in dropped_samples_by_reason.items():
        n_samples = len(samples)
        noun = "samples" if n_samples > 1 else "sample"
        verb = "were" if n_samples > 1 else "was"
        clause = f"{numbers.get(n_samples, 'NA')} {noun} ({', '.join(samples)}) {verb} dropped because {reason}"
        clauses.append(clause)
    
    # join clauses with semicolons
    sentence = "; ".join(clauses) + "." if clauses else ""

    # capitalize first letter in sentence
    if len(sentence) > 0:
        sentence = sentence[0].upper() + sentence[1:]

    return sentence

def _low_count_filter(min_count: int, min_samples: int):
    return f"Genes that did not have at least {min_count} counts in {min_samples} samples were removed."

def _deseq2():
    return f"Pairwise differential expression analysis was performed with DESeq2 (v{version('deseq2')}), fitting a negative-binomial generalized linear model with a ~condition design and testing each contrast with the Wald test."

def _gene_symbol_mapping():
    return "Gene identifiers were mapped to gene symbols using the gene_name attributes in the GTF annotation."

def _gsva_limma():
    return f"Gene-set enrichment scores were computed for the [[GENE SETS, e.g. KEGG, Hallmark, etc.]] gene sets from [[GENE SETS ORIGIN, e.g. MSigDB]] using GSVA (v{version('gsva')}) with a Gaussian kernel and tau = 1 on log2-transformed normalized counts, and differential enrichment between conditions was tested with limma (v{version('limma')})."

def _gsea():
    return f"Gene set enrichment analysis (GSEA) was performed for each contrast on gene lists ranked by the DESeq2 Wald statistic over the same gene sets using fgsea (v{version('fgsea')}), testing gene sets of 15 to 500 genes."

def _fdr_etc():
    return "Multiple testing correction was performed with the Benjamini-Hochberg method (built in to DESeq2 and limma). [[Plotly can be referenced here if you used it for figure generation.]]"

def render_references():
    """Bibliographic references for every tool/package used in the pipeline."""
    references = [
        "fastp: https://doi.org/10.1093/bioinformatics/bty560",
        "FastQC: Andrews S. FastQC: A Quality Control Tool for High Throughput Sequence Data. 2010. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
        "RSeQC: https://doi.org/10.1093/bioinformatics/bts356",
        "MultiQC: https://doi.org/10.1093/bioinformatics/btw354",
        "STAR: https://doi.org/10.1093/bioinformatics/bts635",
        "RSEM: https://doi.org/10.1186/1471-2105-12-323",
        "Picard: Broad Institute. Picard Toolkit. 2019. https://broadinstitute.github.io/picard/",
        "tximport: https://doi.org/10.12688/f1000research.7563.2",
        "DESeq2: https://doi.org/10.1186/s13059-014-0550-8",
        "GSVA: https://doi.org/10.1186/1471-2105-14-7",
        "limma: https://doi.org/10.1093/nar/gkv007",
        "fgsea: https://doi.org/10.1101/060012",
        "Plotly: Plotly Technologies Inc. Collaborative data science. Montréal, QC: Plotly Technologies Inc.; 2015. https://plot.ly",
        "Benjamini-Hochberg procedure: https://doi.org/10.1111/j.2517-6161.1995.tb02031.x",
    ]
    return "\n".join(f"- {ref}" for ref in references)

if __name__ == "__main__":
    main()
