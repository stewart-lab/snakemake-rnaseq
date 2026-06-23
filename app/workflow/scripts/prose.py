import argparse
import csv
import os
import re

_tool_to_version = dict()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--versions-dir", required=True,
                   help="directory of per-env version files (key=value lines)")
    p.add_argument("--exclude", default="",
                   help="exclude CSV (sample, stage, reason); may be empty")
    p.add_argument("--min-count", required=True)
    p.add_argument("--min-samples", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True, help="output methods markdown file")
    args = p.parse_args()

    parse_tool_versions(args.versions_dir)
    dropped_samples_by_reason = read_exclusions_csv(args.exclude)
    paragraph = render_paragraph(dropped_samples_by_reason, args.min_count, args.min_samples, args.fasta, args.gtf)
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
            line = f.readline().strip()

        # Version files come in two shapes: "key=value" (from R packageVersion)
        # and a tool's own --version free text (e.g. "FastQC v0.12.1"). For
        # key=value, drop the key first since the key itself can contain digits
        # (e.g. "deseq2"); then take the first digit and everything after it.
        value = line.split("=", 1)[1] if "=" in line else line
        m = re.search(r"\d.*", value)
        version_number = m.group() if m else None

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

def render_paragraph(dropped_samples_by_reason: dict, min_count: int, min_samples: int, fasta: str, gtf: str):
    sentences = [
        _fastp(),
        _fastqc(), 
        _alignment_reference(fasta, gtf),
        _rsem(),
        # TODO: Picard, etc goes here.
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
    return f"Read quality was assessed before and after trimming with FastQC (v{version('fastqc')})."

def _alignment_reference(fasta: str, gtf: str):
    fasta = os.path.basename(fasta)
    gtf = os.path.basename(gtf)
    return f"The alignment reference was prepared with RSEM using the [[FASTA goes here, e.g. GRCh38 DNA primary assembly; raw filename was {fasta}]] FASTA and the [[GTF goes here, e.g. Ensembl GRCh38.111; raw filename was {gtf}]] GTF annotation."

def _rsem():
    return f"Trimmed reads were aligned to the reference and quantified at the gene level using RSEM (v{version('rsem')}) with STAR (v{version('star')})."

def _multiqc():
    return f"Quality-control metrics from fastp, FastQC, and RSEM were aggregated into a single report with MultiQC (v{version('multiqc')})."

def _tximport():
    return f"Per-sample RSEM gene-level results were imported and assembled into a count matrix with tximport (v{version('tximport')})."

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
    return f"Pairwise differential expression analysis was performed with DESeq2 (v{version('deseq2')})."

def _gene_symbol_mapping():
    return "Gene identifiers were mapped to gene symbols using the gene_name attributes in the GTF annotation."

def _gsva_limma():
    return f"Gene-set enrichment scores were computed for the [[GENE SETS, e.g. KEGG, Hallmark, etc.]] gene sets from [[GENE SETS ORIGIN, e.g. MSigDB]] using GSVA (v{version('gsva')}), and differential enrichment between conditions was tested with limma (v{version('limma')})."

def _gsea():
    return f"In parallel, gene set enrichment analysis (GSEA) was performed for each contrast on the DESeq2 statistic-ranked gene lists over the same gene sets using fgsea (v{version('fgsea')})."

def _fdr_etc():
    return "Multiple testing correction was performed with the Benjamini Hochberg method (built in to DESeq2 and limma). Volcano plots and heatmaps were generated using the Plotly Python library."

def render_references():
    """Bibliographic references for every tool/package used in the pipeline."""
    references = [
        "fastp: Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884-i890.",
        "FastQC: Andrews S. FastQC: A Quality Control Tool for High Throughput Sequence Data. 2010. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
        "MultiQC: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048.",
        "STAR: Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15-21.",
        "RSEM: Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics. 2011;12:323.",
        "tximport: Soneson C, Love MI, Robinson MD. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research. 2015;4:1521.",
        "DESeq2: Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15(12):550.",
        "GSVA: Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. BMC Bioinformatics. 2013;14:7.",
        "limma: Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research. 2015;43(7):e47.",
        "fgsea: Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A. Fast gene set enrichment analysis. bioRxiv. 2021. doi:10.1101/060012.",
        "Plotly: Plotly Technologies Inc. Collaborative data science. Montréal, QC: Plotly Technologies Inc.; 2015. https://plot.ly",
        "Benjamini-Hochberg procedure: Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B. 1995;57(1):289-300.",
    ]
    return "\n".join(f"- {ref}" for ref in references)

if __name__ == "__main__":
    main()
