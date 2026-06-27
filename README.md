[![](https://img.shields.io/appveyor/build/robertmillikin/snakemake-rnaseq/main)](https://ci.appveyor.com/project/robertmillikin/snakemake-rnaseq/history)
[![](https://img.shields.io/docker/v/stewartlab/snakemake-rnaseq?sort=semver&color=blue)](https://hub.docker.com/r/stewartlab/snakemake-rnaseq/tags?page=1&ordering=last_updated)

# snakemake-rna-seq

A [Snakemake](https://snakemake.github.io/) pipeline for bulk paired-end RNA-seq:
from raw FASTQs to differential expression, gene-set enrichment, QC reports, and
an auto-generated methods paragraph. Each step runs in its own conda environment.

## What it does

Trim → align/quantify → QC → build count matrix → filter → explore (PCA) →
differential expression → gene-set analysis (GSVA **and** GSEA) → reporting.

| Step | Output folder | Tool | Description |
|------|---------------|------|-------------|
| 0 | `000_meta` | graphviz | Information about the pipeline |
| 1 | `001_fastp_trimmed` | fastp | Adapter/quality trimming + polyX removal |
| 2 | `002_fastqc` | FastQC | QC of the trimmed reads |
| – | `<genome_dir>/star_reference` | RSEM + STAR | Build the alignment/quantification index |
| – | `<genome_dir>/picard_reference` | Picard | Build the sequence dictionary, refFlat, and rRNA intervals |
| 3 | `003_strandedness` | RSeQC | Strandedness inference (or provided by `samples.csv`) |
| 4 | `004_rsem` | RSEM + STAR | Align and quantify gene-level expression |
| 5 | `005_qc_reports` | Picard + RSEM/STAR | Per-sample QC "report card" (mapping %, reads, genes detected, rRNA, exonic, intergenic, 5'/3' bias) |
| 6 | `006_multiqc` | MultiQC | Aggregate fastp/FastQC/RSeQC/RSEM/Picard QC into one report |
| 7 | `007_raw_counts` | tximport | Assemble the per-sample results into a count matrix |
| 8 | `008_raw_counts_samples_dropped` | pandas | Drop excluded samples (see `exclude.csv`) |
| 9 | `009_raw_counts_filtered` | pandas | Low-count gene filter (all-samples **and** dropped matrices) |
| 10 | `010_pca` | DESeq2/Plotly | VST → PCA + sample clustermap, all samples |
| 11 | `011_pca_samples_dropped` | DESeq2/Plotly | VST → PCA + sample clustermap, samples dropped |
| 12 | `012_deseq2` | DESeq2 | Differential expression per contrast |
| 13 | `013_deseq2_gene_symbols` | gffutils | Map gene IDs → gene names using the GTF |
| 14 | `014_gsva` | GSVA + limma | Per-sample gene-set scores + differential enrichment |
| 15 | `015_gsea` | fgsea | Ranked-list gene-set enrichment per contrast |
| 16 | `016_prose` | – | Auto-generated methods paragraph |
| – | `reports` | Quarto | Render any user-supplied `.qmd` reports |

## Running

The pipeline is packaged as a Docker image. Running the container will 1) create a workflow directory, 
2) copy the script and environment files into the workflow directory, and then 3) run the Snakemake pipeline.

Snakemake will create the conda environments needed for each step in the pipeline.

```bash
# To run the pipeline, pass your config YAML to the container.
# Your inputs (.fastq files, etc.) must be in a volume mounted to the container.

# Tip: give the mounted directory the same path it has on your filesystem, so
# the paths in your config files are findable when viewing logs, etc.

docker run --rm \
  -v /home/myname/projects:/home/myname/projects \
  snakemake-rnaseq:1.0.0 \
  /home/myname/projects/myproject/config/config.yaml
```

### Suggested project directory structure

The following project structure is suggested:

```text
myproject/
├── config/                 # pipeline config
│   ├── config.yaml         #   ← pass THIS file to the container on startup
│   ├── samples.csv
│   ├── design.csv
│   ├── contrasts.csv
│   ├── gene_sets.csv
│   └── exclude.csv
├── resources/              # pipeline inputs
│   ├── reads/              #   raw FASTQs referenced by samples.csv
│   ├── reference/          #   genome_dir: exactly one FASTA + one GTF
│   └── gene_sets/          #   .gmt files referenced by gene_sets.csv
├── reports/                # reports_dir: your Quarto reports, one folder each
│   └── summary/summary.qmd
├── workflow/               # workflow_dir: pipeline code (Snakefile + envs/rules/scripts)
├── results/                # results_dir: pipeline outputs
└── logs/                   # logs_dir: pipeline logs
```

You are free to name and arrange these however you like — the pipeline only
cares about the paths you set in `config.yaml`. The Docker container creates 
the workflow directory at the path you set in `workflow_dir` in the config 
YAML if it doesn't exist yet.

## Inputs

All inputs specified in the config YAML file (see
[`example/config.yaml`](example/config.yaml)). There are also example CSV files in the 
[`example/`](example/) directory.

### `config.yaml`

| Key | Meaning |
|-----|---------|
| `samples` | Path to the sample sheet CSV |
| `design` | Path to the design CSV |
| `contrasts` | Path to the contrasts CSV |
| `gene_sets` | Path to the gene-sets CSV (lists GMT files) |
| `exclude` | Path to the exclusions CSV (optional) |
| `genome_dir` | Path to the directory containing exactly one FASTA and one GTF |
| `filter.min_count` / `filter.min_samples` | Low-count gene filter thresholds |
| `cores` | Number of CPU cores for the pipeline |
| `workflow_dir` | Path to the directory to write pipeline code (Snakefile + envs/rules/scripts) |
| `results_dir` | Path to the directory to write output of each step |
| `reports_dir` | Path to the directory of user-supplied Quarto `.qmd` reports to render (optional; see [Reports](#reports--quarto-optional)) |
| `logs_dir` | Path to the directory to write logs for each step |

Paths may be absolute, or relative to the config directory.
You may wish to use absolute paths so your colleagues can easily find the data
if you send them your config files.

### Sample sheet — [`example/samples.csv`](example/samples.csv)

One row per sample. FASTQs may be `.fastq` or `.fastq.gz`.

```csv
sample,r1,r2
WT_rep1,../reads/WT_rep1_R1.fastq.gz,../reads/WT_rep1_R2.fastq.gz
KO_rep1,../reads/KO_rep1_R1.fastq.gz,../reads/KO_rep1_R2.fastq.gz
```

An optional `strandedness` column lets you set a sample's library strandedness
explicitly; valid values are `forward`, `reverse`, and `unstranded`. Any sample
left blank — or the whole
column omitted — has its strandedness **inferred** with RSeQC
`infer_experiment.py` (step 3), which aligns a subsample of the reads and reads
out the orientation. The strandedness (given or inferred) is passed to RSEM
(`--strandedness`) and Picard (`STRAND_SPECIFICITY`), and a copy of the sample
sheet with the column fully populated is written to
`results/003_strandedness/samples.csv`.

```csv
sample,r1,r2,strandedness
WT_rep1,../reads/WT_rep1_R1.fastq.gz,../reads/WT_rep1_R2.fastq.gz,reverse
KO_rep1,../reads/KO_rep1_R1.fastq.gz,../reads/KO_rep1_R2.fastq.gz,
```

### Design — [`example/design.csv`](example/design.csv)

Maps each sample to a condition (used as the DESeq2 design `~ condition`).

```csv
sample,condition
WT_rep1,WT
KO_rep1,KO
```

### Contrasts — [`example/contrasts.csv`](example/contrasts.csv)

Each row is one comparison; results are named `<treatment>_vs_<control>`.

```csv
control,treatment
WT,KO
```

### Gene sets — [`example/gene_sets.csv`](example/gene_sets.csv)

A single `gmt` column listing GMT files. All listed GMTs are merged for GSVA and
GSEA. Gene sets must use **gene symbols** (matching the GTF-derived names).
Download collections such as the [MSigDB Hallmark sets](https://www.gsea-msigdb.org/gsea/msigdb/)
as `.gmt` files. A tiny illustrative GMT is in
[`example/gene_sets/example_sets.gmt`](example/gene_sets/example_sets.gmt).

```csv
gmt
example/gene_sets/example_sets.gmt
```

### Exclusions — [`example/exclude.csv`](example/exclude.csv) (optional)

Samples to drop before differential expression (e.g. PCA outliers). Every listed
sample must exist in the sample sheet. The `reason` is
included in the methods paragraph. It's typical to run the pipeline once with
a blank `exclude.csv`, inspect the QC reports, then drop poor quality samples
by adding them to `exclude.csv` and re-running the pipeline.

```csv
sample,reason
KO_rep3,it had a low number of reads
```

### Reference (`genome_dir`)

Point `genome_dir` at a directory that contains **exactly one** genome FASTA
(`.fa`/`.fasta`/`.fna`) and **exactly one** annotation GTF (`.gtf`) — both are
auto-detected. The pipeline builds the RSEM/STAR index (`star_reference/`) and
the Picard reference artifacts (`picard_reference/`: sequence dictionary,
refFlat, rRNA intervals) inside this directory. Gene names are read from the
GTF (`gene_name`/`gene`/`gene_symbol`/`Name` attributes), so the pipeline is not
tied to any single organism.

### Reports — Quarto (optional)

Set `reports_dir` to a directory of your own Quarto reports. Each `.qmd` is
expected one level down, i.e. `<reports_dir>/<name>/<name>.qmd`. They are
rendered **last**, after the rest of the pipeline finishes (so a report can read
any pipeline output), and the rendered HTML will be under
`<results_dir>/reports/<name>/<name>.html`. This step is skipped entirely if
`reports_dir` is omitted from the config or contains no `.qmd` files. A minimal
starter report is in [`example/reports/summary/summary.qmd`](example/reports/summary/summary.qmd).
The motivation for this step is to allow producing figures for publication
as part of the normal pipeline run.

## Development

**TODO:**

- Unit tests
- Support single-end read samples
- Handle custom study design (e.g., paired samples)
- Autogenerate exploratory plots
- Pipeline integration test
