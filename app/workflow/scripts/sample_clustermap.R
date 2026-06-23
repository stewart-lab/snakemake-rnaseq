suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(pheatmap)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts", type = "character", help = "Counts matrix CSV (genes x samples)"),
  make_option("--design", type = "character", default = NULL,
              help = "Optional design CSV: sample,condition (for annotation)"),
  make_option("--out",    type = "character", default = "clustermap.png",
              help = "Output plot path [default %default]")
)))

counts <- as.matrix(read.csv(opt$counts, row.names = 1, check.names = FALSE))

# Optional condition annotation, aligned to the count columns.
annotation <- NA
if (!is.null(opt$design)) {
  design <- read.csv(opt$design, stringsAsFactors = FALSE)
  rownames(design) <- design$sample
  design <- design[colnames(counts), , drop = FALSE]
  annotation <- data.frame(condition = factor(design$condition),
                           row.names = colnames(counts))
}

# VST, then sample-to-sample Euclidean distances.
dds <- DESeqDataSetFromMatrix(round(counts), data.frame(row.names = colnames(counts)), ~ 1)
# vsd <- vst(dds, blind = TRUE)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
dists <- dist(t(assay(vsd)))
mat <- as.matrix(dists)

# Clustering uses the distances directly (samples that cluster = similar).
pheatmap(
  mat,
  clustering_distance_rows = dists,
  clustering_distance_cols = dists,
  annotation_col = annotation,
  main = "Sample-to-sample distances (VST)",
  filename = opt$out,
  width = 8, height = 7
)
message(sprintf("[clustermap] wrote %s", opt$out))