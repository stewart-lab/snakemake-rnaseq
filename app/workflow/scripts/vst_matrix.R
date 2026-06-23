#!/usr/bin/env Rscript

# VST-transform a raw count matrix and write the transformed matrix to CSV so
# it can be consumed downstream (e.g. by the Plotly PCA script in Python).

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts", type = "character",
              help = "Counts matrix CSV (genes x samples)"),
  make_option("--out",    type = "character", default = "vst_matrix.csv",
              help = "Output VST matrix CSV (genes x samples) [default %default]")
)))

counts <- as.matrix(read.csv(opt$counts, row.names = 1, check.names = FALSE))

# No design needed for a blind VST; ~ 1 fits an intercept-only model.
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = data.frame(row.names = colnames(counts)),
  design    = ~ 1
)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

write.csv(as.data.frame(assay(vsd)), opt$out, row.names = TRUE)
message(sprintf("[vst] wrote %s (%d genes x %d samples)",
                opt$out, nrow(vsd), ncol(vsd)))
