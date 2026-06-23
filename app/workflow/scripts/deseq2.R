#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts",    type = "character", help = "Counts matrix CSV (genes x samples)"),
  make_option("--design",    type = "character", help = "Design CSV: sample,condition"),
  make_option("--contrasts", type = "character", help = "Contrasts CSV: control,treatment"),
  make_option("--outdir",    type = "character", help = "Output directory")
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Read inputs -------------------------------------------------------------
counts    <- as.matrix(read.csv(opt$counts, row.names = 1, check.names = FALSE))
design    <- read.csv(opt$design, stringsAsFactors = FALSE)
contrasts <- read.csv(opt$contrasts, stringsAsFactors = FALSE)

# Align design rows to count columns and make condition a factor
rownames(design) <- design$sample
design <- design[colnames(counts), , drop = FALSE]
design$condition <- factor(design$condition)
stopifnot(!any(is.na(design$condition)))

# ---- Build and run DESeq2 ----------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = design,
  design    = ~ condition
)
dds <- DESeq(dds)

# ---- Size factors (one CSV: sample,size_factor) ------------------------------
sf <- sizeFactors(dds)
sf_df <- data.frame(sample = names(sf), size_factor = as.numeric(sf))
sf_out <- file.path(opt$outdir, "size_factors.csv")
write.csv(sf_df, sf_out, row.names = FALSE)
message(sprintf("[deseq2] wrote %s", sf_out))

# ---- Normalized counts matrix (keyed by ensembl gene ID) ---------------------
norm_counts <- counts(dds, normalized = TRUE)

write.csv(as.data.frame(norm_counts),
          file.path(opt$outdir, "normalized_counts.csv"), row.names = TRUE)
message(sprintf("[deseq2] wrote %s",
                file.path(opt$outdir, "normalized_counts.csv")))

# ---- One results table per contrast ------------------------------------------
for (i in seq_len(nrow(contrasts))) {
  ctrl <- contrasts$control[i]
  trt  <- contrasts$treatment[i]

  res <- results(dds, contrast = c("condition", trt, ctrl))
  res <- res[order(res$padj), ]
  res_df <- as.data.frame(res)

  # Append normalized counts for this contrast's control + treatment samples.
  ctrl_samples <- rownames(design)[design$condition == ctrl]
  trt_samples  <- rownames(design)[design$condition == trt]
  keep_samples <- c(ctrl_samples, trt_samples)
  res_df <- cbind(res_df,
                  norm_counts[rownames(res_df), keep_samples, drop = FALSE])

  fname <- sprintf("%s_vs_%s.csv", trt, ctrl)
  write.csv(res_df, file.path(opt$outdir, fname), row.names = TRUE)

  message(sprintf("[deseq2] wrote %s (%d genes)",
                  file.path(opt$outdir, fname), nrow(res_df)))
}
