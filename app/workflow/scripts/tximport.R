#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tximport)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--files", type = "character",
              help = "Space- or comma-separated .genes.results files"),
  make_option("--out", type = "character",
              help = "Output counts matrix (TSV)")
)))

files <- strsplit(opt$files, "[, ]+")[[1]]

names(files) <- sub("\\.genes\\.results$", "", basename(files))
stopifnot(all(file.exists(files)))

# RSEM .genes.results are already gene-level, so txIn/txOut = FALSE.
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

df <- data.frame(feature_id = rownames(txi$counts), txi$counts,
                 check.names = FALSE)
write.table(df, opt$out, sep = ",", quote = FALSE, row.names = FALSE)