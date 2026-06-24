#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GSVA)
  library(limma)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--expr",      type = "character",
              help = "Gene-symbol-keyed expression CSV (012_gene_symbols/normalized_counts.csv)"),
  make_option("--design",    type = "character", help = "Design CSV: sample,condition"),
  make_option("--contrasts", type = "character", help = "Contrasts CSV: control,treatment"),
  make_option("--gene_sets", type = "character", help = "CSV listing GMT files (one per row)"),
  make_option("--kcdf",      type = "character", default = "Gaussian",
              help = "GSVA kernel for the empirical CDF: Gaussian, Poisson, or none [default %default]"),
  make_option("--tau",       type = "double", default = 1,
              help = "GSVA tau weighting exponent [default %default]"),
  make_option("--outdir",    type = "character", help = "Output directory")
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Read GMT files and merge into one named list of gene sets ----------------
read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  sets <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    genes <- fields[-(1:2)]              # drop set name + description columns
    unique(genes[nzchar(genes)])
  })
  names(sets) <- vapply(lines, function(line) strsplit(line, "\t")[[1]][1],
                        character(1))
  sets
}

gmt_table <- read.csv(opt$gene_sets, stringsAsFactors = FALSE)
gmt_files <- if ("gmt" %in% names(gmt_table)) gmt_table$gmt else gmt_table[[1]]
stopifnot(length(gmt_files) > 0)
missing_gmt <- gmt_files[!file.exists(gmt_files)]
if (length(missing_gmt)) {
  stop("GMT file(s) not found: ", paste(missing_gmt, collapse = ", "))
}

gene_sets <- list()
for (f in gmt_files) {
  sets <- read_gmt(f)
  dup <- intersect(names(sets), names(gene_sets))
  if (length(dup)) {
    message(sprintf("[gsva_limma] %d duplicate gene set name(s) in %s ignored (kept first): %s",
                    length(dup), basename(f), paste(head(dup, 5), collapse = ", ")))
    sets <- sets[!names(sets) %in% names(gene_sets)]
  }
  gene_sets <- c(gene_sets, sets)
}
message(sprintf("[gsva_limma] loaded %d gene sets from %d GMT file(s)",
                length(gene_sets), length(gmt_files)))

# ---- Read expression (gene symbols x samples) --------------------------------
expr_df <- read.csv(opt$expr, check.names = FALSE)
stopifnot("gene_symbol" %in% names(expr_df))
sample_cols <- setdiff(names(expr_df), c("gene_symbol", "gene_id"))
mat <- as.matrix(expr_df[, sample_cols, drop = FALSE])
symbols <- expr_df$gene_symbol

# Collapse duplicate gene symbols to the row with the highest mean expression.
if (anyDuplicated(symbols)) {
  ord <- order(rowMeans(mat), decreasing = TRUE)
  keep <- sort(ord[!duplicated(symbols[ord])])
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
}
rownames(mat) <- symbols

# GSVA's Gaussian kcdf expects continuous, log-like values; normalized counts
# are right-skewed, so log2-transform before scoring.
mat <- log2(mat + 1)

# ---- Read design and contrasts -----------------------------------------------
design <- read.csv(opt$design, stringsAsFactors = FALSE)
rownames(design) <- design$sample
contrasts <- read.csv(opt$contrasts, stringsAsFactors = FALSE)

# ---- GSVA over all samples (one scores matrix, like normalized_counts) -------
par <- gsvaParam(exprData = mat, geneSets = gene_sets, kcdf = opt$kcdf, tau = opt$tau)
es_all <- gsva(par, verbose = FALSE)

scores_out <- file.path(opt$outdir, "gsva_scores.csv")
write.csv(as.data.frame(es_all), scores_out, row.names = TRUE)
message(sprintf("[gsva_limma] wrote %s (%d gene sets x %d samples)",
                scores_out, nrow(es_all), ncol(es_all)))

# ---- limma per contrast: stats + this contrast's GSVA scores -----------------
for (i in seq_len(nrow(contrasts))) {
  ctrl <- contrasts$control[i]
  trt  <- contrasts$treatment[i]

  ctrl_samples <- rownames(design)[design$condition == ctrl]
  trt_samples  <- rownames(design)[design$condition == trt]
  keep_samples <- intersect(c(ctrl_samples, trt_samples), colnames(es_all))
  es_sub <- es_all[, keep_samples, drop = FALSE]

  # Reference = control, so coef 2 is treatment - control.
  cond <- factor(design[keep_samples, "condition"], levels = c(ctrl, trt))
  fit <- eBayes(lmFit(es_sub, model.matrix(~ cond)))
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "P")

  # Append this contrast's GSVA scores, mirroring the DESeq2 results tables.
  out <- cbind(res, es_sub[rownames(res), , drop = FALSE])

  fname <- sprintf("%s_vs_%s.csv", trt, ctrl)
  write.csv(out, file.path(opt$outdir, fname), row.names = TRUE)
  message(sprintf("[gsva_limma] wrote %s (%d gene sets)",
                  file.path(opt$outdir, fname), nrow(out)))
}
