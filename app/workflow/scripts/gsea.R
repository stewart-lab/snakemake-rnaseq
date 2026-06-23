#!/usr/bin/env Rscript

# Gene Set Enrichment Analysis (GSEA) for a single contrast using fgsea.
# Genes are ranked by a DESeq2 statistic (the Wald statistic by default) and the
# provided gene sets are tested for enrichment among the up/down-regulated genes.

suppressPackageStartupMessages({
  library(optparse)
  library(fgsea)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--deseq",     type = "character",
              help = "Gene-symbol DESeq2 results CSV for one contrast (012_gene_symbols/<contrast>.csv)"),
  make_option("--gene_sets", type = "character", help = "CSV listing GMT files (one per row)"),
  make_option("--out",       type = "character", help = "Output GSEA results CSV"),
  make_option("--rank_by",   type = "character", default = "stat",
              help = "DESeq2 column to rank genes by [default %default]"),
  make_option("--min_size",  type = "integer", default = 15L,
              help = "Minimum gene set size to test [default %default]"),
  make_option("--max_size",  type = "integer", default = 500L,
              help = "Maximum gene set size to test [default %default]")
)))

dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)

# ---- Read GMT files and merge into one named list of gene sets ----------------
gmt_table <- read.csv(opt$gene_sets, stringsAsFactors = FALSE)
gmt_files <- if ("gmt" %in% names(gmt_table)) gmt_table$gmt else gmt_table[[1]]
stopifnot(length(gmt_files) > 0)
missing_gmt <- gmt_files[!file.exists(gmt_files)]
if (length(missing_gmt)) {
  stop("GMT file(s) not found: ", paste(missing_gmt, collapse = ", "))
}

gene_sets <- list()
for (f in gmt_files) {
  sets <- fgsea::gmtPathways(f)
  dup <- intersect(names(sets), names(gene_sets))
  if (length(dup)) {
    message(sprintf("[gsea] %d duplicate gene set name(s) in %s ignored (kept first)",
                    length(dup), basename(f)))
    sets <- sets[!names(sets) %in% names(gene_sets)]
  }
  gene_sets <- c(gene_sets, sets)
}
message(sprintf("[gsea] loaded %d gene sets from %d GMT file(s)",
                length(gene_sets), length(gmt_files)))

# ---- Build the ranked gene list ----------------------------------------------
df <- read.csv(opt$deseq, check.names = FALSE)
stopifnot("gene_symbol" %in% names(df))
if (!opt$rank_by %in% names(df)) {
  stop(sprintf("rank_by column '%s' not found in %s", opt$rank_by, opt$deseq))
}

ranks <- df[[opt$rank_by]]
names(ranks) <- df$gene_symbol
ranks <- ranks[!is.na(ranks)]

# Collapse duplicate gene symbols to their most extreme statistic, then sort
# descending (fgsea expects a ranked, named vector with unique names).
ranks <- ranks[order(abs(ranks), decreasing = TRUE)]
ranks <- ranks[!duplicated(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)
message(sprintf("[gsea] ranked %d unique gene symbols by '%s'",
                length(ranks), opt$rank_by))

# ---- Run fgsea ----------------------------------------------------------------
res <- fgsea(pathways = gene_sets, stats = ranks,
             minSize = opt$min_size, maxSize = opt$max_size)
res <- res[order(res$padj), ]

# leadingEdge is a list column; flatten to a '/'-joined string so it fits a CSV.
res <- as.data.frame(res)
res$leadingEdge <- vapply(res$leadingEdge,
                          function(x) paste(x, collapse = "/"), character(1))

write.csv(res, opt$out, row.names = FALSE)
message(sprintf("[gsea] wrote %s (%d gene sets tested)", opt$out, nrow(res)))
