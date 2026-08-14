#!/usr/bin/env Rscript
# DESeq2 results for one contrast with ashr shrinkage + MA plot.
#
# Port of deseq2.R from snakemake-workflows/rna-seq-star-deseq2 (v3.1.1)
# with identical logic; the snakemake object is replaced by CLI arguments.
#
# Usage: deseq2.R <all.rds> <diffexp.tsv> <ma-plot.svg> <vofs> <base_levels>
#                  <contrast.vof> <contrast.level> <threads> <log>
#   vofs/base_levels: as in deseq2-init.R
#   contrast.vof:     variable_of_interest of the contrast
#   contrast.level:   level_of_interest of the contrast (base level is looked
#                     up from base_levels)
args <- commandArgs(trailingOnly = TRUE)
in_rds <- args[[1]]
out_table <- args[[2]]
out_ma <- args[[3]]
vofs <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
base_levels <- strsplit(args[[5]], ",", fixed = TRUE)[[1]]
contrast_vof <- args[[6]]
contrast_level <- args[[7]]
threads <- as.integer(args[[8]])
log_file <- args[[9]]

dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
log <- file(log_file, open = "wt")
sink(log)
sink(log, type = "message")

library("cli")
library("DESeq2")

parallel <- FALSE
if (threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(threads))
    parallel <- TRUE
}

dds <- readRDS(in_rds)

# check for existence of the contrast's variable_of_interest to provide a
# useful error message
if (!(contrast_vof %in% vofs)) {
  cli_abort(
    c(
      "main.oxoflow: All contrast variables must also exist under diffexp_variables.",
      "x" = "Could not find variable_of_interest: {contrast_vof}",
      " " = "It was not among the configured diffexp_variables:",
      " " = "{vofs}",
      "i" = "Are there any typos in the contrasts' variable_of_interest entries?"
    )
  )
}
base_level <- base_levels[match(contrast_vof, vofs)]
contrast <- c(contrast_vof, contrast_level, base_level)

res <- results(
  dds,
  contrast = contrast,
  parallel = parallel
)
# shrink fold changes for lowly expressed genes
# use ashr so we can use `contrast` as conversion to coef is not trivial
res <- lfcShrink(
  dds,
  contrast = contrast,
  res = res,
  type = "ashr"
)

# sort by p-value
res <- res[order(res$padj), ]

# store results
svg(out_ma)
plotMA(res, ylim = c(-2, 2))
dev.off()

write.table(
  data.frame(
    "gene" = rownames(res),
    res
  ),
  file = out_table,
  row.names = FALSE,
  sep = "\t"
)
