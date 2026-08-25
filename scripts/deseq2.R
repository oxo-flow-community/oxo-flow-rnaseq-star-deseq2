#!/usr/bin/env Rscript
# DESeq2 results for one contrast with ashr shrinkage + MA plot.
#
# Port of deseq2.R from snakemake-workflows/rna-seq-star-deseq2 (v3.1.1)
# with identical logic; the snakemake object is replaced by CLI arguments.
#
# Usage: deseq2.R <all.rds> <diffexp.tsv> <ma-plot.svg> <vofs> <base_levels>
#                  <contrast.vof> <contrast.level> <contrast_exprs>
#                  <contrasts> <contrast_id> <threads> <log>
#   vofs/base_levels: as in deseq2-init.R
#   contrast.vof:     variable_of_interest of the contrast
#   contrast.level:   level_of_interest of the contrast (base level is looked
#                     up from base_levels)
#   contrast_exprs:   semicolon-joined string-form contrast R expressions,
#                     parallel to `contrasts` (upstream diffexp.contrasts
#                     string entries); empty entry = list-form contrast
#   contrasts:        comma-joined contrast ids (to index contrast_exprs)
#   contrast_id:      id of this contrast instance
args <- commandArgs(trailingOnly = TRUE)
in_rds <- args[[1]]
out_table <- args[[2]]
out_ma <- args[[3]]
vofs <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
base_levels <- strsplit(args[[5]], ",", fixed = TRUE)[[1]]
contrast_vof <- args[[6]]
contrast_level <- args[[7]]
contrast_exprs <- strsplit(args[[8]], ";", fixed = TRUE)[[1]]
contrasts_all <- strsplit(args[[9]], ",", fixed = TRUE)[[1]]
contrast_id <- args[[10]]
threads <- as.integer(args[[11]])
log_file <- args[[12]]

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

# string-form contrast specification via list(c(), c()), see ?results docs
# of the DESeq2 package (upstream: a single-character diffexp.contrasts
# entry is eval(parse(text = ...))'d verbatim)
idx <- match(contrast_id, contrasts_all)
contrast_expr <- if (!is.na(idx) && idx <= length(contrast_exprs)) {
    contrast_exprs[[idx]]
} else {
    ""
}
if (nzchar(contrast_expr)) {
    contrast <- eval(parse(text = contrast_expr))
} else {
    # basic case of contrast specification (upstream list form): check for
    # existence of the contrast's variable_of_interest to provide a useful
    # error message
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
}

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
