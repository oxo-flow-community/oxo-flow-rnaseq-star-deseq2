#!/usr/bin/env Rscript
# PCA plot of rlog-transformed counts for one grouping variable.
#
# Port of plot-pca.R from snakemake-workflows/rna-seq-star-deseq2 (v3.1.1)
# with identical logic; the snakemake object is replaced by CLI arguments.
#
# Usage: plot-pca.R <all.rds> <out.svg> <variable> <log>
args <- commandArgs(trailingOnly = TRUE)
in_rds <- args[[1]]
out_svg <- args[[2]]
variable <- args[[3]]
log_file <- args[[4]]

dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
log <- file(log_file, open = "wt")
sink(log)
sink(log, type = "message")

library("DESeq2")

# load deseq2 data
dds <- readRDS(in_rds)

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
svg(out_svg)
plotPCA(counts, intgroup = variable)
dev.off()
