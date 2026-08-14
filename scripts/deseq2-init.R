#!/usr/bin/env Rscript
# Build the DESeq2 dataset, run DESeq and write normalized counts.
#
# Port of deseq2-init.R from snakemake-workflows/rna-seq-star-deseq2 (v3.1.1)
# with identical logic; the snakemake object is replaced by CLI arguments.
#
# Usage: deseq2-init.R <counts.tsv> <samples.tsv> <all.rds> <normcounts.tsv>
#                       <model> <vofs> <base_levels> <batch_effects>
#                       <threads> <log>
#   vofs:         comma-joined variables_of_interest (e.g. treatment_1,treatment_2)
#   base_levels:  comma-joined base level per variable_of_interest (same order)
#   batch_effects: comma-joined batch effects (empty string if none)
#   model:        design formula, e.g. "~condition"; empty builds the default
#                 (additive batch effects + all vof interactions)
args <- commandArgs(trailingOnly = TRUE)
counts_file <- args[[1]]
samples_file <- args[[2]]
out_rds <- args[[3]]
out_norm <- args[[4]]
design_model <- args[[5]]
vofs <- strsplit(args[[6]], ",", fixed = TRUE)[[1]]
base_levels <- strsplit(args[[7]], ",", fixed = TRUE)[[1]]
batch_effects <- strsplit(args[[8]], ",", fixed = TRUE)[[1]]
threads <- as.integer(args[[9]])
log_file <- args[[10]]

dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
log <- file(log_file, open = "wt")
sink(log)
sink(log, type="message")

library(stringr)
library("DESeq2")

parallel <- FALSE
if (threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(threads))
    parallel <- TRUE
}

counts_data <- read.table(
  counts_file,
  header = TRUE,
  row.names = "gene",
  check.names = FALSE
)
counts_data <- counts_data[, order(names(counts_data))]

col_data <- read.table(
  samples_file,
  header = TRUE,
  row.names = "sample_name",
  check.names = FALSE
)
col_data <- col_data[order(row.names(col_data)), , drop = FALSE]

# properly set the base level to the configuration, avoiding
# the default behaviour of choosing the alphabetical minimum level
for (i in seq_along(vofs)) {
  col_data[[vofs[i]]] <- relevel(
    factor(col_data[[vofs[i]]]), base_levels[i]
  )
}

# properly turn all batch effects into factors, even if they are numeric
for (effect in batch_effects) {
  if (str_length(effect) > 0) {
    col_data[[effect]] <- factor(col_data[[effect]])
  }
}

# build up formula with additive batch_effects and all interactions between the
# variables_of_interest
if (str_length(design_model) == 0) {
  batch_effects_str <- str_flatten(batch_effects, " + ")
  if (str_length(batch_effects_str) > 0) {
    batch_effects_str <- str_c(batch_effects_str, " + ")
  }
  vof_interactions <- str_flatten(vofs, " * ")
  design_model <- str_c("~", batch_effects_str, vof_interactions)
}

dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = as.formula(design_model)
)

# remove uninformative columns
dds <- dds[rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel = parallel)

# Write dds object as RDS
saveRDS(dds, file = out_rds)
# Write normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.table(
  data.frame(
    "gene" = rownames(norm_counts),
    norm_counts
  ),
  file = out_norm,
  sep = "\t",
  row.names = FALSE
)
