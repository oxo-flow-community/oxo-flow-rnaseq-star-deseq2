#!/usr/bin/env Rscript
# DESeq2 results for one contrast with ashr shrinkage + MA plot.
#
# Port of deseq2.R from snakemake-workflows/rna-seq-star-deseq2 (v3.1.1)
# with identical logic; the snakemake object is replaced by CLI arguments.
# Both upstream contrast forms are supported: the list-form
# (variable_of_interest / level_of_interest / base_level) and the
# string-form (an arbitrary R expression, eval(parse(...)) — upstream
# deseq2.R's `typeof(contrast_config) == "character"` branch).
#
# Usage: deseq2.R <all.rds> <diffexp.tsv> <ma-plot.svg> <vofs> <base_levels>
#                  <contrast_ids> <contrast_variables> <contrast_levels>
#                  <contrast_expressions> <threads> <log>
#   vofs/base_levels: as in deseq2-init.R
#   contrast_ids:     comma-joined contrast ids (config.contrasts). The id
#                     of THIS instance is derived from <diffexp.tsv>'s
#                     basename (the engine does not substitute scatter
#                     variables inside script fields, so the output-table
#                     filename carries the id).
#   contrast_variables / contrast_levels: comma-joined per-contrast values.
#   contrast_expressions: ";;"-joined per-contrast R expressions; an empty
#                     entry selects the list-form contrast above.
args <- commandArgs(trailingOnly = TRUE)
in_rds <- args[[1]]
out_table <- args[[2]]
out_ma <- args[[3]]
vofs <- strsplit(args[[4]], ",", fixed = TRUE)[[1]]
base_levels <- strsplit(args[[5]], ",", fixed = TRUE)[[1]]
contrast_ids <- strsplit(args[[6]], ",", fixed = TRUE)[[1]]
contrast_variables <- strsplit(args[[7]], ",", fixed = TRUE)[[1]]
contrast_levels <- strsplit(args[[8]], ",", fixed = TRUE)[[1]]
contrast_expressions <- strsplit(args[[9]], ";;", fixed = TRUE)[[1]]
threads <- as.integer(args[[10]])
log_file <- args[[11]]

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

# Identify this instance's contrast: the engine renders scatter variables
# inside script fields literally, so the id is read back from the output
# table's filename (results/diffexp/<id>.diffexp.tsv).
contrast_id <- sub("\\.diffexp\\.tsv$", "", basename(out_table))
i <- match(contrast_id, contrast_ids)
if (is.na(i)) {
    cli_abort(
        c(
            "main.oxoflow: could not match this contrast's id against config.contrasts.",
            "x" = "Contrast id derived from the output file: {contrast_id}",
            " " = "Configured contrasts:",
            " " = "{contrast_ids}",
            "i" = "The output file name must be <contrast-id>.diffexp.tsv, one instance per config.contrasts entry."
        )
    )
}

if (nchar(contrast_expressions[[i]]) > 0) {
    # more complex contrast specification via list(c(), c()), see ?results
    # docs of the DESeq2 package and this tutorial (plus the linked
    # seqanswers thread):
    # https://github.com/tavareshugo/tutorial_DESeq2_contrasts/blob/main/DESeq2_contrasts.md
    contrast <- eval(parse(text = contrast_expressions[[i]]))
} else {
    contrast_vof <- contrast_variables[[i]]
    contrast_level <- contrast_levels[[i]]

    # check for existence of the contrast's variable_of_interest to provide
    # a useful error message
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
