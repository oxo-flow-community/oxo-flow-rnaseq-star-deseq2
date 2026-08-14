#!/usr/bin/env Rscript
# Annotate Ensembl gene ids with gene symbols via biomaRt (with mirror
# fallback). Port of gene2symbol.R from snakemake-workflows/rna-seq-star-
# deseq2 (v3.1.1) with identical logic.
#
# Usage: gene2symbol.R <input.tsv> <output.tsv> <species> <log>
#   species: biomaRt dataset suffix, e.g. "hsapiens" (upstream
#            get_bioc_species_name(): homo_sapiens -> hsapiens)
args <- commandArgs(trailingOnly = TRUE)
in_file <- args[[1]]
out_file <- args[[2]]
species <- args[[3]]
log_file <- args[[4]]

dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
log <- file(log_file, open = "wt")
sink(log)
sink(log, type = "message")

library(biomaRt)
library(tidyverse)
# useful error messages upon aborting
library("cli")

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
while ( class(mart)[[1]] != "Mart" ) {
  mart <- tryCatch(
    {
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(species, "_gene_ensembl"),
        mirror = mart
      )
    },
    error = function(e) {
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        cli_abort(
          str_c(
            "Have tried all 4 available Ensembl biomaRt mirrors ",
            rounds,
            " times. You might have a connection problem, or no mirror is responsive.\n",
            "The last error message was:\n",
            message(e)
          )
        )
      }
      # hop to next mirror
      mart <- switch(mart,
                     useast = "uswest",
                     uswest = "asia",
                     asia = "www",
                     www = {
                       # wait before starting another round through the mirrors,
                       # hoping that intermittent problems disappear
                       Sys.sleep(30)
                       "useast"
                     }
              )
    }
  )
}


df <- read.table(in_file, sep='\t', header=1)

g2g <- biomaRt::getBM(
            attributes = c( "ensembl_gene_id",
                            "external_gene_name"),
            filters = "ensembl_gene_id",
            values = df$gene,
            mart = mart,
            )

annotated <- merge(df, g2g, by.x="gene", by.y="ensembl_gene_id")
annotated$gene <- ifelse(annotated$external_gene_name == '', annotated$gene, annotated$external_gene_name)
annotated$external_gene_name <- NULL
write.table(annotated, out_file, sep='\t', row.names=F)
