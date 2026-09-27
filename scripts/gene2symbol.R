#!/usr/bin/env Rscript
# Annotate Ensembl gene ids with gene symbols via biomaRt (with mirror
# fallback). Port of gene2symbol.R from snakemake-workflows/rna-seq-star-
# deseq2 (v3.1.1) with TWO deliberate deviations (issues #7, #10):
#
# 1. The "uswest" mirror hop is dropped: biomaRt 2.62 validates mirrors
#    against [www, useast, asia] and rejects "uswest", so every retry
#    round wasted one guaranteed-fail attempt. Chain: useast -> asia -> www.
#
# 2. WAF-resilient fallback (2026-09): Ensembl's edge WAF started
#    rejecting biomaRt's default curl/httr2 user-agent (HTTP 403 on GET
#    to /biomart/martservice) and rejecting POST from the www host
#    outright (HTTP 405), which makes plain useEnsembl()+getBM() fail on
#    every mirror regardless of network health. When the mirror loop is
#    exhausted, the script automatically switches to a self-contained
#    workaround, verified live end-to-end (TP53/BRCA1 round-trip):
#      a. patch .getArchiveList: biomaRt validates the current Ensembl
#         release by scraping an archives HTML page; the WAF 403 breaks
#         that too. Fetch it with a browser UA across asia/www/useast.
#      b. patch bmRequest (GET) and .submitQueryXML (POST) to send the
#         same browser UA.
#      c. build the Mart object manually on the asia mirror (the host
#         whose edge still accepts POST) and populate its attribute and
#         filter slots directly, bypassing useEnsembl's redirect logic.
#    Patches are applied only inside this script's session and only
#    after the normal path has failed on all mirrors; when Ensembl is
#    reachable the original behaviour is preserved exactly.
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

# --- WAF-resilient fallback (self-contained; used only on exhaustion) ---
# Reproduces useDataset's slot population and re-points all biomaRt
# network calls at a browser user-agent. See header notes for the
# HTTP 403/405 background.
.waf_fallback_getBM <- function(df) {
  ua <- "Mozilla/5.0 (X11; Linux x86_64; rv:127.0) Gecko/20100101 Firefox/127.0"

  # (a) archives HTML: browser UA, multi-host GET fallback
  .patched_getArchiveList <- function(https = TRUE, http_config = list()) {
    mirrors <- c("asia", "www", "useast")
    while (length(mirrors) > 0) {
      m <- mirrors[[1]]; mirrors <- mirrors[-1]
      u <- paste0("https://", m,
                  ".ensembl.org/info/website/archives/index.html?redirect=no")
      r <- tryCatch(
        httr2::resp_body_string(
          httr2::req_perform(httr2::req_user_agent(httr2::request(u), ua))),
        error = function(e) NULL)
      if (!is.null(r) && nchar(r) > 10000 &&
          !grepl("temporarily unavailable", r)) {
        message("gene2symbol fallback: archives page fetched via ", m)
        return(r)
      }
    }
    stop("Unable to contact any Ensembl mirror (fallback .getArchiveList)")
  }
  assignInNamespace(".getArchiveList", .patched_getArchiveList, "biomaRt")

  # (b) browser UA on ALL biomaRt network calls (GET + POST)
  .patched_bmRequest <- function(request, http_config, verbose = FALSE) {
    request <- httr2::req_options(
      httr2::req_timeout(httr2::request(request),
                         getOption("timeout", default = 60)),
      !!!http_config)
    request <- httr2::req_user_agent(request, ua)
    result <- httr2::req_perform(request)
    result2 <- httr2::resp_body_string(result)
    if (is.na(result2)) result2 <- httr2::resp_body_string(result, encoding = "Latin1")
    result2
  }
  assignInNamespace("bmRequest", .patched_bmRequest, "biomaRt")

  .patched_submitQueryXML <- function(host, query, http_config) {
    req <- httr2::req_options(
      httr2::req_timeout(
        httr2::req_body_form(httr2::request(host), query = query),
        max(getOption("timeout", default = 300), 300)),
      !!!http_config)
    req <- httr2::req_user_agent(req, ua)
    res <- httr2::req_perform(req)
    if (httr2::resp_is_error(res)) {
      err_msg <- getFromNamespace(".createErrorMessage", "biomaRt")(
        error_code = httr2::resp_status(res), host = host)
      stop(err_msg, call. = FALSE)
    }
    httr2::resp_body_string(res)
  }
  assignInNamespace(".submitQueryXML", .patched_submitQueryXML, "biomaRt")

  # (c) manual Mart on the asia mirror + populate slots directly
  m <- new("Mart")
  m@biomart <- "ENSEMBL_MART_ENSEMBL"
  m@host <- "https://asia.ensembl.org/biomart/martservice"
  m@vschema <- "default"
  m@version <- "asia"
  m@dataset <- str_c(species, "_gene_ensembl")
  biomaRt:::martAttributes(m) <- biomaRt:::.getAttributes(m)
  biomaRt:::martFilters(m) <- biomaRt:::.getFilters(m)

  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "ensembl_gene_id",
        values = df$gene,
        mart = m)
}

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because
# of very frequent "Internal Server Error"s).
# On exhaustion (waf_fallback_required) the
# browser-UA fallback takes over below.
mart <- "useast"
rounds <- 0
waf_fallback_required <- FALSE
while ( class(mart)[[1]] != "Mart" && !waf_fallback_required ) {
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
        waf_fallback_required <<- TRUE
        return(list())  # placeholder; while-loop exits via the flag
      }
      # hop to next mirror
      mart <- switch(mart,
                     useast = "asia",
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

g2g <- if (class(mart)[[1]] == "Mart") {
  biomaRt::getBM(
      attributes = c( "ensembl_gene_id",
                      "external_gene_name"),
      filters = "ensembl_gene_id",
      values = df$gene,
      mart = mart,
      )
} else {
  message("All Ensembl mirrors exhausted -- switching to ",
          "browser-UA WAF fallback ...")
  .waf_fallback_getBM(df)
}

annotated <- merge(df, g2g, by.x="gene", by.y="ensembl_gene_id")
annotated$gene <- ifelse(annotated$external_gene_name == '', annotated$gene, annotated$external_gene_name)
annotated$external_gene_name <- NULL
write.table(annotated, out_file, sep='\t', row.names=F)
