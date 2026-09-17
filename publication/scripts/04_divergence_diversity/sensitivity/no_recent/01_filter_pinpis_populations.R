#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: script nuclear_population.tsv mt13_population.tsv output_directory")
excluded <- c("SC", "CH", "LB", "PACH", "FRED")
expected <- c("RS","FG","LG","SL","SR","TL","WB","WK","WT",
              "SAY","AMO","BEA","BOOT","ECHO","GOS","JOE","LAW",
              "MUC","PYE","ROB","SWA","THE")
dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)

filter_table <- function(path, label) {
  x <- fread(path, na.strings = c("NA", "NaN", ""))
  if (!"population" %in% names(x)) stop(label, " lacks population column")
  x[, population := toupper(trimws(population))]
  if (!all(excluded %in% unique(x$population))) stop(label, ": recent populations missing")
  y <- x[!population %chin% excluded]
  if (!setequal(unique(y$population), expected)) stop(label, ": retained set is not the expected 22")
  y
}

nu <- filter_table(args[[1L]], "nuclear")
mt <- filter_table(args[[2L]], "mitochondrial")
fwrite(nu, file.path(args[[3L]], "nuclear_pinpis_population.no_recent.tsv"), sep = "\t")
fwrite(mt, file.path(args[[3L]], "mt13_pinpis_population.no_recent.tsv"), sep = "\t")
audit <- rbind(
  data.table(dataset="nuclear", populations=uniqueN(nu$population), rows=nrow(nu)),
  data.table(dataset="mitochondrial", populations=uniqueN(mt$population), rows=nrow(mt)))
fwrite(audit, file.path(args[[3L]], "population_filter_audit.tsv"), sep = "\t")
print(audit)
