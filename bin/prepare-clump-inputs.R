#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

var_data <- read_tsv(args[[1]], show_col_types = FALSE)

var_data <- var_data |>
  filter(!is.na(pvalue))

for (fid in unique(var_data$feature_id)) {
  out_data <- var_data |>
    filter(feature_id == fid) |>
    select(ID = snp_id, P = pvalue)
  write_tsv(
    out_data,
    paste0("assoc-", fid, ".tsv")
  )
}
