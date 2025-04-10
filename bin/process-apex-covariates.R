#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
t_covs <- read_delim(args[[1]], show_col_types = FALSE)
cell_type <- args[[2]]

apex_covs <- t_covs
apex_covs <- apex_covs[-1, ]
colnames(apex_covs)[[1]] <- "#ID"

write_tsv(
  apex_covs,
  paste0("onek1k-", cell_type, "-all-apex-covs.tsv")
)
