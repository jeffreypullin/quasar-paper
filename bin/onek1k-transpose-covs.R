#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
covs <- read_delim(args[[1]], show_col_types = FALSE)
cell_type <- args[[2]]

t_covs <- covs |>
  pivot_longer(cols = -sample_id) |>
  pivot_wider(names_from = sample_id, values_from = value)
colnames(t_covs)[[1]] <- ""

write_tsv(
  t_covs,
  paste0("onek1k-", cell_type, "-t-all-covs.tsv")
)
