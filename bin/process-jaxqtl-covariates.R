#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
covs <- read_delim(args[[1]])
cell_type <- args[[2]]

jaxqtl_covs <- covs |>
  rename(iid = sample_id) |>
  select(-int)

write_tsv(
  jaxqtl_covs,
  paste0("onek1k-", cell_type, "-all-jaxqtl-covs.tsv")
)
