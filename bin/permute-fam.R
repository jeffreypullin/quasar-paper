#!/usr/bin/env Rscript

library(readr)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
fam <- args[[1]]
chr <- args[[2]]

fam_data <- read_delim(fam, col_names = FALSE)
perm <- order(runif(nrow(fam_data)))
perm_fam_data <- fam_data[perm, ]

write_delim(
  perm_fam_data,
  paste0("permute-", chr, ".fam"),
  col_names = FALSE
)
