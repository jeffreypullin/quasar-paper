#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
pheno_bed <- read_delim(args[[1]], show_col_types = FALSE)
cell_type <- args[[2]]

raw_values <- pheno_bed[, 5:ncol(pheno_bed)]

normalized_values <- apply(raw_values, 1, function(x) {
  as.vector(scale(qnorm(rank(x) / (length(x) + 1))))
})

pheno_bed[, 5:ncol(pheno_bed)] <- t(normalized_values)

write_tsv(
  pheno_bed,
  paste0("onek1k-", cell_type, "-pheno-norm.bed"),
)
