#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
pheno_bed <- read_delim(args[[1]], show_col_types = FALSE)
chr <- args[[2]]
cell_type <- args[[3]]
split_n <- as.numeric(args[[4]])

chr_num <- as.numeric(str_sub(chr, 4))
pheno_bed <- pheno_bed |>
  filter(`#chr` == chr_num)

gene_ids <- pheno_bed$phenotype_id
gene_index <- seq_along(gene_ids)
chunk_index <- ceiling(gene_index / split_n)
for (i in seq_along((unique(chunk_index)))) {
  ind <- unique(chunk_index)[[i]]
  chunk_mask <- which(chunk_index == ind)
  chunk  <- tibble(chunk = gene_ids[chunk_mask])
  chunk_pos <- gene_index[chunk_mask]
  chunk_filename <- paste0(
    "chunk-", min(chunk_pos), "-", max(chunk_pos), ".tsv"
  )
  write_tsv(
    chunk,
    chunk_filename,
    col_names = FALSE
  )
}
