#!/usr/bin/env Rscript

library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

grm <- read_delim(args[[1]], show_col_types = FALSE)
ids <- grm$sample_id
grm_mat <- grm[1:nrow(grm), 2:ncol(grm)]
grm_mat <- as.matrix(grm_mat)

id_col_1 <- character()
id_col_2 <- character()
val_col <- numeric()
n <- 1
for (i in seq_len(nrow(grm_mat))) {
  for (j in seq_len(ncol(grm_mat))) {
    id_col_1[[n]] <- ids[[i]]
    id_col_2[[n]] <- ids[[j]]
    val_col[[n]] <- grm_mat[i, j]
    n <- n + 1
  }
}

sparse_grm <- tibble(
  `#id1` = id_col_1,
  `id2` = id_col_2,
  grm = val_col
)

write_tsv(sparse_grm, "onek1k-sparse-grm.tsv")
