#!/usr/bin/env Rscript

library(readr)

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[[1]]

grm_id <- read_tsv("king_ibd_out.king.id")
grm_mat <- read_tsv("king_ibd_out.king", col_names = FALSE)

tmp <- 2 * grm_mat
tmp <- cbind(grm_id$IID, tmp)
colnames(tmp) <- c("sample_id", grm_id$IID)
tmp <- as.data.frame(tmp)

write_tsv(tmp, paste0(prefix, "-grm.tsv"))
