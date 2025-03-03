#!/usr/bin/env Rscript

library(readr)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
prefix <- args[[1]]

grm_id <- read_tsv("king_ibd_out.king.id")
raw_grm <- as.matrix(read_tsv("king_ibd_out.king", col_names = FALSE))

grm <- as.matrix(nearPD(2 * raw_grm, corr = TRUE)$mat)
tmp <- cbind(grm_id$IID, grm)
colnames(tmp) <- c("sample_id", grm_id$IID)
tmp <- as.data.frame(tmp)

write_tsv(tmp, paste0(prefix, "-grm.tsv"))
