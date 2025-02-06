#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
raw_feat_anno_data_path <- args[[1]]

raw_feat_annot_data <- read_tsv(raw_feat_anno_data_path)

feat_anno_data <- raw_feat_annot_data |>
  select(-ind) |>
  mutate(chrom = str_sub(chrom, 4, -1)) |>
  filter(chrom %in% 1:22)

write_tsv(feat_anno_data, "feat-anno-data.tsv")
