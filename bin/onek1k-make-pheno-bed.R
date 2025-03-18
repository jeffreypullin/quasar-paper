#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
pheno <- read_delim(args[[1]], show_col_types = FALSE)
# FIXME: Issue with inconsistent column name/row delimiter.
feat_anno <- read_delim(
  args[[2]],
  col_names = FALSE,
  skip = 1,
  show_col_types = FALSE
)
colnames(feat_anno) <- c("feature_id", "chrom", "start", "end")
cell_type <- args[[3]]
pb_type <- args[[4]]

pheno_bed <- inner_join(feat_anno, pheno, by = "feature_id") |>
  rename(phenotype_id = feature_id, `#chr` = chrom) |>
  relocate(phenotype_id, .after = end) |>
  arrange(`#chr`, start)

write_tsv(
  pheno_bed,
  paste0("onek1k-", cell_type, "-pheno-", pb_type, ".bed")
)
