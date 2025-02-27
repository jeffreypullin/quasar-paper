#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
expr_covs <- read_delim(args[[1]])
geno_pcs <- read_delim(args[[2]])
cell_type <- args[[3]]

all_covs <- expr_covs |>
  left_join(geno_pcs, by = "sample_id") |>
  mutate(
    int = 1,
    sex = sex - 1,
  ) |>
  select(int, sex, age, paste0("PC_", 1:5), paste0("geno_pc", 1:5)) |>
  rename_with(
    ~ paste0("expr_pc", str_sub(.x, -1), recycle0 = TRUE),
    starts_with("PC_")
  ) |>
  rename_with(
    ~ paste0("geno_pc", str_sub(.x, -1), recycle0 = TRUE),
    starts_with("geno_pc")
  )

write_tsv(
  all_covs,
  paste0("onek1k-", cell_type, "-all-covs.tsv")
)

