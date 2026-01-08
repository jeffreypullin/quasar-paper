#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

id_file <- args[1]
cov_file <- args[2]
pheno_file <- args[3]
grm_file <- args[4]
pb_type <- args[5]

ids <- read_csv(id_file, col_names = FALSE)[[1]]

covs <- read_tsv(cov_file)
filt_covs <- filter(covs, sample_id %in% ids)
write_tsv(filt_covs, "cov-n100.tsv")

pheno <- read_tsv(pheno_file)
# Filter the individuals and to the first 20 genes on chromosome 22.
filt_pheno <- pheno |>
  select(all_of(c("#chr", "start", "end", "phenotype_id", ids))) |>
  filter(`#chr` == 22) |>
  slice_head(n = 20)
write_tsv(filt_pheno, paste0(pb_type, "-pheno-n100.bed"))

grm <- read_tsv(grm_file)
filt_grm <- grm |>
  select(all_of(c("sample_id", ids))) |>
  filter(sample_id %in% ids) |>
  mutate(sample_id = factor(sample_id, levels = ids)) %>%
  arrange(sample_id)

write_tsv(filt_grm, "grm-n100.tsv")
