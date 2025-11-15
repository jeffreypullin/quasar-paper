#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(stringr)
library(purrr)

source("/home/jp2045/quasar-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)
print(args[1])
phenos <- to_r_vec(args[1])

n_indivs <- map_dbl(phenos, function(x) ncol(read_tsv(x)) - 4)

file_name <- basename(phenos)
cell_type <- str_extract(file_name, "(?<=onek1k-).*?(?=-pheno-)")
pb_type <- str_extract(file_name, "(?<=-pheno-).*?(?=\\.bed)")

n_indiv_data <- tibble(
  cell_type = cell_type,
  pb_type = pb_type,
  n_indiv = n_indivs
) |>
  filter(pb_type == "sum") |>
  select(-pb_type)

write_tsv(n_indiv_data, "n-indiv.tsv")
