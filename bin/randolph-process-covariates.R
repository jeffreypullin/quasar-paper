#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

raw_cov_data_path <- args[[1]]

raw_cov_data <- read_delim(args[[1]])
cov_data <- raw_cov_data |>
  filter(str_sub(sample_id, -2, -1) ==  "NI") |>
  mutate(ethnicity = if_else(ethnicity == "EUR", 0, 1)) |>
  mutate(int = 1) |>
  select(
    sample_id,
    int,
    age_scale = age_Scale,
    ethnicity
  )

write_tsv(cov_data, "cov-data.tsv")
