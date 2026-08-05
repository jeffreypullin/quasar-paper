#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(glue)
})

args <- commandArgs(trailingOnly = TRUE)
chr <- args[[1]]
cell_type <- args[[2]]
variant_file <- args[[3]]

clump_files <- list.files(pattern = "*.clumps")

clumps_data <- tibble(file_name = clump_files) |>
  mutate(feature_id = str_extract(file_name, "ENSG[0-9]+")) |>
  rowwise() |>
  mutate(data = list(read_tsv(file_name, show_col_types = FALSE))) |>
  unnest(cols = c(data)) |>
  select(-file_name)

model <- str_extract(variant_file, glue("(?<={cell_type}-).*?(?=-quasar)"))

write_tsv(
  clumps_data, 
  paste0(chr, "-", cell_type, "-", model, "-harmonised-clumps.tsv")
)
