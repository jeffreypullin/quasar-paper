#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(arrow)
  library(patchwork)
  library(tidyr)
  library(purrr)
  library(forcats)
  library(glue)
})

#FIXME: Make this more reproducible.
source("/home/jp2045/quasar-paper/code/plot-utils.R")
args <- commandArgs(trailingOnly = TRUE)

quasar_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-cis)")),
  ) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(quasar_file)$pvalue)) |>
  ungroup() |>
  unnest(cols = pvalue)

tensorqtl_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_parquet(tensorqtl_file)$pval_nominal)) |>
  ungroup() |>
  unnest(cols = pvalue)

jaxqtl_data <- tibble(jaxqtl_file = to_r_vec(args[3])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=jaxqtl-).*?(?=-chr)"),
  ) |>
  mutate(method = "jaxqtl") |>
  rowwise() |>
  mutate(pvalue = list(read_parquet(jaxqtl_file)$pval_nominal)) |>
  ungroup() |>
  unnest(cols = pvalue)

apex_data <- tibble(apex_file = to_r_vec(args[4])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)"),
    cell_type = str_extract(apex_file, "(?<=apex-).*?(?=-chr)"),
    method = "apex"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(apex_file)$pval)) |>
  ungroup() |>
  unnest(cols = pvalue)

other_plot_data <- bind_rows(
  jaxqtl_data |>
    select(method, cell_type, pvalue),
  tensorqtl_data |>
    select(method, cell_type, pvalue),
  apex_data |>
    select(method, cell_type, pvalue)
)

p <- other_plot_data |>
  ggplot(aes(pvalue)) +
  geom_histogram() +
  facet_wrap(~method + cell_type)

ggsave("plot-other-fdr.pdf", p)

p <- quasar_data |>
  ggplot(aes(pvalue)) +
  geom_histogram() +
  facet_wrap(~method + cell_type)

#ggsave("plot-quasar-fdr.pdf", p)
