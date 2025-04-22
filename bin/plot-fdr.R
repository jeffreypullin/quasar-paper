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
})

#FIXME: Make this more reproducible.
source("/home/jp2045/quasar-paper/code/plot-utils.R")
args <- commandArgs(trailingOnly = TRUE)

quasar_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(chr = str_extract(quasar_file, "chr[0-9]+")) |>
  mutate(model = str_extract(quasar_file, "(?<=-)[^-]+(?=-cis)")) |>
  mutate(cell_type = "B IN") |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model)

jaxqtl_data <- tibble(jaxqtl_file = to_r_vec(args[2])) |>
  mutate(chr = str_extract(jaxqtl_file, "chr[0-9]+")) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = "chr")

plot_data <- jaxqtl_data |>
  unnest(cols = jaxqtl_file) |>
  rowwise() |>
  mutate(pvalue = list(read_parquet(jaxqtl_file)$pval_nominal)) |>
  ungroup() |>
  unnest(cols = pvalue)

p <- plot_data |>
  ggplot(aes(pvalue)) +
  geom_histogram()

ggsave("plot-jaxqtl-fdr.pdf", p)

plot_data <- quasar_data |>
  mutate(ind = seq_len(n()), .by = method) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(quasar_file)$pvalue)) |>
  ungroup() |>
  unnest(cols = pvalue)

p <- plot_data |>
  ggplot(aes(pvalue)) +
  geom_histogram() +
  facet_wrap(~method)

ggsave("plot-quasar-fdr.pdf", p)
