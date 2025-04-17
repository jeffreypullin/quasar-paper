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

plot_data <- quasar_data |>
  rowwise() |>
  mutate(n_sig = read_tsv(quasar_file, show_col_types = FALSE) |>
    mutate(sig = pvalue < 1e-6) |>
    count(sig) |>
    pull(n) |>
    pluck(2)
 ) |>
 ungroup() |>
 summarise(n_sig = sum(n_sig), .by = method) |>
 print()

power_plot <- plot_data |>
  mutate(method = fct_reorder(factor(method), n_sig)) |>
  ggplot(aes(method, n_sig)) +
  geom_col() +
  coord_flip() +
  theme_bw()

ggsave("plot-power.pdf", power_plot)
