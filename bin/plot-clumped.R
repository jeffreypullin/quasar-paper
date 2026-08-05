#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tidyr)
  library(forcats)
  library(glue)
  library(data.table)
})

source("/home/jp2045/quasar-paper/code/plot-utils.R")
args <- commandArgs(trailingOnly = TRUE)

quasar_clump_data <- tibble(quasar_clump_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_clump_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_clump_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_clump_file, glue("(?<={cell_type}-).*?(?=-harmonised)")),
  ) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method != "quasar-p_glm") |>
  filter(method != "quasar-nb_glmm") |>
  rowwise() |>
  mutate(clump_data = list(read_tsv(quasar_clump_file, show_col_types = FALSE))) |>
  ungroup() |>
  unnest(cols = c(clump_data))

plot_data <- quasar_clump_data |>
  filter(method %in% c("quasar-nb_glm-apl", "quasar-lm")) |>
  mutate(method = if_else(method == "quasar-lm", "LM", "NB-GLM (APL)")) |>
  summarise(
    n_signals = n_distinct(ID),
    .by = c(method, cell_type, feature_id)
  ) |>
  summarise(
    n_signals = sum(n_signals),
    .by = c(method, cell_type)
  ) |>
  mutate(cell_type = fct_reorder(cell_type, n_signals, .fun = sum))

p_total <- plot_data |>
  ggplot(aes(cell_type, n_signals, fill = method)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = c("#99DDFF", "#44BB99")) +
  labs(
    x = "Cell type",
    y = "Number of independent signals",
    fill = "Method"
  ) + 
  theme_bw()

ggsave(
  "clumped-n-independent-eqtls-plot.pdf",
  plot = p_total,
  width = 10,
  height = 5
)
