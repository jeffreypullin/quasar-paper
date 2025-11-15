#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
  library(tidyr)
  library(purrr)
  library(glue)
  library(forcats)
  library(scales)
  library(data.table)
})

# FIXME: Make this more reproducible.
source("/home/jp2045/quasar-paper/code/plot-utils.R")
args <- commandArgs(trailingOnly = TRUE)

compute_qq_data_quasar <- function(files) {

  pvalue <- unlist(map(files, function(x) read_tsv(x)$pvalue))
  n <- length(pvalue)
  y_pvalue <- sort(pvalue)
  x_pvalue <- 1:n / (n + 1)
  log_x_pvalue <- -log10(x_pvalue)
  x_bin <- cut(
    log_x_pvalue,
    breaks = seq(0, 6, by = 0.1),
    include.lowest = TRUE
  )

  tibble(x_bin, y_pvalue) |>
    summarise(log_y_pvalue = -log10(mean(y_pvalue)), .by = x_bin) |>
    mutate(log_x_bin_mid = seq(0.05, 5.95, by = 0.1)[as.numeric(x_bin)])
}

filt_perm_variant_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-quasar)")),
  ) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method == "quasar-nb_glm-apl") |>
  summarise(file_list = list(quasar_file), .by = c(method, cell_type)) |>
  rowwise() |>
  mutate(qq_data = list(compute_qq_data_quasar(file_list))) |>
  ungroup() |>
  filter(cell_type == "CD4 NC") |>
  mutate(method = factor(method_lookup[method])) |>
  select(-file_list) |>
  unnest(cols = qq_data) |>
  mutate(type = "CoV filtered")

no_filt_perm_variant_data <- tibble(quasar_file = to_r_vec(args[2])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-quasar)")),
  ) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method == "quasar-nb_glm-apl") |>
  summarise(file_list = list(quasar_file), .by = c(method, cell_type)) |>
  rowwise() |>
  mutate(qq_data = list(compute_qq_data_quasar(file_list))) |>
  ungroup() |>
  filter(cell_type == "CD4 NC") |>
  mutate(method = factor(method_lookup[method])) |>
  select(-file_list) |>
  unnest(cols = qq_data) |>
  mutate(type = "Unfiltered")

perm_variant_data <- bind_rows(
    filt_perm_variant_data, 
    no_filt_perm_variant_data
)

p <- perm_variant_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue, colour = type)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  scale_colour_manual(values = c("#228833", "#AA3377")) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)",
    colour = ""
  ) +
  theme_jp()

ggsave(
  "plot-filter.pdf",
  p,
  width = 10,
  height = 6
)
