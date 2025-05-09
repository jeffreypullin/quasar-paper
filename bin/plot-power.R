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
  mutate(n_sig = sum(read_tsv(quasar_file)$pvalue < 1e-6)) |>
  ungroup()

tensorqtl_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(n_sig = sum(read_parquet(tensorqtl_file)$pval_nominal < 1e-6)) |>
  ungroup()

jaxqtl_data <- tibble(jaxqtl_file = to_r_vec(args[3])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=jaxqtl-).*?(?=-chr)"),
  ) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = c("chr", "cell_type")) |>
  mutate(method = "jaxqtl") |>
  rowwise() |>
  mutate(n_sig = sum(map_dbl(
    jaxqtl_file,
    function(x)  sum(read_parquet(x)$pval_nominal < 1e-6)))
  ) |>
  ungroup()

apex_data <- tibble(apex_file = to_r_vec(args[4])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)"),
    cell_type = str_extract(apex_file, "(?<=apex-).*?(?=-chr)"),
    method = "apex"
  ) |>
  rowwise() |>
  mutate(n_sig = sum(read_tsv(apex_file)$pval < 1e-6)) |>
  ungroup()

plot_data <- bind_rows(
  apex_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  quasar_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  jaxqtl_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  tensorqtl_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
) |>
  filter(method != "quasar-p_glm") |>
  mutate(method = fct_reorder(factor(method), n_sig))

power_plot <- plot_data |>
  ggplot(aes(method, n_sig, fill = cell_type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  theme_bw()

ggsave("plot-power.pdf", power_plot)
