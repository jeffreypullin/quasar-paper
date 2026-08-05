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
  library(data.table)
})

# FIXME: Make this more reproducible.
source("/home/jp2045/quasar-paper/code/plot-utils.R")
args <- commandArgs(trailingOnly = TRUE)

compute_genomic_inflation <- function(file, type) {
  
  if (type == "quasar") {
    read_f <- function(x) read_tsv(x, show_col_types = FALSE)$pvalue
  } else if (type %in% c("tensorqtl", "jaxqtl")) {
    read_f <- function(x) read_parquet(x)$pval_nominal
  } else if (type == "apex") {
    read_f <- function(x) read_tsv(x)$pval
  }

  pvalue <- read_f(file)
  n <- length(pvalue)
  x_pvalue <- 1:n / (n + 1)
  obs_chisq <- qchisq(pvalue, df = 1, lower.tail = FALSE)
  exp_chisq <- qchisq(x_pvalue, df = 1, lower.tail = FALSE)
  lambda <- median(obs_chisq) / median(exp_chisq)
  lambda
}

quasar_variant_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-quasar)")),
  ) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method != "quasar-nb_glmm") |>
  rowwise() |>
  mutate(gen_inf = compute_genomic_inflation(quasar_file, "quasar")) |>
  ungroup() |>
  summarise(gen_inf = median(gen_inf), .by = c(cell_type, method)) |>
  mutate(
    method = factor(method_lookup[method]),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

tensorqtl_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(gen_inf = compute_genomic_inflation(tensorqtl_file, "tensorqtl")) |>
  ungroup() |>
  summarise(gen_inf = median(gen_inf), .by = c(cell_type, method)) |>
  mutate(
    method = factor(method_lookup[method]),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

jaxqtl_data <- tibble(jaxqtl_file = to_r_vec(args[3])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=jaxqtl-).*?(?=-chr)"),
  ) |>
  mutate(method = "jaxqtl") |>
  rowwise() |>
  mutate(gen_inf = compute_genomic_inflation(jaxqtl_file, "jaxqtl")) |>
  ungroup() |>
  summarise(gen_inf = median(gen_inf), .by = c(cell_type, method)) |>
  mutate(
    method = factor(method_lookup[method]),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

apex_data <- tibble(apex_file = to_r_vec(args[4])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)"),
    cell_type = str_extract(apex_file, "(?<=apex-).*?(?=-chr)"),
    method = "apex"
  ) |>
  rowwise() |>
  mutate(gen_inf = compute_genomic_inflation(apex_file, "apex")) |>
  ungroup() |>
  summarise(gen_inf = median(gen_inf), .by = c(cell_type, method)) |>
  mutate(
    method = factor(method_lookup[method]),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

variant_data <- bind_rows(
  quasar_variant_data,
  tensorqtl_data,
  jaxqtl_data,
  apex_data
) |>
  mutate(gen_inf = round(gen_inf, 2))

write_tsv(variant_data, "genomic-inflation.tsv")
