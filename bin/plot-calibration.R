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

compute_qq_data <- function(files, type, method) {

  if (type == "quasar") {
    read_f <- function(x) read_tsv(x)$pvalue
  } else if (type %in% c("tensorqtl", "jaxqtl")) {
    read_f <- function(x) read_parquet(x)$pval_nominal
  } else if (type == "apex") {
    read_f <- function(x) read_tsv(x)$pval
  }

  pvalue <- unlist(map(files, read_f))
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

# Variant level analysis.

quasar_variant_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-quasar)")),
  ) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method != "quasar-nb_glmm") |>
  summarise(file_list = list(quasar_file), .by = c(method, cell_type)) |>
  rowwise() |>
  mutate(qq_data = list(compute_qq_data(file_list, "quasar", method))) |>
  ungroup() |>
  mutate(
    method = factor(method_lookup[method]),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  ) |>
  select(-file_list) |>
  unnest(cols = qq_data)

p <- quasar_variant_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue, colour = cell_type)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  facet_wrap(~ method) +
  scale_colour_manual(values = cell_type_cols) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)",
    colour = "Cell type"
  ) +
  theme_jp()

ggsave(
  "/home/jp2045/quasar-paper/output/plot-quasar-variant-fdr.pdf",
  p,
  width = 10,
  height = 8
)

tensorqtl_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  summarise(file_list = list(tensorqtl_file), .by = c(method, cell_type)) |>
  rowwise() |>
  mutate(qq_data = list(compute_qq_data(file_list, "tensorqtl", method))) |>
  ungroup()

jaxqtl_data <- tibble(jaxqtl_file = to_r_vec(args[3])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=jaxqtl-).*?(?=-chr)"),
  ) |>
  mutate(method = "jaxqtl") |>
  summarise(file_list = list(jaxqtl_file), .by = c(method, cell_type)) |>
  rowwise() |>
  mutate(qq_data = list(compute_qq_data(file_list, "jaxqtl", method))) |>
  ungroup()

apex_data <- tibble(apex_file = to_r_vec(args[4])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)"),
    cell_type = str_extract(apex_file, "(?<=apex-).*?(?=-chr)"),
    method = "apex"
  ) |>
  summarise(file_list = list(apex_file), .by = c(method, cell_type)) |>
  rowwise() |>
  mutate(qq_data = list(compute_qq_data(file_list, "apex", method))) |>
  ungroup()

other_variant_data <- bind_rows(jaxqtl_data, tensorqtl_data, apex_data) |>
  mutate(
   method = factor(method_lookup[method]),
   cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  ) |>
  select(-file_list) |>
  unnest(cols = qq_data)

p1 <- quasar_variant_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue, colour = cell_type)) +
  geom_point(alpha = 0.8) + 
  geom_abline(linetype = "dashed") + 
  facet_wrap(~ method, ncol = 2) + 
  coord_cartesian(ylim = c(0, 7.5)) +
  scale_colour_manual(values = cell_type_cols) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log1o(p-value)",
    colour = "Cell type"
  ) +
  theme_jp()

p2 <- other_variant_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue, colour = cell_type)) +
  geom_point(alpha = 0.8) + 
  geom_abline(linetype = "dashed") + 
  facet_wrap(~ method, ncol = 1) + 
  scale_colour_manual(values = cell_type_cols) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log10(p-value)",
    colour = "Cell type"
  ) +
  theme_jp()

p <- p1 + p2 + 
  plot_layout(guides = "collect", widths = c(2, 1)) + 
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    legend.direction = "vertical",
    legend.position = "right"
  )

ggsave(
  "plot-variant-calibration.pdf", 
  p,
  width = 12,
  height = 10
)

# Gene level analysis.

quasar_gene_data <- tibble(quasar_file = to_r_vec(args[5])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-quasar)")),
  ) |>
  filter(model != "p_glm") |>
  filter(model != "nb_glmm") |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(quasar_file)$pvalue)) |>
  ungroup() |>
  unnest(cols = pvalue) |>
  mutate(
    y_pvalue = sort(pvalue),
    x_pvalue = 1:n() / (n() + 1),
    log_x_pvalue = -log10(x_pvalue),
    x_bin = cut(log_x_pvalue, breaks = seq(0, 3, by = 0.05), include.lowest = TRUE),
    .by = c(method, cell_type)
  ) |>
  summarise(log_y_pvalue = -log10(mean(y_pvalue)), .by = c(method, cell_type, x_bin)) |>
  mutate(log_x_bin_mid = seq(0.05, 2.95, by = 0.05)[as.numeric(x_bin)]) |>
  mutate(
   method = factor(method_lookup[method]),
   cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

tensorqtl_gene_data <- tibble(tensorqtl_file = to_r_vec(args[6])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(tensorqtl_file)$pval_beta)) |>
  ungroup() |>
  unnest(cols = pvalue)

jaxqtl_gene_data <- tibble(jaxqtl_file = to_r_vec(args[7])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=jaxqtl-).*?(?=-chr)"),
  ) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = c("chr", "cell_type")) |>
  mutate(method = "jaxqtl") |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(jaxqtl_file)$pval_beta)) |>
  ungroup() |>
  unnest(cols = pvalue)

apex_gene_data <- tibble(apex_file = to_r_vec(args[8])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)"),
    cell_type = str_extract(apex_file, "(?<=apex-).*?(?=-chr)"),
    method = "apex"
   ) |>
   rowwise() |>
   mutate(pvalue = list(read_tsv(apex_file)$egene_pval)) |>
   ungroup() |>
   unnest(cols = pvalue)

other_plot_data <- bind_rows(
  tensorqtl_gene_data |>
    select(method, cell_type, pvalue),
  jaxqtl_gene_data |>
    select(method, cell_type, pvalue),
  apex_gene_data |>
    select(method, cell_type, pvalue)
) |>
  mutate(
    y_pvalue = sort(pvalue),
    x_pvalue = 1:n() / (n() + 1),
    log_x_pvalue = -log10(x_pvalue),
    x_bin = cut(log_x_pvalue, breaks = seq(0, 3, by = 0.05), include.lowest = TRUE),
    .by = c(method, cell_type)
  ) |>
  summarise(log_y_pvalue = -log10(mean(y_pvalue)), .by = c(method, cell_type, x_bin)) |>
  mutate(log_x_bin_mid = seq(0.05, 2.95, by = 0.05)[as.numeric(x_bin)]) |>
  mutate(
   method = factor(method_lookup[method]),
   cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

p1 <- quasar_gene_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue, colour = cell_type)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  facet_wrap(~ method, ncol = 2) +
  coord_cartesian(ylim = c(0, 4.5)) +
  scale_colour_manual(values = cell_type_cols) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log1o(p-value)",
    colour = "Cell type"
  ) +
  theme_jp()

p2 <- other_plot_data |>
  ggplot(aes(log_x_bin_mid, log_y_pvalue, colour = cell_type)) +
  geom_point(alpha = 0.8) +
  geom_abline(linetype = "dashed") +
  facet_wrap(~ method, ncol = 1) +
  coord_cartesian(ylim = c(0, 4.5)) +
  scale_colour_manual(values = cell_type_cols) +
  labs(
    x = "Expected -log10(p-value)",
    y = "Observed -log1o(p-value)",
    colour = "Cell type"
  ) +
  theme_jp()

p <- p1 + p2 +
  plot_layout(guides = "collect", widths = c(2, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    legend.direction = "vertical",
    legend.position = "right"
  )

ggsave(
  "plot-gene-calibration.pdf",
  p,
  width = 12,
  height = 10
)
