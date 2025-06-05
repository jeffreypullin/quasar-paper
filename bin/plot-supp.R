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
  library(qvalue)
})

# FIXME: Make this more reproducible.
source("/home/jp2045/quasar-paper/code/plot-utils.R")
args <- commandArgs(trailingOnly = TRUE)

# Plot effect of APL on phi.

phi_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-cis)")),
  ) |> 
  filter(model %in% c("nb_glm", "nb_glm-apl")) |>
  rowwise() |>
  mutate(phi = list(fread(quasar_file, select = c("phi", "feature_id")) |>
                      summarise(phi = first(phi), .by = feature_id) |>
                      pull(phi))
  ) |>
  ungroup()

p <- phi_data |>
  unnest(cols = phi) |>
  select(-c(quasar_file, chr)) |>
  filter(phi < 10) |>
  mutate(use_apl = if_else(model == "nb_glm", "MLE", "APL")) |>
  ggplot(aes(cell_type, log10(phi), fill = use_apl)) + 
  scale_fill_manual(values = c("#009988", "#EE7733")) +
  geom_boxplot() + 
  labs(
    x = "Cell type",
    y = "log10(phi)"
  ) + 
  theme_jp()

ggsave(
  "plot-phi.pdf", 
  p,
  width = 10,
  height = 8
)

tensorqtl_acat_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_parquet(tensorqtl_file, col_select = c("phenotype_id", "pval_nominal")) |>
    summarise(pvalue = acat(pval_nominal), .by = phenotype_id) |>
    pull(pvalue))
  ) |>
  ungroup() |>
  unnest(cols = pvalue)

tensorqtl_gene_data <- tibble(tensorqtl_file = to_r_vec(args[3])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(tensorqtl_file)$pval_beta)) |>
  ungroup() |>
  unnest(cols = pvalue)

perm_p <- tensorqtl_gene_data |>
  ggplot(aes(pvalue)) + 
  geom_histogram(breaks = seq(-0.05, 1.05, by = 0.05)) + 
  coord_cartesian(xlim = c(0, 1)) + 
  labs(
    x = "Permutation p-value",
    y = ""
  ) + 
  theme_jp()

acat_p <- tensorqtl_acat_data |>
  ggplot(aes(pvalue)) + 
  geom_histogram(breaks = seq(-0.05, 1.05, by = 0.05)) + 
  coord_cartesian(xlim = c(0, 1)) + 
  labs(
    x = "ACAT p-value",
    y = ""
  ) + 
  theme_jp()

p <- perm_p + acat_p + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(
  "plot-null-hist.pdf", 
  p,
  width = 10,
  height = 8
)
