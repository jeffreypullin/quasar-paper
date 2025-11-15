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

quasar_variant_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-quasar)")),
  ) |>
  filter(cell_type %in% c("Plasma", "B IN", "CD4 NC")) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method != "quasar-nb_glmm") |>
  filter(method != "quasar-p_glm") |>
  expand_grid(alpha_level = c(1e-2, 1e-3)) |>
  rowwise() |>
  mutate(pvalues = list(fread(quasar_file, select = "pvalue")$pvalue)) |>
  mutate(
    n = length(pvalues),
    n_alpha = sum(pvalues < alpha_level),
    prop = n_alpha / n
  ) |>
  ungroup()

tensorqtl_variant_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  expand_grid(alpha_level = c(1e-2, 1e-3)) |>
  rowwise() |>
  mutate(pvalues = list(read_parquet(tensorqtl_file, col_select = "pval_nominal")$pval_nominal)) |>
  mutate(
    n = length(pvalues),
    n_alpha = sum(pvalues < alpha_level),
    prop = n_alpha / n
  ) |>
  ungroup()

jaxqtl_variant_data <- tibble(jaxqtl_file = to_r_vec(args[3])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=jaxqtl-).*?(?=-chr)"),
  ) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = c("chr", "cell_type")) |>
  mutate(method = "jaxqtl") |>
  expand_grid(alpha_level = c(1e-2, 1e-3)) |>
  rowwise() |>
  mutate(n = sum(map_dbl(
    jaxqtl_file,
    function(x) length(read_parquet(x, col_select = "pval_nominal")$pval_nominal)
  ))) |>
  mutate(n_alpha = sum(map_dbl(
    jaxqtl_file,
    function(x) sum(read_parquet(x, col_select = "pval_nominal")$pval_nominal < alpha_level)
  ))) |>
  mutate(prop = n_alpha / n) |>
  ungroup()

variant_plot_data <- bind_rows(
  quasar_variant_data |>
    summarise(
      prop = mean(prop), n = sum(n),
      .by = c(method, cell_type, alpha_level)
    ),
  tensorqtl_variant_data |>
    summarise(
      prop = mean(prop), n = sum(n),
      .by = c(method, cell_type, alpha_level)
    ),
  jaxqtl_variant_data |>
    summarise(
      prop = mean(prop), n = sum(n),
      .by = c(method, cell_type, alpha_level)),
) |>
  mutate(
    method = fct_reorder(factor(method_lookup[method]), prop),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  ) |>
  select(-ends_with("file"))

data_hline <- tibble(
  alpha_level = c(1e-2, 1e-3),
)

p1 <- variant_plot_data |>
  ggplot(aes(method, prop, fill = cell_type)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_hline(
    data = data_hline,
    aes(yintercept = alpha_level),
    linetype = "dotted"
  ) +
  scale_fill_manual(values = cell_type_cols) +
  facet_wrap(~alpha_level, scales = "free") +
  coord_flip() +
  labs(
    y = "Proportion of variants",
    x = "Method",
    fill = "Cell type"
  ) +
  theme_jp_vgrid() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "right"
  )

# Gene level analysis

quasar_gene_data <- tibble(quasar_file = to_r_vec(args[5])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-quasar)")),
  ) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method != "quasar-p_glm") |>
  filter(method != "quasar-nb_glmm") |>
  expand_grid(alpha_level = c(0.1, 0.01)) |>
  rowwise() |>
  mutate(pvalues = list(read_tsv(quasar_file)$pvalue)) |>
  mutate(
    n = length(pvalues),
    n_alpha = sum(pvalues < alpha_level),
    prop = n_alpha / n
  ) |>
  ungroup()

tensorqtl_gene_data <- tibble(tensorqtl_file = to_r_vec(args[6])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  expand_grid(alpha_level = c(0.1, 0.01)) |>
  rowwise() |>
  mutate(pvalues = list(read_tsv(tensorqtl_file)$pval_beta)) |>
  mutate(
    n = length(pvalues),
    n_alpha = sum(pvalues < alpha_level),
    prop = n_alpha / n
  ) |>
  ungroup()

jaxqtl_gene_data <- tibble(jaxqtl_file = to_r_vec(args[7])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=jaxqtl-).*?(?=-chr)"),
  ) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = c("chr", "cell_type")) |>
  mutate(method = "jaxqtl") |>
  expand_grid(alpha_level = c(0.1, 0.01)) |>
  rowwise() |>
  mutate(pvalues = list(read_tsv(jaxqtl_file)$pval_beta)) |>
  mutate(
    n = length(pvalues),
    n_alpha = sum(pvalues < alpha_level),
    prop = n_alpha / n
  ) |>
  ungroup()

gene_plot_data <- bind_rows(
  quasar_gene_data |>
    filter(cell_type %in% c("Plasma", "B IN", "CD4 NC")) |>
    summarise(
      prop = mean(prop), n = sum(n),
      .by = c(method, cell_type, alpha_level)
    ),
  tensorqtl_gene_data |>
    summarise(
      prop = mean(prop), n = sum(n),
      .by = c(method, cell_type, alpha_level)
    ),
  jaxqtl_gene_data |>
    summarise(
      prop = mean(prop), n = sum(n),
      .by = c(method, cell_type, alpha_level)
    ),
) |>
  mutate(
    method = fct_reorder(factor(method_lookup[method]), prop),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  ) |>
  select(-ends_with("file"))

data_hline <- tibble(
  alpha_level = c(0.1, 0.01)
)

p2 <- gene_plot_data |>
  ggplot(aes(method, prop, fill = cell_type)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_hline(
    data = data_hline,
    aes(yintercept = alpha_level),
    linetype = "dotted"
  ) +
  scale_fill_manual(values = cell_type_cols) +
  facet_wrap(~alpha_level, scales = "free") +
  coord_flip() +
  labs(
    y = "Proportion of genes",
    x = "Method",
    fill = "Cell type"
  ) +
  theme_jp_vgrid() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    legend.position = "right"
  )

p <- p1 / p2 +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    legend.direction = "vertical",
    legend.position = "right"
  )

ggsave(
  "plot-fdr.pdf",
  p,
  width = 12,
  height = 10
)