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
  library(ggupset)
})

# FIXME: Make this more reproducible.
source("/home/jp2045/quasar-paper/code/plot-utils.R")
args <- commandArgs(trailingOnly = TRUE)

# Variant level (cis-nomianal) data.

quasar_variant_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-cis)")),
  ) |> 
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method != "quasar-p_glm") |>
  rowwise() |>
  mutate(n_sig = sum(fread(quasar_file, select = "pvalue")$pvalue < 1e-6)) |>
  ungroup()

tensorqtl_variant_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(n_sig = sum(read_parquet(tensorqtl_file, col_select = "pval_nominal")$pval_nominal < 1e-6)) |>
  ungroup()

jaxqtl_variant_data <- tibble(jaxqtl_file = to_r_vec(args[3])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=jaxqtl-).*?(?=-chr)"),
  ) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = c("chr", "cell_type")) |>
  mutate(method = "jaxqtl") |>
  rowwise() |>
  mutate(n_sig = sum(map_dbl(
    jaxqtl_file,
    function(x)  sum(read_parquet(x, col_select = "pval_nominal")$pval_nominal < 1e-6, na.rm = TRUE)))
  ) |>
  ungroup()

apex_variant_data <- tibble(apex_file = to_r_vec(args[4])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)"),
    cell_type = str_extract(apex_file, "(?<=apex-).*?(?=-chr)"),
    method = "apex"
  ) |>
  rowwise() |>
  mutate(n_sig = sum(fread(apex_file, select = "pval")$pval < 1e-6)) |>
  ungroup()

# Gene level (cis) data.

quasar_gene_data <- tibble(quasar_file = to_r_vec(args[5])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-cis)")),
  ) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method != "quasar-p_glm") |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(quasar_file)$pvalue)) |>
  ungroup() |>
  unnest(cols = pvalue) |>
  summarise(
    n_sig = sum(p.adjust(pvalue, method = "BH") < 0.01), 
    .by = c(method, cell_type)
  )

tensorqtl_gene_data <- tibble(tensorqtl_file = to_r_vec(args[6])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(
    n_sig = sum(p.adjust(read_tsv(tensorqtl_file)$pval_beta, method = "BH") < 0.01)
  ) |>
  ungroup()

tensorqtl_gene_acat_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl_acat"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_parquet(tensorqtl_file) |>
    summarise(pvalue = acat(pval_nominal), .by = phenotype_id) |>
    pull(pvalue)
  )) |>
  unnest(cols = pvalue) |>
  ungroup() |>
  summarise(
    n_sig = sum(p.adjust(pvalue, method = "BH") < 0.01), 
    .by = c(method, cell_type)
  )

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
  unnest(cols = pvalue) |>
  summarise(
    n_sig = sum(p.adjust(pvalue, method = "BH") < 0.01, na.rm = TRUE), 
    .by = c(method, cell_type)
  )

apex_gene_data <- tibble(apex_file = to_r_vec(args[8])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)"),
    cell_type = str_extract(apex_file, "(?<=apex-).*?(?=-chr)"),
    method = "apex"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(apex_file)$egene_pval)) |>
  ungroup() |>
  unnest(cols = pvalue) |>
  summarise(
    n_sig = sum(p.adjust(pvalue, method = "BH") < 0.01), 
    .by = c(method, cell_type)
  )

# Variant-level analysis.

variant_plot_data <- bind_rows(
  quasar_variant_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  tensorqtl_variant_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  jaxqtl_variant_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  apex_variant_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type))
) |>
  mutate(
    method = fct_reorder(factor(method_lookup[method]), n_sig),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

variant_power_plot <- variant_plot_data |>
  ggplot(aes(method, n_sig, fill = cell_type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = cell_type_cols) +
  guides(colour = guide_legend(ncol = 1)) +
  coord_flip() +
  labs(
    y = "Number of eSNPs",
    x = "Method"
  ) + 
  theme_jp_vgrid() + 
  theme(legend.direction = "vertical")

# Gene level analysis.

gene_plot_data <- bind_rows(
  quasar_gene_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  tensorqtl_gene_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  jaxqtl_gene_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  apex_gene_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  tensorqtl_gene_acat_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type))
) |>
  mutate(
    method = fct_reorder(factor(method_lookup[method]), n_sig),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

gene_power_plot <- gene_plot_data |>
  ggplot(aes(method, n_sig, fill = cell_type)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = cell_type_cols) +
  guides(colour = guide_legend(ncol = 1)) +
  coord_flip() +
  labs(
    y = "Number of eGenes",
    x = ""
  ) + 
  theme_jp_vgrid() + 
  theme(legend.direction = "vertical")

# Upset plot

upset_plot_data <- tibble(quasar_file = to_r_vec(args[5])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-cis)")),
  ) |> 
  mutate(method = paste0("quasar-", model)) |>
  select(-model) |>
  filter(method != "quasar-p_glm") |>
  filter(cell_type == "B IN") |>
  rowwise() |>
  mutate(data = list(fread(quasar_file, select = c("pvalue", "feature_id")))) |>
  ungroup() |>
  summarise(data = list(bind_rows(data)), .by = method) |>
  rowwise() |>
  mutate(gene = list(data |> 
    mutate(pvalue_adjust = p.adjust(pvalue, method = "BH")) |>
    filter(pvalue_adjust < 0.01) |>
    pull(feature_id))
  ) |>
  select(-data) |>
  ungroup() |>
  unnest(cols = gene) |>
  mutate(method = method_lookup[method]) |>
  summarise(method = list(method), .by = gene)

upset_plot <- upset_plot_data |>
  ggplot(aes(x = method)) +
  geom_bar(fill = "#BBBBBB") +
  scale_x_upset(n_intersections = 10) + 
  labs(
    y = "Number of significant genes",
    x = "Combination of methods"
  ) +
  theme_jp() + 
  theme(axis.title = element_text(family = "Helvetica", size = 18, color = "#222222")) + 
  theme_combmatrix(
    combmatrix.label.text = element_text(family = "Helvetica", size = 18, color = "#222222")
  )

power_plot <- variant_power_plot + gene_power_plot + upset_plot +
  plot_layout(guides = "collect", nrow = 2) +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    legend.direction = "vertical"
  )

ggsave(
  "plot-power.pdf", 
  power_plot,
  width = 12,
  height = 12
)
