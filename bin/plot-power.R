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

count_data <- fread(args[9]) |> 
  as_tibble() |>
  select(-c(`#chr`, start, end)) |>
  mutate(sum = rowSums(across(!starts_with("phenotype_id")))) |>
  select(gene = phenotype_id, sum)

upset_plot_data <- tibble(quasar_file = to_r_vec(args[5])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-cis)")),
  ) |> 
  filter(model != "p_glm") |>
  filter(cell_type == "B IN") |>
  rowwise() |>
  mutate(data = list(fread(quasar_file, select = c("pvalue", "feature_id")))) |>
  ungroup() |>
  summarise(data = list(bind_rows(data)), .by = model) |>
  rowwise() |>
  mutate(gene = list(data |> 
    mutate(pvalue_adjust = p.adjust(pvalue, method = "BH")) |>
    filter(pvalue_adjust < 0.01) |>
    pull(feature_id))
  ) |>
  select(-data) |>
  ungroup() |>
  unnest(cols = gene) |>
  mutate(model = model_lookup[model]) |>
  summarise(model = list(model), .by = gene) |>
  left_join(count_data, by = "gene")

upset_plot <- upset_plot_data |>
  ggplot(aes(x = model)) +
  geom_bar(fill = "#BBBBBB") +
  scale_x_upset(n_intersections = 10) + 
  labs(
    y = "Number of eGenes",
    x = "Combination of methods"
  ) +
  theme_jp() + 
  theme(axis.title = element_text(family = "Helvetica", size = 18, color = "#222222")) + 
  theme_combmatrix(
    combmatrix.label.text = element_text(family = "Helvetica", size = 18, color = "#222222")
  )

total_boxplot <- upset_plot_data |> 
  rowwise() |>
  mutate(type = case_when(
    length(model) == 5 ~ "all",
    length(model) == 3 && 
    "P GLMM" %in% model && 
    "NB GLM" %in% model &&
    "NB GLM, APL" %in% model ~ "count",
    .default = "other"
  )) |>
  ungroup() |>
  filter(type %in% c("all", "count")) |>
  mutate(type = if_else(type == "all", "All models", "Count models only")) |>
  ggplot(aes(factor(type), log10(sum))) +
  geom_boxplot() +
  coord_flip() + 
  labs(
    y = "log10(Total counts)",
    x = "Type of eGene"
  ) +
  theme_jp_vgrid()

power_plot <- ((variant_power_plot + gene_power_plot) + plot_layout(guides = "collect")) | (upset_plot / total_boxplot) + 
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    legend.direction = "vertical"
  )

write_csv(gene_plot_data, "gene-power-data.csv")
write_csv(variant_plot_data, "variant-power-data.csv")

ggsave(
  "plot-power.pdf", 
  power_plot,
  width = 15,
  height = 12
)
