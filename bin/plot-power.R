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
  library(scales)
  library(qvalue)
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
  filter(method != "quasar-nb_glmm") |>
  rowwise() |>
  mutate(n_sig = sum(fread(quasar_file, select = "pvalue")$pvalue < 5e-6)) |>
  ungroup()

tensorqtl_variant_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(n_sig = sum(read_parquet(tensorqtl_file, col_select = "pval_nominal")$pval_nominal < 5e-6)) |>
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
    function(x)  sum(read_parquet(x, col_select = "pval_nominal")$pval_nominal < 5e-6, na.rm = TRUE)))
  ) |>
  ungroup()

apex_variant_data <- tibble(apex_file = to_r_vec(args[4])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)"),
    cell_type = str_extract(apex_file, "(?<=apex-).*?(?=-chr)"),
    method = "apex"
  ) |>
  rowwise() |>
  mutate(n_sig = sum(fread(apex_file, select = "pval")$pval < 5e-6)) |>
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
  filter(method != "quasar-nb_glmm") |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(quasar_file)$pvalue)) |>
  ungroup() |>
  unnest(cols = pvalue) |>
  summarise(
    n_sig = sum(p.adjust(pvalue, method = "BH") < 0.05),
    .by = c(method, cell_type)
  )

tensorqtl_gene_data <- tibble(tensorqtl_file = to_r_vec(args[6])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(tensorqtl_file)$pval_beta)) |>
  unnest(cols = pvalue) |>
  ungroup() |>
  summarise(
    n_sig = sum(p.adjust(pvalue, method = "BH") < 0.05),
    .by = c(method, cell_type)
  )

tensorqtl_gene_qvalue_data <- tibble(tensorqtl_file = to_r_vec(args[6])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(tensorqtl_file)$pval_beta)) |>
  unnest(cols = pvalue) |>
  ungroup() |>
  summarise(
    n_sig = sum(qvalue(pvalue)$qvalues < 0.05),
    .by = c(method, cell_type)
  )

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
    n_sig = sum(p.adjust(pvalue, method = "BH") < 0.05),
    .by = c(method, cell_type)
  )

tensorqtl_gene_acat_qvalue_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl_acat_qvalue"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_parquet(tensorqtl_file) |>
    summarise(pvalue = acat(pval_nominal), .by = phenotype_id) |>
    pull(pvalue)
  )) |>
  unnest(cols = pvalue) |>
  ungroup() |>
  summarise(
    n_sig = sum(qvalue(pvalue)$qvalues < 0.05),
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
    n_sig = sum(qvalue(pvalue)$qvalues < 0.05, na.rm = TRUE),
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
    n_sig = sum(p.adjust(pvalue, method = "BH") < 0.05), 
    .by = c(method, cell_type)
  )

# Variant-level analysis.

variant_plot_data <- bind_rows(
  quasar_variant_data |>
    filter(cell_type %in% c("Plasma", "B IN", "CD4 NC")) |>
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
  geom_col() +
  scale_fill_manual(values = cell_type_cols) +
  scale_y_continuous(labels = label_number(scale_cut = cut_long_scale())) +
  guides(colour = guide_legend(ncol = 1)) +
  coord_flip() +
  labs(
    y = "Number of eQTLs",
    x = "Method"
  ) + 
  facet_wrap(~cell_type, scales = "free_x") +
  guides(fill = "none") + 
  theme_jp_vgrid() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# Gene level analysis.

gene_plot_data <- bind_rows(
  quasar_gene_data |>
    filter(cell_type %in% c("Plasma", "B IN", "CD4 NC")) |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  jaxqtl_gene_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  apex_gene_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)),
  tensorqtl_gene_qvalue_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)) |>
    mutate(method = "tensorqtl")
) |>
  mutate(
    method = fct_reorder(factor(method_lookup[method]), n_sig),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

cor_data <- bind_rows(
  tensorqtl_gene_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)) |>
    mutate(method = "Perm, BH"),
  tensorqtl_gene_acat_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)) |>
    mutate(method = "ACAT, BH"),
  tensorqtl_gene_qvalue_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)) |>
    mutate(method = "Perm, q-value"),
  tensorqtl_gene_acat_qvalue_data |>
    summarise(n_sig = sum(n_sig), .by = c(method, cell_type)) |>
    mutate(method = "ACAT, q-value"),
) |>
  mutate(
    method = fct_reorder(factor(method), n_sig),
    cell_type = factor(cell_type, levels = c("Plasma", "B IN", "CD4 NC"))
  )

p1 <- cor_data |>
  ggplot(aes(method, n_sig, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = cell_type_cols) +
  guides(colour = guide_legend(ncol = 1)) +
  coord_flip() +
  labs(
    y = "Number of eGenes",
    x = ""
  ) +
  facet_wrap(~cell_type, scales = "free_x") +
  guides(fill = "none") +
  theme_jp_vgrid() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

acat_upset_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl_acat"
  ) |>
  rowwise() |>
  mutate(data = list(read_parquet(tensorqtl_file) |>
    summarise(pvalue = acat(pval_nominal), .by = phenotype_id)
  )) |>
  mutate(
    pvalue = list(data$pvalue),
    phenotype_id = list(data$phenotype_id)
  ) |>
  unnest(cols = c(pvalue, phenotype_id)) |>
  ungroup() |>
  mutate(
    pvalue_acat_bh = p.adjust(pvalue, method = "BH"),
    pvalue_acat_q = qvalue(pvalue)$qvalues,
    .by = cell_type
  ) |>
  select(phenotype_id, cell_type, pvalue_acat_bh, pvalue_acat_q)

perm_upset_data <- tibble(tensorqtl_file = to_r_vec(args[6])) |>
  mutate(
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  rowwise() |>
  mutate(pvalue = list(read_tsv(tensorqtl_file)$pval_beta)) |>
  mutate(phenotype_id = list(read_tsv(tensorqtl_file)$phenotype_id)) |>
  unnest(cols = c(pvalue, phenotype_id)) |>
  ungroup() |>
  mutate(
    pvalue_perm_bh = p.adjust(pvalue, method = "BH"),
    pvalue_perm_q = qvalue(pvalue)$qvalue,
    .by = cell_type
  ) |>
  select(phenotype_id, cell_type, pvalue_perm_bh, pvalue_perm_q)

upset_data <- acat_upset_data |>
  left_join(perm_upset_data, by = c("cell_type", "phenotype_id")) |>
  filter(cell_type == "CD4 NC") |>
  select(-cell_type) |>
  pivot_longer(
    cols = -phenotype_id,
    values_to = "pvalue",
    names_to = "method"
  ) |>
  filter(pvalue < 0.05) |>
  mutate(method = pvalue_method_lookup[method]) |>
  summarise(method = list(method), .by = phenotype_id)

p2 <- upset_data |>
  ggplot(aes(x = method)) +
  geom_bar(fill = "#BBBBBB") +
  scale_x_upset(n_intersections = 10) +
  labs(
    y = "Number of eGenes",
    x = "Combination of p-value methods"
  ) +
  theme_jp() +
  theme(axis.title = element_text(family = "Helvetica", size = 18, color = "#222222")) + 
  theme_combmatrix(
    combmatrix.label.text = element_text(family = "Helvetica", size = 18, color = "#222222")
  )

p <- p1 + p2 + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(
  "plot-cor-method-power.pdf",
  p,
  width = 12,
  height = 10
)

gene_power_plot <- gene_plot_data |>
  ggplot(aes(method, n_sig, fill = cell_type)) +
  geom_col() +
  scale_fill_manual(values = cell_type_cols) +
  guides(colour = guide_legend(ncol = 1)) +
  coord_flip() +
  labs(
    y = "Number of eGenes",
    x = ""
  ) + 
  facet_wrap(~cell_type, scales = "free_x") + 
  guides(fill = "none") + 
  theme_jp_vgrid() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

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
    filter(pvalue_adjust < 0.05) |>
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

total_data <- upset_plot_data |> 
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
  mutate(type = if_else(type == "all", "All models", "Count models only"))

wilcox.test(sum ~ type, data = total_data)

total_boxplot <- total_data |>
  ggplot(aes(factor(type), log10(sum))) +
  geom_boxplot() +
  coord_flip() + 
  labs(
    y = "log10(Total counts)",
    x = "Type of eGene"
  ) +
  theme_jp_vgrid()

power_plot <- (variant_power_plot + gene_power_plot) +
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    legend.direction = "vertical"
  )

ggsave(
  "plot-comparison-power.pdf", 
  power_plot,
  width = 12,
  height = 8
)

write_csv(gene_plot_data, "gene-power-data.csv")
write_csv(variant_plot_data, "variant-power-data.csv")

power_supp_plot <- upset_plot + total_boxplot +
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    legend.direction = "vertical"
  )

ggsave(
  "plot-supp-power.pdf",
  power_supp_plot,
  width = 12,
  height = 10
)

overall_variant_data <- quasar_variant_data |>
  filter(method %in% c("quasar-lm", "quasar-nb_glm-apl")) |>
  mutate(method = if_else(method == "quasar-lm", "LM", "NB-GLM (APL)")) |>
  summarise(n_sig = sum(n_sig), .by = c(method, cell_type))

p1 <- overall_variant_data |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig)) |>
  ggplot(aes(cell_type, n_sig, fill = method)) +
  geom_col(position = "dodge2") +
  scale_y_continuous(labels = label_number(scale_cut = cut_long_scale())) +
  scale_fill_manual(values = c("#99DDFF", "#44BB99")) +
  labs(
    x = "Cell type",
    y = "Number of eQTLs",
    fill = "Method"
  ) + 
  theme_bw()

overall_gene_data <- quasar_gene_data |>
  filter(method %in% c("quasar-lm", "quasar-nb_glm-apl")) |>
  mutate(method = if_else(method == "quasar-lm", "LM", "NB-GLM (APL)")) |>
  summarise(n_sig = sum(n_sig), .by = c(method, cell_type))

p2 <- overall_gene_data |>
  mutate(cell_type = fct_reorder(factor(cell_type), n_sig)) |>
  ggplot(aes(cell_type, n_sig, fill = method)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = c("#99DDFF", "#44BB99")) +
  labs(
    x = "Cell type",
    y = "Number of eGenes",
    fill = "Method"
  ) +
  theme_bw()

p <- (p1 / p2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    plot.tag = element_text(size = 18),
    legend.direction = "vertical",
    legend.position = "right"
  )

write_csv(overall_variant_data, "overall-variant-data.csv")
write_csv(overall_gene_data, "overall-gene-data.csv")

ggsave(
  "plot-overall-power.pdf", 
  p,
  width = 12,
  height = 10
)
