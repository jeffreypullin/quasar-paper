#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(arrow)
  library(patchwork)
  library(tidyr)
  library(glue)
  library(data.table)
})

# FIXME: Make this more reproducible.
source("/home/jp2045/quasar-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

set.seed(20251107)

quasar_files_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-quasar)")),
  ) |>
  filter(cell_type == "B IN")

tensorqtl_files_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    chr = paste0("chr", str_extract(tensorqtl_file, "(?<=\\.)[0-9]+(?=\\.parquet)")),
    cell_type = str_extract(tensorqtl_file, "(?<=onek1k-).*?(?=\\.cis)"),
    method = "tensorqtl"
  ) |>
  filter(cell_type == "B IN")

jaxqtl_files_data <- tibble(jaxqtl_file = to_r_vec(args[3])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=jaxqtl-).*?(?=-chr)"),
  ) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = c("chr", "cell_type")) |>
  filter(cell_type == "B IN")

apex_files_data <- tibble(apex_file = to_r_vec(args[4])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)"),
    cell_type = str_extract(apex_file, "(?<=apex-).*?(?=-chr)"),
    method = "apex"
  ) |>
  filter(cell_type == "B IN")

# Comparison of linear model methods.
lm_qtl_data <- left_join(
  quasar_files_data |>
    filter(model == "lm"),
  tensorqtl_files_data,
  by = "chr"
)

lm_plot_data_list <- list()
for (i in seq_len(nrow(lm_qtl_data))) {

  quasar_data <- fread(lm_qtl_data$quasar_file[[i]])
  tensorqtl_data <- read_parquet(lm_qtl_data$tensorqtl_file[[i]])

  joined_data <- inner_join(
    quasar_data,
    tensorqtl_data,
    by = c("snp_id" = "variant_id", "feature_id" = "phenotype_id")
  ) |>
    mutate(
      z_quasar = beta / se,
      z_tensorqtl = slope / slope_se
    ) |>
    slice_sample(n = 100)

  lm_plot_data_list[[i]] <- joined_data
}

lm_plot_data <- bind_rows(!!!lm_plot_data_list) |>
  as_tibble()

# Comparison of Negative Binomial GLM model methods.
glm_qtl_data <- left_join(
  quasar_files_data |>
    filter(model == "nb_glm"),
  jaxqtl_files_data,
  by = "chr"
)

glm_plot_data_list <- list()
for (i in seq_len(nrow(glm_qtl_data))) {

  quasar_data <- fread(lm_qtl_data$quasar_file[[i]])
  jaxqtl_data <- bind_rows(
    lapply(glm_qtl_data$jaxqtl_file[[i]], read_parquet)
  )

  joined_data <- inner_join(
    quasar_data,
    jaxqtl_data,
    by = c("snp_id" = "snp", "feature_id" = "phenotype_id")
  ) |>
    mutate(
      z_quasar = beta / se,
      z_jaxqtl = slope / slope_se
    ) |>
    slice_sample(n = 100)

  glm_plot_data_list[[i]] <- joined_data
}

glm_plot_data <- bind_rows(!!!glm_plot_data_list) |>
  as_tibble()

# Comparison of linear mixed model methods.
lmm_qtl_data <- left_join(
  quasar_files_data |>
    filter(model == "lmm"),
  apex_files_data,
  by = "chr"
)

lmm_plot_data_list <- list()
for (i in seq_len(nrow(lmm_qtl_data))) {

  quasar_data <- fread(lm_qtl_data$quasar_file[[i]])
  apex_data <- read_tsv(
    lmm_qtl_data$apex_file[[i]],
    show_col_types = FALSE
  ) |>
  mutate(apex_variant_id = paste0(`#chrom`, ":", pos, ref, "-", alt)) |>
  rename(apex_beta = beta, apex_se = se, apex_pval = pval)
  
  joined_data <- inner_join(
    quasar_data,
    apex_data,
    by = c("snp_id" = "apex_variant_id", "feature_id" = "gene")
  ) |>
    mutate(
      z_quasar = beta / se,
      z_apex = apex_beta / apex_se
    ) |>
    slice_sample(n = 100)

  lmm_plot_data_list[[i]] <- joined_data
}

lmm_plot_data <- bind_rows(!!!lmm_plot_data_list) |>
  as_tibble()

plot_data <- bind_rows(
  lm_plot_data |>
    mutate(model = "LM"),
  glm_plot_data |>
    mutate(model = "NB-GLM"),
  lmm_plot_data |>
    mutate(model = "LMM")
) |>
  select(model, z_quasar, z_apex, z_tensorqtl, z_jaxqtl) |>
  pivot_longer(cols = -c(model, z_quasar))

p <- plot_data |>
  mutate(model = factor(model, levels = c("LM", "NB-GLM", "LMM"))) |>
  ggplot(aes(z_quasar, value)) +
  geom_point() +
  geom_abline() +
  labs(
    x = "quasar z-score",
    y = "Comparison method z-score"
  ) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(-10, 10)) +
  facet_wrap(~model, ncol = 1) +
  theme_jp() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

saveRDS(p, "plot-concordance.rds")

ggsave(
  "plot-concordance.pdf",
  p,
  width = 12,
  height = 8
)