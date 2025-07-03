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

quasar_files_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-cis)")),
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

pvalue_beta_vec <- numeric(0)
gene_vec <- character(0)
for (i in seq_len(nrow(lm_qtl_data))) {

  quasar_data <- fread(lm_qtl_data$quasar_file[[i]]) |>
     mutate(variant_id = paste0(chrom, ":", pos, alt, "-", ref))
  tensorqtl_data <- read_parquet(lm_qtl_data$tensorqtl_file[[i]])

  features <- unique(quasar_data$feature_id)
  for (j in seq_along(features)) {

    feature <- features[[j]]

    quasar_filtered <- quasar_data |>
      filter(feature_id == feature)

    tensorqtl_filtered <- tensorqtl_data |>
      filter(phenotype_id == feature)

     if (nrow(tensorqtl_filtered) == 0 || nrow(quasar_filtered) == 0) {
      next
    } 

    fit_data <- quasar_filtered |>
      left_join(tensorqtl_filtered, by = "variant_id") |>
      mutate(
        log_pvalue = -log10(pvalue),
        log_pval_nominal = -log10(pval_nominal)
      )

    if (any(is.na(fit_data$log_pvalue)) || any(is.na(fit_data$log_pval_nominal))) {
      next
    }
    if (nrow(fit_data) < 10) {
        next
    }

    lm_fit <- cor.test( ~ log_pvalue + log_pval_nominal, data = fit_data)
    
    pvalue_beta_vec <- c(pvalue_beta_vec, unname(lm_fit$estimate))
    gene_vec <- c(gene_vec, fit_data$feature_id[[1]])
  }
}

lm_plot_data <- tibble(
  pvalue_beta = pvalue_beta_vec,
  gene = gene_vec, 
  model = "LM"
)

# Comparison of Negative Binomial GLM model methods.
glm_qtl_data <- left_join(
  quasar_files_data |>
    filter(model == "nb_glm"),
  jaxqtl_files_data,
  by = "chr"
)

pvalue_beta_vec <- numeric(0)
gene_vec <- character(0)
for (i in seq_len(nrow(glm_qtl_data))) {

  quasar_data <- read_tsv(
    glm_qtl_data$quasar_file[[i]],
    show_col_types = FALSE
  )

  jaxqtl_data <- bind_rows(
    lapply(glm_qtl_data$jaxqtl_file[[i]], read_parquet)
  )

  features <- unique(quasar_data$feature_id)
  for (j in seq_along(features)) {

    feature <- features[[j]]

    quasar_filtered <- quasar_data |>
      filter(feature_id == feature) |>
      mutate(variant_id = paste0(chrom, ":", pos, alt, "-", ref))

    jaxqtl_filtered <- jaxqtl_data |>
      filter(phenotype_id == feature)

    if (nrow(jaxqtl_filtered) == 0 || nrow(quasar_filtered) == 0) {
      next
    }

    data <- quasar_filtered |>
      left_join(jaxqtl_filtered, by = join_by(variant_id == snp)) |>
      mutate(
        log_pvalue = -log10(pvalue),
        log_pval_nominal = -log10(pval_nominal)
      )
    
    if (any(is.na(data$log_pvalue)) || any(is.na(data$log_pval_nominal))) {
        next
    }
    if (nrow(data) < 10) {
        next
    }

    glm_fit <- cor.test(~ log_pvalue + log_pval_nominal, data = data)
    pvalue_beta_vec <- c(pvalue_beta_vec, unname(glm_fit$estimate))
    gene_vec <- c(gene_vec, data$feature_id[[1]])
  }
}

glm_plot_data <- tibble(
  pvalue_beta = pvalue_beta_vec,
  gene = gene_vec,
  model = "NB-GLM"
)

# Comparison of linear mixed model methods.
lmm_qtl_data <- left_join(
  quasar_files_data |>
    filter(model == "lmm"),
  apex_files_data,
  by = "chr"
)

pvalue_beta_vec <- numeric(0)
gene_vec <- character(0)
for (i in seq_len(nrow(lmm_qtl_data))) {

  quasar_data <- read_tsv(
    lmm_qtl_data$quasar_file[[i]],
    show_col_types = FALSE
  )

  apex_data <- read_tsv(
    lmm_qtl_data$apex_file[[i]],
    show_col_types = FALSE
  )
  
  features <- unique(quasar_data$feature_id)
  for (j in seq_along(features)) {

    feature <- features[[j]]

    quasar_filtered <- quasar_data |>
      filter(feature_id == feature) |>
      mutate(variant_id = paste0(chrom, ":", pos, alt, "-", ref))

    apex_filtered <- apex_data |>
      filter(gene == feature) |>
      mutate(apex_variant_id = paste0(`#chrom`, ":", pos, ref, "-", alt)) |>
      rename(apex_beta = beta, apex_se = se, apex_pval = pval)

    if (nrow(apex_filtered) == 0 || nrow(quasar_filtered) == 0) {
      next
    } 

    fit_data <- quasar_filtered |>
      left_join(apex_filtered, by = join_by(variant_id == apex_variant_id)) |>
      mutate(
        log_pvalue = -log10(pvalue),
        log_apex_pval = -log10(apex_pval)
      )
    
    if (any(is.na(fit_data$log_pvalue)) || any(is.na(fit_data$log_apex_pval))) {
        next
    }
    if (nrow(fit_data) < 10) {
        next
    }

    lmm_fit <- cor.test(~ log_pvalue + log_apex_pval, data = fit_data)
    
    pvalue_beta_vec <- c(pvalue_beta_vec, unname(lmm_fit$estimate))
    gene_vec <- c(gene_vec, fit_data$feature_id[[1]])
  }
}

lmm_plot_data <- tibble(
  pvalue_beta = pvalue_beta_vec,
  gene = gene_vec, 
  model = "LMM"
)

plot_data <- bind_rows(
  lm_plot_data, 
  glm_plot_data, 
  lmm_plot_data
)

p <- plot_data |>
  filter(!is.na(pvalue_beta)) |>
  filter(pvalue_beta < 2 & pvalue_beta > 0) |>
  mutate(model = factor(model, levels = c("LM", "NB-GLM", "LMM"))) |>
  ggplot(aes(factor(model), pvalue_beta)) +
  geom_boxplot() +
  coord_flip(ylim = c(0, 1)) +
  labs(
    y = "log10 p-value correlation",
    x = "Statistical model"
  ) + 
  theme_jp_vgrid()

saveRDS(p, "plot-concordance.rds")

ggsave(
  "plot-concordance.pdf",
  p,
  width = 12,
  height = 8
)