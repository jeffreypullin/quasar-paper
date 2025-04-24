#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(arrow)
  library(patchwork)
  library(tidyr)
})

# FIXME: Make this more reproducible.
source("/home/jp2045/quasar-paper/code/plot-utils.R")

args <- commandArgs(trailingOnly = TRUE)

quasar_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(chr = str_extract(quasar_file, "chr[0-9]+")) |>
  mutate(model = str_extract(quasar_file, "(?<=-)[^-]+(?=-cis)")) |>
  relocate(chr, quasar_file)
old_quasar_data <- quasar_data

tensorqtl_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    chr = paste0("chr", str_extract(tensorqtl_file, "[0-9]+(?=\\.parquet)"))
  ) |>
  relocate(chr, tensorqtl_file)

jaxqtl_data <- tibble(jaxqtl_file = to_r_vec(args[3])) |>
  mutate(chr = str_extract(jaxqtl_file, "chr[0-9]+")) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = "chr") |>
  relocate(chr, jaxqtl_file)

apex_data <- tibble(apex_file = to_r_vec(args[4])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\.cis)")
  ) |>
  relocate(chr, apex_file)

# Comparison of linear model methods.
lm_qtl_data <- left_join(
  quasar_data |>
    filter(model == "lm"),
  tensorqtl_data,
  by = "chr"
)

for (i in seq_len(nrow(lm_qtl_data))) {

  quasar_data <- read_tsv(
    lm_qtl_data$quasar_file[[i]],
    show_col_types = FALSE
  )
  tensorqtl_data <- read_parquet(
    lm_qtl_data$tensorqtl_file[[i]]
  )

  first_feature <- quasar_data$feature_id[[1]]

  quasar_filtered <- quasar_data |>
    filter(feature_id == first_feature) |>
    mutate(variant_id = paste0(chrom, ":", pos, alt, "-", ref))

  tensorqtl_filtered <- tensorqtl_data |>
    filter(phenotype_id == first_feature)

  plot_data <- quasar_filtered |>
    left_join(tensorqtl_filtered, by = "variant_id")

  p1 <- plot_data |>
    ggplot(aes(x = beta, y = slope)) +
    geom_point() +
    geom_abline()
  p2 <- plot_data |>
    ggplot(aes(x = se, y = slope_se)) +
    geom_point() +
    geom_abline()
  p3 <- plot_data |>
    ggplot(aes(x = -log10(pvalue), y = -log10(pval_nominal))) +
    geom_point() +
    geom_abline()

  p <- p1 + p2 + p3
  ggsave("plot-lm-concordance.pdf", p)
  break
}

# Comparison of Negative Binomial GLM model methods.
quasar_data <- old_quasar_data
glm_qtl_data <- left_join(
  quasar_data |>
    filter(model == "nb_glm"),
  jaxqtl_data,
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

    lm_fit <- lm(log_pvalue ~ log_pval_nominal, data = data)
    print(lm_fit)
    pvalue_beta_vec <- c(pvalue_beta_vec, coef(lm_fit)[[2]])
    gene_vec <- c(gene_vec, data$feature_id[[1]])

  #p1 <- plot_data |>
  #  ggplot(aes(x = beta, y = slope)) +
  #  geom_point() +
  #  geom_abline()
  #p2 <- plot_data |>
  #  ggplot(aes(x = se, y = slope_se)) +
  #geom_point() +
  #  geom_abline()
  #p3 <- plot_data |>
  #  ggplot(aes(x = -log10(pvalue), y = -log10(pval_nominal))) +
  #  geom_point() +
  #  geom_abline()
  }
}

plot_data <- tibble(
  pvalue_beta = pvalue_beta_vec,
  gene = gene_vec
)

p <- plot_data |>
  filter(!is.na(pvalue_beta)) |>
  filter(pvalue_beta > 0) |>
  ggplot(aes(pvalue_beta)) +
  geom_histogram(binwidth = 0.05) +
  coord_cartesian(xlim = c(0.5, 1.5))

ggsave("plot-glm-concordance.pdf", p)

# Comparison of linear mixed model methods.
quasar_data <- old_quasar_data
lmm_qtl_data <- left_join(
  quasar_data |>
    filter(model == "lmm"),
  apex_data,
  by = "chr"
)

for (i in seq_len(nrow(lmm_qtl_data))) {

  quasar_data <- read_tsv(
    lmm_qtl_data$quasar_file[[i]],
    show_col_types = FALSE
  )

  apex_data <- read_tsv(
    lmm_qtl_data$apex_file[[i]],
    show_col_types = FALSE
  )

  first_feature <- quasar_data$feature_id[[1]]

  quasar_filtered <- quasar_data |>
    filter(feature_id == first_feature) |>
    mutate(variant_id = paste0(chrom, ":", pos, alt, "-", ref))
  
  apex_filtered <- apex_data |>
    filter(gene == first_feature) |>
    mutate(apex_variant_id = paste0(`#chrom`, ":", pos, ref, "-", alt)) |>
    rename(apex_beta = beta, apex_se = se, apex_pval = pval)
  
  plot_data <- quasar_filtered |>
    left_join(apex_filtered, by = join_by(variant_id == apex_variant_id))
  
  p1 <- plot_data |>
    ggplot(aes(x = beta, y = apex_beta)) +
    geom_point() +
    geom_abline()
  p2 <- plot_data |>
    ggplot(aes(x = se, y = apex_se)) +
    geom_point() +
    geom_abline()
  p3 <- plot_data |>
    ggplot(aes(x = -log10(pvalue), y = -log10(apex_pval))) +
    geom_point() +
    geom_abline()

  p <- p1 + p2 + p3
  ggsave("plot-lmm-concordance.pdf", p)
  break
}
