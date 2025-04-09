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

jaxqtl_data <- tibble(jaxqtl_file = to_r_vec(args[3])) |>
  mutate(chr = str_extract(jaxqtl_file, "chr[0-9]+")) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = "chr") |>
  relocate(chr, jaxqtl_file)

tensorqtl_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(
    chr = paste0("chr", str_extract(tensorqtl_file, "[0-9]+(?=\\.parquet)"))
  ) |>
  relocate(chr, tensorqtl_file)

quasar_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(chr = str_extract(quasar_file, "chr[0-9]+")) |>
  mutate(model = str_extract(quasar_file, "(?<=-)[^-]+(?=-cis)")) |>
  relocate(chr, quasar_file)
old_quasar_data <- quasar_data

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

# Comparison of Negative Binomial model methods.
quasar_data <- old_quasar_data
glm_qtl_data <- left_join(
  quasar_data |>
    filter(model == "glm"),
  jaxqtl_data,
  by = "chr"
)

for (i in seq_len(nrow(lm_qtl_data))) {

  quasar_data <- read_tsv(
    glm_qtl_data$quasar_file[[i]],
    show_col_types = FALSE
  )

  jaxqtl_data <- bind_rows(
    lapply(glm_qtl_data$jaxqtl_file[[i]], read_parquet)
  )

  read_parquet(glm_qtl_data$jaxqtl_file[[i]][[1]]) |>
    select(chrom, slope, converged, alpha) |>
    print()

  first_feature <- quasar_data$feature_id[[1]]

  quasar_filtered <- quasar_data |>
    filter(feature_id == first_feature) |>
    mutate(variant_id = paste0(chrom, ":", pos, alt, "-", ref))

  jaxqtl_filtered <- jaxqtl_data |>
    filter(phenotype_id == first_feature)

  plot_data <- quasar_filtered |>
    left_join(jaxqtl_filtered, by = join_by(variant_id == snp))

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
  ggsave("plot-glm-concordance.pdf", p)
  break
}
