#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(arrow)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
quasar_files <- read_csv(
    args[[1]], 
    col_names = "quasar_file", 
    show_col_types = FALSE
)
tensorqtl_files <- read_csv(
    args[[2]], 
    col_names = "tensorqtl_file", 
    show_col_types = FALSE
)

qtl_output_data <- left_join(
  tensorqtl_files |>
    mutate(chr = paste0("chr", str_extract(tensorqtl_file, "[0-9]+(?=\\.parquet)"))),
  quasar_files |>
    mutate(chr = str_extract(quasar_file, "chr[0-9]+")), 
  by = "chr"
) |>
  relocate(chr, tensorqtl_file)

for (i in seq_len(nrow(qtl_output_data))) {
  
  quasar_data <- read_tsv(
    qtl_output_data$quasar_file[[i]], 
    show_col_types = FALSE
  )
  tensorqtl_data <- read_parquet(
    qtl_output_data$tensorqtl_file[[i]]
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
    geom_point()
    geom_abline()
  p2 <- plot_data |>  
    ggplot(aes(x = se, y = slope_se)) +
    geom_point()
    geom_abline()
  p3 <- plot_data |>  
    ggplot(aes(x = -log10(pvalue), y = -log10(pval_nominal))) +
    geom_point()
    geom_abline()
  p <- p1 + p2 + p3 
  ggsave("plot.pdf", p)
  break
}


#jaxqtl_data <- read_tsv(jaxqtl_path)
#quasar_data <- read_tsv(quasar_path)
#tensorqtl_data <- read_parquet(tensorqtl_path)

#jaxqtl_data |>
#    select(beta_converged, opt_status, slope, slope_se)

#plot_data <- quasar_data |>
#  filter(feature_id == "ENSG00000153575") |>
#  mutate(variant_id = paste0(chrom, ":", pos, alt, "-", ref)) |>
#  left_join(tensorqtl_data, by = join_by(feature_id == phenotype_id, variant_id))

#plot_data |>
#  ggplot(aes(beta, slope)) + 
#  geom_point() + 
#  geom_abline()

#plot_data |>
#  ggplot(aes(se, slope_se)) + 
#  geom_point()

#plot_data |>
#  ggplot(aes(-log10(pvalue), -log10(pval_nominal))) + 
#  geom_point()

#plot_data |>
#  select(variant_id, beta, slope)

