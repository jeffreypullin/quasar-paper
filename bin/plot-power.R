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

args <- commandArgs(trailingOnly = TRUE)

to_r_vec <- function(str) {
    str <- str_replace(str, "\\[", "\\(")
    str <- str_replace(str, "\\]", "\\)")
    str <- str_replace_all(str, "([^,\\(\\)]+)", "'\\1'")
    str <- paste0("c", str)
    out <- eval(parse(text = str))
    out <- str_trim(out)
    out
}

qtl_output_data <- left_join(
  tibble(tensorqtl_file = to_r_vec(args[2])) |>
    mutate(chr = paste0("chr", str_extract(tensorqtl_file, "[0-9]+(?=\\.parquet)"))),
   tibble(quasar_file = to_r_vec(args[1])) |>
    mutate(chr = str_extract(quasar_file, "chr[0-9]+")), 
  by = "chr"
) |>
  relocate(chr, tensorqtl_file)


n_sig_quasar <- numeric()
n_sig_tensorqtl <- numeric()
chr <- character()
for (i in seq_len(nrow(qtl_output_data))) {
  
  quasar_data <- read_tsv(
    qtl_output_data$quasar_file[[i]], 
    show_col_types = FALSE
  )
  tensorqtl_data <- read_parquet(
    qtl_output_data$tensorqtl_file[[i]]
  )

  n_sig_quasar[[i]] <- sum(quasar_data$pvalue < 1e-6)
  n_sig_tensorqtl[[i]] <- sum(tensorqtl_data$pval_nominal < 1e-6)
  chr[[i]] <- qtl_output_data$chr[[i]]
}

p <- tibble(chr, n_sig_quasar, n_sig_tensorqtl) |>
  pivot_longer(cols = -chr, names_to = "method", values_to = "n_sig") |>
  ggplot(aes(method, n_sig)) + 
  geom_col() + 
  facet_wrap(~chr, scales = "free_y")

ggsave("plot-power.pdf", p)