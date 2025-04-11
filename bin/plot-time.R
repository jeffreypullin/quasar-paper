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


read_time <- function(path) {
  data <- read.delim(path)
  raw_str <- colnames(data)
  str <- str_extract(raw_str, "(?<=real\\.).*")
  as.numeric(str)
}

quasar_data <- tibble(quasar_file = to_r_vec(args[1])) |>
  mutate(chr = str_extract(quasar_file, "chr[0-9]+")) |>
  mutate(model = str_extract(quasar_file, "(?<=-)[^-]+(?=-time)")) |>
  rowwise() |>
  mutate(time = read_time(quasar_file))

p <- quasar_data |>
  ggplot(aes(model, time)) +
  geom_col() +
  facet_wrap(~chr)

ggsave("plot-time.pdf", p)
