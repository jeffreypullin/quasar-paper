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
  mutate(time = read_time(quasar_file)) |>
  ungroup() |>
  mutate(cell_type = "B IN") |>
  summarise(time = sum(time), .by = c(cell_type, model)) |>
  mutate(method = paste0("quasar-", model), type = "both") |>
  select(-model)

tensorqtl_cis_nominal_data <- tibble(tensorqtl_file = args[2]) |>
  rowwise() |>
  mutate(time = read_time(tensorqtl_file)) |>
  ungroup() |>
  mutate(cell_type = "B IN") |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "tensorqtl", type = "cis_nominal")

tensorqtl_cis_data <- tibble(tensorqtl_file = args[3]) |>
  rowwise() |>
  mutate(time = read_time(tensorqtl_file)) |>
  ungroup() |>
  mutate(cell_type = "B IN") |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "tensorqtl", type = "cis")

jaxqtl_cis_nominal_data <- tibble(jaxqtl_file = to_r_vec(args[4])) |>
  mutate(chr = str_extract(jaxqtl_file, "chr[0-9]+")) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = "chr") |>
  rowwise() |>
  mutate(time = sum(map_dbl(jaxqtl_file, read_time))) |>
  ungroup() |>
  mutate(cell_type = "B IN") |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "jaxqtl", type = "cis_nominal")

jaxqtl_cis_data <- tibble(jaxqtl_file = to_r_vec(args[5])) |>
  mutate(chr = str_extract(jaxqtl_file, "chr[0-9]+")) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = "chr") |>
  rowwise() |>
  mutate(time = sum(map_dbl(jaxqtl_file, read_time))) |>
  ungroup() |>
  mutate(cell_type = "B IN") |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "jaxqtl", type = "cis")

apex_data <- tibble(apex_file = to_r_vec(args[6])) |>
  mutate(chr = str_extract(apex_file, "chr[0-9]+(?=\\-time)")) |>
  rowwise() |>
  mutate(time = read_time(apex_file)) |>
  ungroup() |>
  mutate(cell_type = "B IN") |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "apex", type = "both")

plot_data <- bind_rows(
  apex_data,
  quasar_data,
  jaxqtl_cis_data,
  jaxqtl_cis_nominal_data,
  tensorqtl_cis_data,
  tensorqtl_cis_nominal_data
)

time_plot <- plot_data |>
  mutate(min = time / 60) |>
  ggplot(aes(method, min, fill = type)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~cell_type) +
  theme_bw()

ggsave("plot-time.pdf", time_plot)
