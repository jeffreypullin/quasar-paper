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
  library(glue)
  library(forcats)
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
  mutate(
    chr = str_extract(quasar_file, "chr[0-9]+"),
    cell_type = str_extract(quasar_file, "chr[0-9]+-(.*?)-", group = 1),
    model = str_extract(quasar_file, glue("(?<={cell_type}-).*?(?=-time)")),
  ) |> 
  rowwise() |>
  mutate(time = read_time(quasar_file)) |>
  ungroup() |>
  summarise(time = sum(time), .by = c(cell_type, model)) |>
  mutate(method = paste0("quasar-", model)) |>
  select(-model)

tensorqtl_cis_nominal_data <- tibble(tensorqtl_file = to_r_vec(args[2])) |>
  mutate(cell_type = str_extract(tensorqtl_file, "(?<=/)([^/]+)(?=-time)")) |>
  rowwise() |>
  mutate(time = read_time(tensorqtl_file)) |>
  ungroup() |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "tensorqtl_cis_nominal")

tensorqtl_cis_data <- tibble(tensorqtl_file = to_r_vec(args[3])) |>
  mutate(cell_type = str_extract(tensorqtl_file, "(?<=/)([^/]+)(?=-time)")) |>
  rowwise() |>
  mutate(time = read_time(tensorqtl_file)) |>
  ungroup() |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "tensorqtl_cis")

jaxqtl_cis_nominal_data <- tibble(jaxqtl_file = to_r_vec(args[4])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file,  "(?<=/)([^/]+)(?=-chr)"),
  ) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = c("chr", "cell_type")) |>
  rowwise() |>
  mutate(time = sum(map_dbl(jaxqtl_file, read_time))) |>
  ungroup() |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "jaxqtl_cis_nominal")

jaxqtl_cis_data <- tibble(jaxqtl_file = to_r_vec(args[5])) |>
  mutate(
    chr = str_extract(jaxqtl_file, "chr[0-9]+"),
    cell_type = str_extract(jaxqtl_file, "(?<=/)([^/]+)(?=-chr)"),
  ) |>
  summarise(jaxqtl_file = list(jaxqtl_file), .by = c("chr", "cell_type")) |>
  rowwise() |>
  mutate(time = sum(map_dbl(jaxqtl_file, read_time))) |>
  ungroup() |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "jaxqtl_cis")

apex_data <- tibble(apex_file = to_r_vec(args[6])) |>
  mutate(
    chr = str_extract(apex_file, "chr[0-9]+(?=\\-time)"),
    cell_type = str_extract(apex_file,  "(?<=/)([^/]+)(?=-chr)")
  ) |>
  rowwise() |>
  mutate(time = read_time(apex_file)) |>
  ungroup() |>
  summarise(time = sum(time), .by = cell_type) |>
  mutate(method = "apex")

model_type_lookup <- c(
  "jaxqtl_cis" = "GLM",
  "jaxqtl_cis_nominal" = "GLM",
  "tensorqtl_cis" = "LM",
  "tensorqtl_cis_nominal" = "LM",
  "apex" = "LMM",
  "quasar-p_glm" = "GLM",
  "quasar-nb_glm" = "GLM",
  "quasar-nb_glm-apl" = "GLM",
  "quasar-lm" = "LM",
  "quasar-lmm" = "LMM",
  "quasar-p_glmm" = "GLMM",
  "quasar-nb_glmm" = "GLMM"
)

plot_data <- bind_rows(
  apex_data,
  quasar_data,
  jaxqtl_cis_data,
  jaxqtl_cis_nominal_data,
  tensorqtl_cis_data,
  tensorqtl_cis_nominal_data
) |>
  mutate(
    model_type = model_type_lookup[method],
    method = fct_reorder(factor(method_lookup[method]), time, .desc = TRUE),
    cell_type = factor(cell_type, levels = c("CD4 NC", "B IN", "Plasma"))
  )

write_csv(plot_data, "time-data.csv")

time_plot <- plot_data |>
  mutate(min = time / 60) |>
  ggplot(aes(method, min, fill = cell_type)) +
  geom_col(position = "dodge2") +
  coord_flip() +
  facet_wrap(~model_type, scales = "free") + 
  scale_fill_manual(values = cell_type_cols) +
  labs(
    y = "Time (min)",
    fill = "Cell type",
    x = "Method"
  ) + 
  theme_jp_vgrid()

ggsave(
  "plot-time.pdf", 
  time_plot,
  width = 10,
  height = 6
)
