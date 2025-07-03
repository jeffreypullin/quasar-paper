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

time_plot <- readRDS(args[[1]])
concordance_plot <- readRDS(args[[2]])

p <- concordance_plot + time_plot +
  plot_layout(guides = "collect", widths = c(1, 2.5)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 18))

ggsave(
  "figure-1.pdf",
  p,
  width = 12,
  height = 8
)
