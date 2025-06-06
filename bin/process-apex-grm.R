#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
rel <- read.table(args[[1]], sep = "\t")
ids <- read.table(args[[2]], sep = "\t")

rel <- as_tibble(rel)
id_vec <- ids$V2
colnames(rel) <- id_vec
rel <- mutate(rel, id = id_vec) 

out <- rel |> pivot_longer(cols = -id)
names(out) <- c("#id1", "id2", "grm")

write.table(
    out, 
    file = "apex-grm.tsv",
    sep = "\t", 
    quote = FALSE,
    row.names = FALSE, 
    col.names = TRUE
)
