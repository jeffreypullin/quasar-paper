#!/usr/bin/env Rscript

library(readr)
library(Matrix)
library(vcfR)

args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[[1]]
chr <- args[[2]]

vcf <- read.vcfR(vcf_file, verbose = FALSE)
col_names <- colnames(vcf@gt)

n <- length(col_names)

perm_ind <- sample(n - 1)
vcf@gt[, 2:n] <- vcf@gt[, perm_ind + 1]

write.vcf(vcf, paste0("permute-", chr, ".vcf.gz"))
