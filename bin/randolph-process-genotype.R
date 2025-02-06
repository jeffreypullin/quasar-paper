#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Append missing snp_id column name to start of header.
#system("sed -i '1s/^/snp_id\t/' genotypes.txt")

args <- commandArgs(trailingOnly = TRUE)
ped_path <- args[[1]]

raw_geno_data <- read_tsv(ped_path)

n_samples <- ncol(raw_geno_data) - 1
n_snps <- nrow(raw_geno_data)

snp_id <- raw_geno_data$snp_id
chr <- str_split_i(snp_id, "_", 1)
bp <- str_split_i(snp_id, "_", 2)
ref <- str_split_i(snp_id, "_", 3)
alt <- str_split_i(snp_id, "_", 4)

map <- data.frame(
  CHR = chr,
  SNP = snp_id,
  CM = 0,
  BP = bp
)

sample_id <- colnames(raw_geno_data)[2:ncol(raw_geno_data)]

ped_first6 <- data.frame(
  FID = paste0(sample_id, "_NI"),
  IID = paste0(sample_id, "_NI"),
  FatherID = 0,
  MotherID = 0,
  Sex = 1,
  Phenotype = 1
)

alleles <- data.frame(
  "0" = paste(ref, ref),
  "1" = paste(ref, alt),
  "2" = paste(alt, alt)
)

tmp <- vector("list", length = n_snps)
for (i in 1:n_snps) {
  print(i)
  allele_lookup <- c(
    "-9" = "0 0",
    "0" = paste(ref[[i]], ref[[i]]),
    "1" = paste(ref[[i]], alt[[i]]),
    "2" = paste(alt[[i]], alt[[i]])
  )
  tmp[[i]] <- allele_lookup[as.character(raw_geno_data[i, 2:(n_samples + 1)])]
}
tmp <- lapply(tmp, unname)
ped_snps <- do.call("rbind", tmp)

ped <- cbind(ped_first6, t(ped_snps))

file.remove("randolph.ped")
file.remove("randolph.map")

write.table(
  ped,
  "randolph.ped",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t"
)

write.table(
  map,
  "randolph.map",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE,
  sep = "\t"
)