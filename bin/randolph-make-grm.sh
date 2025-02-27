#!/usr/bin/env bash

module load plink

plink2 --bfile $1 --maf 0.05 --hwe 1e-6 --make-bed --out randolph_filtered --threads 2
plink2 --bfile randolph_filtered --indep-pairwise 250 50 0.2 --out randolph_pruning_info --threads 2
plink2 --bfile randolph_filtered --extract randolph_pruning_info.prune.in --make-king square --out king_ibd_out --threads 2
