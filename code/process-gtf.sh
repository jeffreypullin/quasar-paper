#! /bin/bash

python code/collapse-annotation.py \
 data/onek1k/Homo_sapiens.GRCh37.82.gtf.gz \
 data/onek1k/Homo_sapiens.GRCh37.82.genes.gtf \
 --collapse_only

python code/gtf-to-tss.py

zcat data/onek1k/Homo_sapiens.GRCh37.82.bed.gz | awk -F'\t' '$1 ~ /(^[1-9]$)|(^1[0-9]$)|(^2[012]$)/ {print $4, $1, $2, $3}' > data/onek1k/Homo_sapiens.GRCh37.82.bed
sed -i "1s/^/feature_id\tchrom\tstart\tend\n/" data/onek1k/Homo_sapiens.GRCh37.82.bed
