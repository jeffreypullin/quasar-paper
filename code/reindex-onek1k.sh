#! /bin/bash

for chr in {1..22}; do
    bcftools index --tbi --force data/onek1k/chr${chr}.dose.filtered.R2_0.8.vcf.gz
done
