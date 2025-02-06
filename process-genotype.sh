#! /bin/bash

module load plink
export UNZIP_DISABLE_ZIPBOMB_DETECTION=TRUE

cd data
unzip -o ONEK1K_genotype_data_VCF.zip
plink2 --vcf OneK1K_AllChr.vcf --recode --make-bed --out OneK1K_AllChr --threads 2

