#!/bin/bash

# This script downloads the randolph data and extract the relevant component data.
cd data/randolph

wget -O randolph-data.tar.gz "https://zenodo.org/records/4273999/files/inputs.tar.gz?download=1"

tar -xvzf randolph-data.tar.gz inputs/3_eQTL_mapping/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt --strip-components=2
tar -xvzf randolph-data.tar.gz inputs/3_eQTL_mapping/genotypes.txt --strip-components=2
tar -xvzf randolph-data.tar.gz inputs/3_eQTL_mapping/GRCh38.92_gene_positions.txt --strip-components=2
tar -xvzf randolph-data.tar.gz inputs/1_calculate_pseudobulk/mergedAllCells_withCellTypeIdents_CLEAN.rds --strip-components=2
