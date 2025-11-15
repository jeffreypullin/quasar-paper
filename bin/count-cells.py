#! /usr/bin/env python3

import sys
import scanpy as sc

adata = sc.read_h5ad(sys.argv[1])

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_genes=100)

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

exclude_labels = ['Erythrocytes', 'Platelets']
adata = adata[~adata.obs['cell_label'].isin(exclude_labels)].copy()

cell_counts = adata.obs['cell_label'].value_counts().sort_index()
cell_counts_df = cell_counts.rename_axis('cell_type').reset_index(name='n_cells')
cell_counts_df.to_csv("n-cells.tsv", sep="\t", index=False)
