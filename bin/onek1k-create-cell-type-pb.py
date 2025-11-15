#! /usr/bin/env python3

import sys
import scanpy as sc
import pandas as pd

def flatten(xss):
    return [x for xs in xss for x in xs]

onek1k = sc.read_h5ad(sys.argv[1])
cell_label = sys.argv[2]
pb_type = sys.argv[3]

sc.pp.filter_genes(onek1k, min_cells=3)
sc.pp.filter_cells(onek1k, min_genes=100)

onek1k.layers["counts"] = onek1k.X.copy()
sc.pp.normalize_total(onek1k)
sc.pp.log1p(onek1k)

cell_type_subset = onek1k[onek1k.obs.cell_label == cell_label, :]

pbs = []
for indiv in cell_type_subset.obs.individual.unique():
    indiv_cell_subset = cell_type_subset[cell_type_subset.obs['individual'] == indiv]
    
    rep_adata = sc.AnnData(X = indiv_cell_subset.X.mean(axis = 0),
                           var = indiv_cell_subset.var[[]])
    rep_adata.layers['counts'] = indiv_cell_subset.layers['counts'].sum(axis = 0)
    rep_adata.obs_names = [indiv]
    rep_adata.obs['individual'] = indiv
    rep_adata.obs['sex'] = indiv_cell_subset.obs['sex'].iloc[0]
    rep_adata.obs['age'] = indiv_cell_subset.obs['age'].iloc[0]
    pbs.append(rep_adata)

pb = sc.concat(pbs)

# Compute PCA.
sc.tl.pca(pb)

# We only save the first 20 PCs.
pc_data = pb.obsm['X_pca'][0:, 0:20]

pc_df = pd.DataFrame(
  data=pc_data[0:,0:],
  index=pb.obs_names,
  columns=[f"PC_{i}" for i in range(1, 21)]
)

# Write out covariates.
pb_cov = pb.obs[['sex', 'age']]
pb_cov = pd.concat([pb_cov, pc_df], axis=1)
pb_cov.index.name = 'sample_id'
pb_cov.to_csv(f"onek1k-{cell_label}-cov.tsv", sep = "\t")

# Write out phenotype values.

# Filter to genes where >10% of counts are non-zero.
counts_matrix = pb.layers['counts']
nonzero_fracs = (counts_matrix > 0).mean(axis=0)
keep_genes = nonzero_fracs > 0.1
pb = pb[:, keep_genes]

gene_ids = pb.var.index.tolist()
indiv_ids = pb.obs['individual'].tolist()

first_line = indiv_ids
first_line.insert(0, 'feature_id')
file_name = f"onek1k-{cell_label}-pheno-{pb_type}.txt"

# Write out mean file.
if pb_type == "mean":
    with open(file_name, 'w') as fmean:
        fmean.write('\t'.join(first_line) + '\n') 
        for id in gene_ids:
            gene_res_list = flatten(pb[:, id].X.toarray().tolist())
            line = gene_res_list
            line.insert(0, id)
            fmean.write('\t'.join(map(str, line)) + '\n') 
# Write out sum file.
elif pb_type.startswith("sum"):
    with open(file_name, 'w') as fsum:
        fsum.write('\t'.join(first_line) + '\n') 
        for id in gene_ids:
            gene_res_list = flatten(pb[:, id].layers['counts'].toarray().tolist())
            line = gene_res_list
            line.insert(0, id)
            fsum.write('\t'.join(map(str, line)) + '\n') 

