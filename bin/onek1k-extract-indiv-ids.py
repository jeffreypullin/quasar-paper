#! /usr/bin/env python3

import sys
import scanpy as sc

onek1k = sc.read_h5ad(sys.argv[1])
indiv_ids = onek1k.obs.individual.unique().tolist()

with open(r'onek1k-indiv-ids.txt', 'w') as fp:
    for id in indiv_ids:
        fp.write("%s\n" % id)

