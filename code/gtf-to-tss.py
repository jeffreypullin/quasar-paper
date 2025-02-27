
import qtl.io
import pandas as pd

annotation_gtf="data/onek1k/Homo_sapiens.GRCh37.82.genes.gtf"
gtf=qtl.io.gtf_to_tss_bed(annotation_gtf, feature='gene')
gtf.to_csv("data/onek1k/Homo_sapiens.GRCh37.82.bed.gz",index=False, sep="\t")
