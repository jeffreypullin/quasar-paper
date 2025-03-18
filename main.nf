#!/usr/bin/env nextflow

params.onek1k_raw_single_cell_data = file("data/onek1k/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz")
params.onek1k_supp_tables = file("data/onek1k/science.abf3041_tables_s6_to_s19/science.abf3041_tables_s6_to_s19.xlsx")
params.onek1k_feature_annot = file("data/onek1k/Homo_sapiens.GRCh37.82.bed")
params.onek1k_cell_types = [
    'CD4 NC', 'CD4 ET', 'CD4 SOX4', 'CD8 ET', 'CD8 NC', 'CD8 S100B', 'NK', 
    'NK R', 'Plasma', 'B Mem', 'B IN', 'Mono C', 'Mono NC', 'DC'
]
params.pb_types = ['mean', 'sum']

workflow {

    cell_types = channel.fromList(params.onek1k_cell_types)
    pb_types = channel.fromList(params.pb_types)
    pb_spec = cell_types
        .combine(pb_types)

    cell_type_pb = ONEK1K_CREATE_CELLTYPE_PB(
        pb_spec,
        params.onek1k_raw_single_cell_data,
        params.onek1k_supp_tables
    )
    vcf_files = channel.fromFilePairs(
        "data/onek1k/chr*.dose.filtered.R2_0.8.vcf.gz",
        size: 1, 
        flat: true
    )
    chrs = vcf_files.map{ it -> it[0] }
    ids = ONEK1K_EXTRACT_INDIV_IDS(params.onek1k_raw_single_cell_data)
    bed_files = ONEK1K_CONVERT_VCF_TO_BED(vcf_files, ids)
    all_bed = ONEK1K_CONCAT_GENOTYPE_DATA(bed_files.map({ it[1] }).collect())

    geno_pcs = ONEK1K_COMPUTE_GENOTYPE_PCS(all_bed)
    grm = ONEK1K_CREATE_GRM(all_bed)

    covs = ONEK1K_PROCESS_COVARIATES(cell_type_pb.filter({it -> it[1] == "mean"}), geno_pcs) 

    quasar_input = cell_type_pb 
      .filter({ x -> x[1] == "mean" })
      .map({it -> [it[0], it[1], it[3]]})
      .join(covs)
      .combine(bed_files)
    
    // Run quasar.
    quasar_out = RUN_QUASAR(
      quasar_input, 
      grm,
      params.onek1k_feature_annot 
    )
    quasar_out.filter({it -> it[1] == "CD4 NC"}).view()

    pheno_bed = ONEK1K_MAKE_PHENO_BED(cell_type_pb, params.onek1k_feature_annot)
    t_covs = ONEK1K_TRANSPOSE_COVS(covs)

    // Run TensorQTL.
    tensorqtl_input = pheno_bed
      .filter({ x -> x[1] == "mean" })
      .join(t_covs)
      .combine(all_bed)

    tensorqtl_out = RUN_TENSORQTL(tensorqtl_input)
    
    // Run jaxQTL.
    jaxqtl_spec = pheno_bed
      .filter({ x -> x[1] == "sum" })
      .combine(chrs)
      .map({ it -> [it[3], it[0], it[2]] })
    
    chunks = ONEK1K_MAKE_GENELIST_CHUNKS(jaxqtl_spec, 100)
    split_chunks = chunks
        .flatMap({ chr, cell_type, chunks ->
            chunks.collect { chunk ->
                tuple(chr, cell_type, chunk)
            }
        })
        .map({ it -> [it, it[2].getBaseName()].flatten() } )
        // FIXME: How do other cases occur?
        .filter(it -> it[3].startsWith("chunk"))
        .filter(it -> !it[2].startsWith("chunk"))

    split_chunks.filter(it -> it[2].startsWith("chunk")).view()

    jaxqtl_input = jaxqtl_spec
        .combine(split_chunks, by: [0, 1])
        .combine(bed_files, by: 0)
        .combine(chrs.combine(covs), by: [0, 1])
    
    RUN_JAXQTL(jaxqtl_input)
}
 
// OneK1K data.
process ONEK1K_EXTRACT_INDIV_IDS {
    conda "$projectDir/envs/scanpy.yaml"

    input: val raw_sc_data
    output: path("onek1k-indiv-ids.txt")

    script: 
    """
    onek1k-extract-indiv-ids.py $raw_sc_data
    """
}

process ONEK1K_CREATE_CELLTYPE_PB {
    conda "$projectDir/envs/scanpy.yaml"
    label "small"

    input:
        tuple val(cell_type), val(pb_type) 
        val raw_sc_data
        val supp_tables
    output: tuple val(cell_type), val(pb_type),
        path("onek1k-${cell_type}-cov.tsv"), 
        path("onek1k-${cell_type}-pheno-${pb_type}.txt")

    script:
    """
    onek1k-create-cell-type-pb.py $raw_sc_data $supp_tables "${cell_type}" $pb_type
    """
}

process ONEK1K_CONVERT_VCF_TO_BED {
    conda "$projectDir/envs/cli.yaml"
    label "tiny"

    input: 
        tuple val(chr), val(vcf)
        val sample_file
    output: tuple val(chr), path("onek1k-${chr}.bed")

    shell:
    '''
    bcftools view -S !{sample_file} !{vcf} > tmp.vcf
    plink2 --vcf tmp.vcf \
        --maf 0.05 --hwe 1e-6 \
        --set-all-var-ids '@:#$r-$a' \
        --make-bed --out onek1k-!{chr} --const-fid
    '''
}

process ONEK1K_CONCAT_GENOTYPE_DATA {
    conda "$projectDir/envs/cli.yaml"
    label "tiny"

    input: val bed_files
    output: path("onek1k-all.bed")

    script:
    """
    printf "%s\\n" $bed_files | tr -d '[],' | sort -V -t/ -k9 > all_bed_files.txt
    awk -F. '{print \$1".bed", \$1".bim", \$1".fam"}' all_bed_files.txt > all_files.txt
    head all_files.txt
    plink --merge-list all_files.txt --make-bed --out onek1k-all
    """
}

process ONEK1K_COMPUTE_GENOTYPE_PCS {
    conda "$projectDir/envs/cli.yaml"

    input: val bed
    output: path("onek1k_geno_pcs.txt")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    plink2 --bfile $prefix --maf 0.05 --hwe 1e-6 --make-bed --out ${prefix}_filtered --threads 2
    plink2 --bfile ${prefix}_filtered --indep-pairwise 50000 200 0.05 --out onek1k_pruned_variants --threads 2 --const-fid 
    plink2 --bfile ${prefix}_filtered --extract onek1k_pruned_variants.prune.in --make-bed --out onek1k_pruned --const-fid 
    plink2 --bfile onek1k_pruned --pca 10

    cat plink2.eigenvec | \
      sed '1s/IID/sample_id/' | \
      sed '1s/PC/geno_pc/g' | \
      cut -f2- > onek1k_geno_pcs.txt
    """
}

process ONEK1K_CREATE_GRM {
    conda "$projectDir/envs/cli.yaml"
    
    input: val bed
    output: path("onek1k-grm.tsv")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    plink2 --bfile $prefix --maf 0.05 --hwe 1e-6 --make-bed --out onek1k_filtered --threads 2
    plink2 --bfile onek1k_filtered --indep-pairwise 250 50 0.2 --out onek1k_pruning_info --threads 2
    plink2 --bfile onek1k_filtered --extract onek1k_pruning_info.prune.in --make-king square --out king_ibd_out --threads 2
    process-grm.R onek1k
    """
}

process ONEK1K_PROCESS_COVARIATES {

    input: 
        tuple val(cell_type), val(pb_type), val(covs), val(pheno)
        val geno_pcs 
    output: tuple val(cell_type), path("onek1k-${cell_type}-all-covs.tsv")
    
    script:
    """
    onek1k-process-covariates.R "$covs" "$geno_pcs" "$cell_type"
    """
}

process ONEK1K_MAKE_PHENO_BED {

    input:
        tuple val(cell_type), val(pb_type), val(covs), val(pheno)
        val feat_anno
    output: tuple val(cell_type), val(pb_type), path("onek1k-${cell_type}-pheno-${pb_type}.bed")

    script:
    """
    onek1k-make-pheno-bed.R "$pheno" "$feat_anno" "$cell_type" "$pb_type"
    """
}

process ONEK1K_TRANSPOSE_COVS {

    input: tuple val(cell_type), val(covs)
    output: tuple val(cell_type), path("onek1k-${cell_type}-t-all-covs.tsv")

    script:
    """
    onek1k-transpose-covs.R "$covs" "$cell_type"
    """
}

process RUN_TENSORQTL {
    conda "$projectDir/envs/tensorqtl.yaml"
    label "tiny"

    input:
      tuple val(cell_type), val(pb_type), val(pheno), val(covs), val(bed)
    output: tuple val(cell_type), path("onek1k-${cell_type}.cis_qtl.txt.gz") 

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    python3 -m tensorqtl "$prefix" "$pheno" "onek1k-${cell_type}" \
      --covariates "$covs" \
      --mode cis
    """
}

process ONEK1K_MAKE_GENELIST_CHUNKS {
    input: 
        tuple val(chr), val(cell_type), val(pheno_bed)
        val split_n
    output: tuple val(chr), val(cell_type), path("chunk-*.tsv")

    script:
    """
    onek1k-make-genelist-chunks.R "$pheno_bed" $chr "$cell_type" $split_n
    """
}

process RUN_JAXQTL {
    conda "$projectDir/envs/jaxqtl.yaml"
    label "tiny"

    input: tuple val(chr), val(cell_type), val(pheno), val(chunk), val(chunk_id), val(bed), val(covs)
    output: tuple val(chr), val(cell_type), path("jaxqtl-${cell_type}_${chr}_${chunk_id}.cis_score.tsv.gz")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    jaxqtl \
      --geno "$prefix" \
      --covar "$covs" \
      --pheno "$pheno" \
      --model "NB" \
      --mode "cis" \
      --window 500000 \
      --genelist $chunk \
      --test-method "score" \
      --nperm 1000 \
      --addpc 0 \
      --standardize \
      --out "jaxqtl-${cell_type}_${chr}_${chunk_id}"
    """
}

process RUN_QUASAR {
    label "tiny"

    input: 
      tuple val(cell_type), val(pb_type), val(pheno), val(covs), val(chr), val(bed)
      val grm 
      val feat_anno
    output: tuple val(chr), val(cell_type),
        path("${chr}-${cell_type}-lm-cis-region.txt.gz"), 
        path("${chr}-${cell_type}-lm-cis-variant.txt.gz")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    /home/jp2045/quasar/build/quasar -b $prefix \
      -f "$feat_anno" \
      -p "$pheno" \
      -c "$covs" \
      -g "$grm" \
      -o "${chr}-${cell_type}-lm" \
      -m lm \
      --verbose
    gzip "${chr}-${cell_type}-lm-cis-variant.txt"
    gzip "${chr}-${cell_type}-lm-cis-region.txt"
    """
}
