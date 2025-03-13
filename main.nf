#!/usr/bin/env nextflow

params.onek1k_raw_single_cell_data = file("data/onek1k/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz")
params.onek1k_supp_tables = file("data/onek1k/science.abf3041_tables_s6_to_s19/science.abf3041_tables_s6_to_s19.xlsx")
params.onek1k_feature_annot = file("data/onek1k/Homo_sapiens.GRCh37.82.bed")
params.onek1k_cell_types = [
    'CD4 NC', 'CD4 ET', 'CD4 SOX4', 'CD8 ET', 'CD8 NC', 'CD8 S100B', 'NK', 
    'NK R', 'Plasma', 'B Mem', 'B IN', 'Mono C', 'Mono NC', 'DC'
]

workflow {

    // OneK1K data.
    cell_types = channel.fromList(params.onek1k_cell_types)
    cell_type_pb = ONEK1K_CREATE_CELLTYPE_PB(
        params.onek1k_raw_single_cell_data,
        params.onek1k_supp_tables,
        cell_types
    )
    ids = ONEK1K_EXTRACT_INDIV_IDS(params.onek1k_raw_single_cell_data)
    vcf_files = channel.fromFilePairs(
        "data/onek1k/chr*.dose.filtered.R2_0.8.vcf.gz",
        size: 1, 
        flat: true
    )
    chrs = vcf_files.map({ it.first() })
    bed_files = ONEK1K_CONVERT_VCF_TO_BED(vcf_files, ids)
    all_bed = ONEK1K_CONCAT_GENOTYPE_DATA(bed_files.map({ it[1] }).collect())
    geno_pcs = ONEK1K_COMPUTE_GENOTYPE_PCS(all_bed)
    grm = ONEK1K_CREATE_GRM(all_bed)
    all_covs = ONEK1K_PROCESS_COVARIATES(cell_type_pb, geno_pcs) 

    linear_models = channel.of('lm', 'lmm')
    generalised_models = channel.of('glm', 'glmm')

    mean_pb = cell_type_pb
      .map({ [it[0], it[2]] })
      .join(all_covs)
      .combine(linear_models)
    sum_pb = cell_type_pb
      .map({ [it[0], it[3]] })
      .join(all_covs)
      .combine(generalised_models)

    quasar_input = mean_pb
      .concat(sum_pb)
      .combine(bed_files)

    RUN_QUASAR(
      quasar_input, 
      grm,
      params.onek1k_feature_annot 
    )
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
        val raw_sc_data
        val supp_tables
        val cell_type
    output: tuple val(cell_type), 
        path("onek1k-${cell_type}-cov.tsv"), 
        path("onek1k-${cell_type}-pheno-mean.txt"),
        path("onek1k-${cell_type}-pheno-sum.txt")

    script:
    """
    onek1k-create-cell-type-pb.py $raw_sc_data $supp_tables "${cell_type}"
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
        tuple val(cell_type), val(covs), val(mean_pheno), val(sum_pheno)
        val geno_pcs 
    output: tuple val(cell_type), path("onek1k-${cell_type}-all-covs.tsv")
    
    script:
    """
    onek1k-process-covariates.R "$covs" "$geno_pcs" "$cell_type"
    """
}

process RUN_QUASAR {
    label "tiny"
    cache false

    input: 
      tuple val(cell_type), val(pheno), val(covs), val(model), val(chr), val(bed)
      val grm 
      val feat_anno
    output: tuple val(chr), val(cell_type), val(model),
        path("${chr}-${cell_type}-${model}-cis-region.txt.gz"), 
        path("${chr}-${cell_type}-${model}-cis-variant.txt.gz")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    /home/jp2045/quasar/build/quasar -b $prefix \
      -f "$feat_anno" \
      -p "$pheno" \
      -c "$covs" \
      -g "$grm" \
      -o "${chr}-${cell_type}-${model}" \
      -m $model \
      --verbose
    gzip "${chr}-${cell_type}-${model}-cis-variant.txt"
    gzip "${chr}-${cell_type}-${model}-cis-region.txt"
    """
}
