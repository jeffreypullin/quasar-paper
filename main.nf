#!/usr/bin/env nextflow

// Randolph.
params.randolph_raw_cov_data = file("data/randolph/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt")
params.randolph_raw_feat_anno_data = file("data/randolph/GRCh38.92_gene_positions.txt")
params.randolph_raw_genotype_data = file("data/randolph/genotypes.txt")
params.randolph_raw_single_cell_data = file("data/randolph/mergedAllCells_withCellTypeIdents_CLEAN.rds")

// OneK1K.
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
    cell_type_subsets = ONEK1K_CREATE_CELLTYPE_PB(
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
    all_covs = ONEK1K_PROCESS_COVARIATES(cell_type_subsets, geno_pcs) 

    quasar_input = cell_type_subsets
      .map({ [it[0], it[2]] })
      .join(all_covs)
      .combine(bed_files)
      .first()
    
    RUN_QUASAR(
      quasar_input, 
      grm,
      params.onek1k_feature_annot 
    )

    // Randolph data.
    chrs = channel.from(1..22)
    RANDOLPH_PROCESS_COVARIATES(params.randolph_raw_cov_data)
    RANDOLPH_PROCESS_FEAT_ANNO_DATA(params.randolph_raw_feat_anno_data)
    randolph_geno_data = RANDOLPH_PROCESS_GENOTYPE_DATA(params.randolph_raw_genotype_data)
    //randolph_grm = RANDOLPH_CREATE_GRM(randolph_geno_data)
    RANDOLPH_SPLIT_GENOTYPE_DATA(chrs.merge(randolph_geno_data))
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
    output: tuple val(cell_type), path("onek1k-${cell_type}-cov.tsv"), path("onek1k-${cell_type}-pheno-mean.txt")

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
        tuple val(cell_type), val(covs), val(pheno)
        val geno_pcs 
    output: tuple val(cell_type), path("onek1k-${cell_type}-all-covs.tsv")
    
    script:
    """
    onek1k-process-covariates.R "$covs" "$geno_pcs" "$cell_type"
    """
}

process RUN_QUASAR {

    input: 
      tuple val(cell_type), val(pheno), val(covs), val(chr), val(bed)
      val grm 
      val feat_anno
    output: tuple val(chr), val(cell_type), path("todo")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    /home/jp2045/quasar/build/quasar -b $prefix \
      -f "$feat_anno" \
      -p "$pheno" \
      -c "$covs" \
      -g "$grm" \
      -m lmm \
      --verbose
    """
}




// Randolph data.
process RANDOLPH_PROCESS_COVARIATES {

    input: file raw_cov_data
    output: path("cov-data.tsv")

    script: 
    """
    sed -i '1s/^/sample_id,/' ${raw_cov_data}
    randolph-process-covariates.R ${raw_cov_data}
    """
}

process RANDOLPH_PROCESS_FEAT_ANNO_DATA {

    input: val raw_feat_anno_data
    output: path("feat-anno-data.tsv")

    script:
    """
    sed '1s/.*/ind\tfeature_id\tchrom\tstart\tend/' ${raw_feat_anno_data} > tmp-feat-anno-data.tsv
    randolph-process-feat-anno.R tmp-feat-anno-data.tsv
    """
}

process RANDOLPH_PROCESS_GENOTYPE_DATA {

    input: val randolph_genotype_data
    output: tuple path("randolph.bed"), path("randolph.map")

    // First convert the txt file to a ped file then into a bed.
    script:
    """
    module load plink/1.9
    sed '1s/^/snp_id\t/' ${randolph_genotype_data} > tmp-geno.txt
    randolph-process-genotype.R tmp-geno.txt
    plink --file randolph --make-bed --out randolph
    """
}

process RANDOLPH_CREATE_GRM {

    input: tuple val(randolph_bed), val(randolph_map)
    output: path("randolph_grm.tsv")

    script:
    """
    randolph-make-grm.sh ${randolph_bed.getParent().toString() + '/' + randolph_bed.getSimpleName()}
    randolph-process-grm.R
    """
}

process RANDOLPH_SPLIT_GENOTYPE_DATA {

    input: tuple val(chr), val(randolph_bed), val(randolph_map)
    output: path("randolph_${chr}.bed")

    script:
    """
    module load plink
    plink2 --bfile ${randolph_bed.getParent().toString() + '/' + randolph_bed.getSimpleName()} \
      --chr ${chr} \
      --make-bed \
      --out randolph_${chr}
    """
}

/*
process RANDOLPH_PROCESS_SINGLE_CELL_DATA {
    # TODO
}
*/

/*
process RUN_QUASAR {

    # TODO
    script:
    """
    ./quasar -b /home/jp2045/quasar/randolph-data/randolph_1 \
    -f /home/jp2045/quasar/randolph-data/randolph_feat_anno.tsv \
    -p /home/jp2045/quasar/randolph-data/b_corrected_exp.tsv  \
    -c /home/jp2045/quasar/randolph-data/randolph_cov.tsv \
    -g /home/jp2045/quasar/randolph-data/randolph_grm.tsv \
    -m lmm \
    --verbose
    """

}
*/
