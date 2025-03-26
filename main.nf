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
    filt_vcf_files = ONEK1K_FILTER_VCF(vcf_files, ids)
    bed_files = ONEK1K_CONVERT_VCF_TO_BED(filt_vcf_files)
    all_bed = ONEK1K_CONCAT_BED_FILES(bed_files.map({ it[1] }).collect())

    geno_pcs = ONEK1K_COMPUTE_GENOTYPE_PCS(all_bed)
    grm = ONEK1K_CREATE_GRM(all_bed)

    covs = ONEK1K_PROCESS_COVARIATES(cell_type_pb.filter({it -> it[1] == "mean"}), geno_pcs) 

    quasar_input = cell_type_pb 
      .filter({ x -> x[1] == "mean" })
      .map({it -> [it[0], it[1], it[3]]})
      .join(covs)
      .combine(bed_files)
      .filter( { it -> it[0] == "B IN" })
    
    // Run quasar.
    quasar_out = RUN_QUASAR(
      quasar_input, 
      grm,
      params.onek1k_feature_annot 
    )

    pheno_bed = ONEK1K_MAKE_PHENO_BED(cell_type_pb, params.onek1k_feature_annot)
    t_covs = ONEK1K_TRANSPOSE_COVS(covs)

    // Run TensorQTL.
    tensorqtl_input = pheno_bed
      .filter({ x -> x[1] == "mean" })
      .join(t_covs)
      .combine(all_bed)
      .filter( { it -> it[0] == "B IN"})

    //RUN_TENSORQTL_CIS(tensorqtl_input)
    tensorqtl_cis_nominal_out = RUN_TENSORQTL_CIS_NOMINAL(tensorqtl_input)
    
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

    jaxqtl_input = jaxqtl_spec
        .combine(split_chunks, by: [0, 1])
        .combine(bed_files, by: 0)
        .combine(chrs.combine(covs), by: [0, 1])
        .filter( { it -> it[1] == "B IN"})
        .map( { it -> [it, it[0].substring(3).toInteger()].flatten()})

    RUN_JAXQTL_CIS(jaxqtl_input)
    RUN_JAXQTL_CIS_NOMINAL(jaxqtl_input)

    // Run apex.
    pheno_bed_gz = COMPRESS_AND_INDEX_BED(pheno_bed)
    sparse_grm = TO_SPARSE_GRM_FORMAT(grm)

    apex_input = pheno_bed_gz
        .filter({ x -> x[1] == "mean" })
        .combine(chrs)
        .map({ it -> [it[4], it[0], it[1], it[2], it[3]] })
        .combine(filt_vcf_files, by: 0)
        .combine(chrs.combine(t_covs), by: [0, 1])
        .combine(sparse_grm)
        .filter( { it -> it[1] == "B IN"})
    
    RUN_APEX(apex_input) 

    grouped_quasar_variant = quasar_out
        .map( { it -> [it[1], it[0], it[3]]})
        .groupTuple()

    PLOT_CONCORDANCE(grouped_quasar_variant, tensorqtl_cis_nominal_out)
    PLOT_POWER(grouped_quasar_variant, tensorqtl_cis_nominal_out)
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

process ONEK1K_FILTER_VCF {

    conda "$projectDir/envs/cli.yaml"

    input: 
        tuple val(chr), val(vcf)
        val sample_file
    output: tuple val(chr), path("filtered-${chr}.vcf.gz"), path("filtered-${chr}.vcf.gz.tbi")

    shell:
    """
    bcftools view -S $sample_file $vcf > tmp.vcf
    vcftools --vcf tmp.vcf \
        --maf 0.05 \
        --hwe 1e-6 \
        --recode \
        --recode-INFO-all \
        --stdout \
        --stdout | bgzip -c > "filtered-${chr}.vcf.gz"
    tabix -p vcf "filtered-${chr}.vcf.gz"
    """
}

process ONEK1K_CONVERT_VCF_TO_BED {
    conda "$projectDir/envs/cli.yaml"
    label "tiny"

    input: tuple val(chr), val(vcf), val(vcf_tbi)
    output: tuple val(chr), path("onek1k-${chr}.bed")

    shell:
    '''
    plink2 --vcf !{vcf} \
        --set-all-var-ids '@:#$r-$a' \
        --make-bed --out onek1k-!{chr} --const-fid
    '''
}

process ONEK1K_CONCAT_BED_FILES {
    conda "$projectDir/envs/cli.yaml"
    label "tiny"

    input: val bed_files
    output: path("onek1k-all.bed")

    script:
    """
    printf "%s\\n" $bed_files | tr -d '[],' | sort -V -t/ -k9 > all_bed_files.txt
    awk -F. '{print \$1".bed", \$1".bim", \$1".fam"}' all_bed_files.txt > all_files.txt
    plink --keep-allele-order --merge-list all_files.txt --make-bed --out onek1k-all
    """
}

process ONEK1K_COMPUTE_GENOTYPE_PCS {
    conda "$projectDir/envs/cli.yaml"

    input: val bed
    output: path("onek1k_geno_pcs.txt")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    plink2 --bfile $prefix --indep-pairwise 50000 200 0.05 --out onek1k_pruned_variants --threads 2 --const-fid 
    plink2 --bfile $prefix --extract onek1k_pruned_variants.prune.in --make-bed --out onek1k_pruned --const-fid 
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
    plink2 --bfile $prefix --indep-pairwise 250 50 0.2 --out onek1k_pruning_info --threads 2
    plink2 --bfile $prefix --extract onek1k_pruning_info.prune.in --make-king square --out king_ibd_out --threads 2
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

process RUN_TENSORQTL_CIS {
    conda "$projectDir/envs/tensorqtl.yaml"
    label "tensorqtl"

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

process RUN_TENSORQTL_CIS_NOMINAL {
    conda "$projectDir/envs/tensorqtl.yaml"
    label "tensorqtl"

    input:
      tuple val(cell_type), val(pb_type), val(pheno), val(covs), val(bed)
    output: tuple val(cell_type), path("onek1k-${cell_type}.cis_qtl_pairs.*.parquet") 

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    python3 -m tensorqtl "$prefix" "$pheno" "onek1k-${cell_type}" \
      --covariates "$covs" \
      --mode cis_nominal
    """
}

process ONEK1K_MAKE_GENELIST_CHUNKS {
    input: 
        tuple val(chr), val(cell_type), val(pheno_bed)
        val split_n
    output: tuple val(chr), val(cell_type), path("chunk-*.tsv", arity: "1..*")

    script:
    """
    onek1k-make-genelist-chunks.R "$pheno_bed" $chr "$cell_type" $split_n
    """
}

process RUN_JAXQTL_CIS {
    conda "$projectDir/envs/jaxqtl.yaml"
    label "jaxqtl"

    input: tuple val(chr), val(cell_type), val(pheno), val(chunk), val(chunk_id), val(bed), val(covs), val(chr_num)
    output: tuple val(chr), val(cell_type), path("jaxqtl-${cell_type}-${chr}-${chunk_id}.cis_score.tsv.gz")

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
      --out "jaxqtl-${cell_type}-${chr}-${chunk_id}"
    """
}

process RUN_JAXQTL_CIS_NOMINAL {
    conda "$projectDir/envs/jaxqtl.yaml"
    label "jaxqtl"

    input: tuple val(chr), val(cell_type), val(pheno), val(chunk), val(chunk_id), val(bed), val(covs), val(chr_num)
    output: tuple val(chr), val(cell_type), path("jaxqtl-${cell_type}-${chr}-${chunk_id}.cis_qtl_pairs.${chr_num}.score.parquet")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    jaxqtl \
      --geno "$prefix" \
      --covar "$covs" \
      --pheno "$pheno" \
      --model "NB" \
      --mode "nominal" \
      --window 500000 \
      --genelist $chunk \
      --test-method "score" \
      --nperm 1000 \
      --addpc 0 \
      --standardize \
      --out "jaxqtl-${cell_type}-${chr}-${chunk_id}"
    """
}

process COMPRESS_AND_INDEX_BED {
    conda "$projectDir/envs/cli.yaml"

    input: tuple val(chr), val(cell_type), val(bed)
    output: tuple val(chr), val(cell_type), 
        path("${cell_type}-${chr}-bed.gz"),
        path("${cell_type}-${chr}-bed.gz.tbi")

    script:
    """ 
    bgzip "$bed" -o "${cell_type}-${chr}-bed.gz"
    tabix -p bed "${cell_type}-${chr}-bed.gz"
    """
}

process TO_SPARSE_GRM_FORMAT {

    input: val grm
    output: path("onek1k-sparse-grm.tsv")

    script: 
    """
    to-sparse-grm-format.R $grm
    """
}

process RUN_APEX {
    label "apex"

    input: tuple val(chr), val(cell_type), val(pb_type), val(pheno_bed_gz), val(pheno_bed_gz_tbi), 
        val(vcf), val(vcf_tbi), val(covs), val(sparse_grm)
    output: tuple val(chr), val(cell_type), 
        path("apex-${cell_type}-${chr}.cis_gene_table.txt.gz"),
        path("apex-${cell_type}-${chr}.cis_sumstats.txt.gz")

    script: 
    """
    /home/jp2045/quasar-paper/apex/bin/apex cis \
        --vcf "$vcf" \
        --bed "$pheno_bed_gz" \
        --cov "$covs" \
        --grm "$sparse_grm" \
        --prefix "apex-${cell_type}-${chr}"
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
      -m  lm \
      --verbose
    gzip "${chr}-${cell_type}-lm-cis-variant.txt"
    gzip "${chr}-${cell_type}-lm-cis-region.txt"
    """
}

process PLOT_CONCORDANCE {
    publishDir "output"

    input:
        tuple val(cell_type), val(chrs), val(variants_list) 
        tuple val(cell_type), val(pairs_list)
    output: path("plot.pdf")

    script:
    """
    echo $variants_list | tr ',' '\n' | tr -d '[]' | sed 's/.*/"&"/' > quasar_files.txt
    echo $pairs_list | tr ',' '\n' | tr -d '[]' | sed 's/.*/"&"/' > tensorqtl_files.txt
    plot-concordance.R quasar_files.txt tensorqtl_files.txt
    """
}

process PLOT_POWER {
    publishDir "output"

    input:
        tuple val(cell_type), val(chrs), val(variants_list) 
        tuple val(cell_type), val(pairs_list)
    output: path("plot.pdf")

    script:
    """
    plot-power.R "${variants_list.collect()}" "${pairs_list.collect()}"
    """
}

