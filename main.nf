#!/usr/bin/env nextflow

params.onek1k_raw_single_cell_data = file("data/onek1k/OneK1K_cohort_gene_expression_matrix_14_celltypes.h5ad.gz")
params.onek1k_supp_tables = file("data/onek1k/science.abf3041_tables_s6_to_s19/science.abf3041_tables_s6_to_s19.xlsx")
params.onek1k_feature_annot = file("data/onek1k/Homo_sapiens.GRCh37.82.bed")
params.onek1k_cell_types = [
    'CD4 NC', 'CD4 ET', 'CD4 SOX4', 'CD8 ET', 'CD8 NC', 'CD8 S100B', 'NK', 
    'NK R', 'Plasma', 'B Mem', 'B IN', 'Mono C', 'Mono NC', 'DC'
]
params.pb_types = ['mean', 'sum']
params.quasar_spec = "data/quasar-spec.tsv"

include { RUN_APEX; RUN_QUASAR; RUN_JAXQTL_CIS; RUN_JAXQTL_CIS_NOMINAL; RUN_TENSORQTL_CIS; RUN_TENSORQTL_CIS_NOMINAL } from './modules/methods'
include { RUN_QUASAR as RUN_QUASAR_PERM } from './modules/methods'
include { RUN_JAXQTL_CIS_NOMINAL as RUN_JAXQTL_CIS_NOMINAL_PERM} from './modules/methods'
include { RUN_TENSORQTL_CIS_NOMINAL as RUN_TENSORQTL_CIS_NOMINAL_PERM} from './modules/methods'
include { RUN_APEX as RUN_APEX_PERM} from './modules/methods'

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
    pheno_bed = ONEK1K_MAKE_PHENO_BED(cell_type_pb, params.onek1k_feature_annot)

    // Run quasar.
    quasar_spec = channel
      .fromPath(params.quasar_spec)
      .splitCsv()
      .combine(chrs)
      .combine(cell_types)
      .map({ it -> [it[2], it[3], it[0], it[1]] })

    quasar_input = pheno_bed
      .combine(chrs)
      .map({ it -> [it[3], it[0], it[1], it[2]] })
      .combine(bed_files, by: 0)
      .combine(chrs.combine(covs), by: [0, 1])
      .combine(grm)
      .combine(quasar_spec, by: [0, 1, 2])
      .filter( { it -> it[1] == "B IN" || it[1] == "Plasma" || it[1] == "CD4 NC"})
      .filter( { it -> it[7] != "nb_glmm"})
      .filter( { it -> it[7] != "nb_glmm-apl"})

    quasar = RUN_QUASAR(quasar_input)

    // Run TensorQTL.
    norm_pheno_bed = NORMALISE_PHENO_BED(pheno_bed
        .filter( {it -> it[1] == "mean"} )
        .map( { it -> [it[0], it[2]] } )
    )
    t_covs = ONEK1K_TRANSPOSE_COVS(covs)

    tensorqtl_input = norm_pheno_bed
      .join(t_covs)
      .combine(all_bed)
      .filter( { it -> it[0] == "B IN" || it[0] == "Plasma" || it[0] == "CD4 NC"})

    tensorqtl_cis = RUN_TENSORQTL_CIS(tensorqtl_input)
    tensorqtl_cis_nominal = RUN_TENSORQTL_CIS_NOMINAL(tensorqtl_input)
    
    // Run jaxQTL.
    jaxqtl_spec = pheno_bed
      .filter({ x -> x[1] == "sum" })
      .combine(chrs)
      .map({ it -> [it[3], it[0], it[2]] })
    
    chunks = ONEK1K_MAKE_GENELIST_CHUNKS(jaxqtl_spec, 100)
    jaxqtl_covs = PROCESS_JAXQTL_COVARIATES(covs)
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
        .combine(chrs.combine(jaxqtl_covs), by: [0, 1])
        .filter( { it -> it[1] == "B IN" || it[1] == "Plasma" || it[1] == "CD4 NC"})
        .map( { it -> [it, it[0].substring(3).toInteger()].flatten()})

    jaxqtl_cis = RUN_JAXQTL_CIS(jaxqtl_input)
    jaxqtl_cis_nominal = RUN_JAXQTL_CIS_NOMINAL(jaxqtl_input)

    // Run apex.
    pheno_bed_gz = COMPRESS_AND_INDEX_BED(pheno_bed)
    apex_covs = PROCESS_APEX_COVS(t_covs)
    sparse_grm = TO_SPARSE_GRM_FORMAT(grm)

    apex_input = pheno_bed_gz
        .filter({ x -> x[1] == "mean" })
        .combine(chrs)
        .map({ it -> [it[4], it[0], it[1], it[2], it[3]] })
        .combine(filt_vcf_files, by: 0)
        .combine(chrs.combine(apex_covs), by: [0, 1])
        .combine(sparse_grm)
        .filter( { it -> it[1] == "B IN" || it[1] == "Plasma" || it[1] == "CD4 NC"})
    
    apex = RUN_APEX(apex_input) 

    // Run permutation analysis
    rep_bed_files = bed_files
      .filter({ it -> it[0] == "chr22"})
      .combine(channel.of(1..10))

    permute_bed_files = PERMUTE_BED(rep_bed_files)

    rep_vcf_files = filt_vcf_files
      .filter({ it -> it[0] == "chr22"})
      .combine(channel.of(1..10))

    permute_vcf_files = PERMUTE_VCF(rep_vcf_files)
    
    permute_quasar_input = quasar_input
      .filter( {it -> it[0] == "chr22" })
      .combine(permute_bed_files, by: 0)
      .map({ it -> [it[0], it[1], it[2], it[3], it[8], it[5], it[6], it[7]]})

    quasar_perm = RUN_QUASAR_PERM(permute_quasar_input)

    permute_jaxqtl_input = jaxqtl_input
      .filter({it[0] == "chr22"})
      .combine(permute_bed_files, by: 0)
      .map({ it -> [it[0], it[1], it[2], it[3], it[4], it[8], it[6], it[7]]})

    jaxqtl_cis_nominal_perm = RUN_JAXQTL_CIS_NOMINAL_PERM(permute_jaxqtl_input)

    permute_tensorqtl_input = tensorqtl_input
      .combine(permute_bed_files)
      .map({ it -> [it[0], it[1], it[2], it[5]]})

    tensorqtl_cis_nominal_perm = RUN_TENSORQTL_CIS_NOMINAL_PERM(permute_tensorqtl_input)

    permute_apex_input = apex_input
      .filter( { it -> it[0] == "chr22"})
      .combine(permute_vcf_files)
      .map( {it -> [it[0], it[1], it[2], it[3], it[4], it[10], it[11], it[7], it[8]]})
    
    apex_perm = RUN_APEX_PERM(permute_apex_input)

    // Reformat output.
    ind_channel = channel.of(1)

    quasar_grouped = ind_channel
        .combine(quasar)
        .groupTuple()

    jaxqtl_cis_nominal_grouped = ind_channel
        .combine(jaxqtl_cis_nominal)
        .groupTuple() 
    
    jaxqtl_cis_grouped = ind_channel
        .combine(jaxqtl_cis)
        .groupTuple() 

    apex_grouped = ind_channel
        .combine(apex)
        .groupTuple()

    tensorqtl_cis_nominal_grouped = ind_channel
        .combine(tensorqtl_cis_nominal)
        .groupTuple()
        .map( { it -> [it[0], it[1], it[2].flatten(), it[3]]})
    
    tensorqtl_cis_grouped = ind_channel
        .combine(tensorqtl_cis)
        .groupTuple()

    // Permutation data.
    quasar_perm_grouped = ind_channel
        .combine(quasar_perm)
        .groupTuple()

    jaxqtl_cis_nominal_perm_grouped = ind_channel
        .combine(jaxqtl_cis_nominal_perm)
        .groupTuple()

    tensorqtl_cis_nominal_perm_grouped  = ind_channel
        .combine(tensorqtl_cis_nominal_perm)
        .groupTuple()
        .map( { it -> [it[0], it[1], it[2].flatten(), it[3]]})

    apex_perm_grouped = ind_channel
        .combine(apex_perm)
        .groupTuple()

    // Make plots.
    //PLOT_CONCORDANCE(
    //    quasar_grouped, 
    //    tensorqtl_cis_nominal_grouped,
    //    jaxqtl_cis_nominal_grouped,
    //    apex_grouped
    //)
    
    PLOT_POWER(
        quasar_grouped,
        tensorqtl_cis_nominal_grouped,
        jaxqtl_cis_nominal_grouped,
        apex_grouped
    )     
    
    PLOT_TIME(
        quasar_grouped,
        tensorqtl_cis_nominal_grouped,
        tensorqtl_cis_grouped,
        jaxqtl_cis_nominal_grouped,
        jaxqtl_cis_grouped,
        apex_grouped
    )

    PLOT_FDR(
       quasar_perm_grouped,
       tensorqtl_cis_nominal_perm_grouped,
       jaxqtl_cis_nominal_perm_grouped,
       apex_perm_grouped
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

process PROCESS_JAXQTL_COVARIATES {

    input: tuple val(cell_type), val(covs)
    output: tuple val(cell_type), path("onek1k-${cell_type}-all-jaxqtl-covs.tsv")

    script: 
    """
    process-jaxqtl-covariates.R "$covs" "$cell_type"
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

process NORMALISE_PHENO_BED {

   input: tuple val(cell_type), val(pheno_bed)
   output: tuple val(cell_type), path("onek1k-${cell_type}-pheno-norm.bed")

   script:
   """
   onek1k-norm-pheno-bed.R "$pheno_bed" "$cell_type"
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

process PERMUTE_BED {

    input: tuple val(chr), val(plink_bed), val(ind)
    output: tuple val(chr), path("permute-${chr}.bed")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    """
    permute-fam.R ${prefix}.fam ${chr}
    cp ${prefix}.bim ./"permute-${chr}.bim"
    cp ${prefix}.bed ./"permute-${chr}.bed"
    """
}

process PERMUTE_VCF {
    conda "$projectDir/envs/cli.yaml"

    input: tuple val(chr), val(vcf), val(tbi), val(ind)
    output: tuple val(chr), path("permute-${chr}.vcf.gz"), path("permute-${chr}.vcf.gz.tbi")

    script: 
    """
    permute-vcf.R $vcf $chr
    gunzip "permute-${chr}.vcf.gz"
    bgzip "permute-${chr}.vcf"
    tabix -p vcf "permute-${chr}.vcf.gz"
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

process PROCESS_APEX_COVS {

    input: tuple val(cell_type), val(covs)
    output: tuple val(cell_type), path("onek1k-${cell_type}-all-apex-covs.tsv")

    script: 
    """
    process-apex-covariates.R "$covs" "$cell_type"
    """
}

process PLOT_CONCORDANCE {
    publishDir "output"

    input:
        tuple val(cell_type), val(chrs), val(quasar_pairs_list), val(quasar_time)
        tuple val(cell_type), val(tensorqtl_pairs_list), val(tensorqtl_time)
        tuple val(cell_type), val(chrs), val(jaxqtl_pairs_list), val(jaxqtl_time)
        tuple val(cell_type), val(chrs), val(apex_pairs_list), val(apex_time)
    output: tuple path("plot-lm-concordance.pdf"), 
        path("plot-glm-concordance.pdf"), 
        path("plot-lmm-concordance.pdf")

    script:
    """
    plot-concordance.R \
        "${quasar_pairs_list.collect()}" \
        "${tensorqtl_pairs_list.collect()}" \
        "${jaxqtl_pairs_list.collect()}" \
        "${apex_pairs_list.collect()}"
    """
}

process PLOT_POWER {
    publishDir "output"

    input: 
        tuple val(ind), val(cell_type), val(chrs), val(quasar_region), val(quasar_pairs_list), val(quasar_time)
        tuple val(ind), val(cell_type), val(tensorqtl_pairs_list), val(tensorqtl_time)
        tuple val(ind), val(cell_type), val(chrs), val(jaxqtl_pairs_list), val(jaxqtl_cis_nominal_time)
        tuple val(ind), val(cell_type), val(chrs), val(apex_region_list), val(apex_pairs_list), val(apex_time)
    output: path("plot-power.pdf")

    script:
    """
    plot-power.R \
        "${quasar_pairs_list.collect()}" \
        "${tensorqtl_pairs_list.collect()}" \
        "${jaxqtl_pairs_list.collect()}" \
        "${apex_pairs_list.collect()}"
    """
}

process PLOT_TIME {
    publishDir "output"

    input:
        tuple val(ind), val(cell_type), val(chrs), val(quasar_region), val(quasar_pairs_list), val(quasar_time)
        tuple val(ind), val(cell_type), val(tensorqtl_pairs_list), val(tensorqtl_cis_nominal_time)
        tuple val(ind), val(cell_type), val(tensorqtl_cis_list), val(tensorqtl_cis_time)
        tuple val(ind), val(cell_type), val(chrs), val(jaxqtl_pairs_list), val(jaxqtl_cis_nominal_time)
        tuple val(ind), val(cell_type), val(chrs), val(jaxqtl_cis_list), val(jaxqtl_cis_time)
        tuple val(ind), val(cell_type), val(chrs), val(apex_region_list), val(apex_pairs_list), val(apex_time)
    output: path("plot-time.pdf")

    script:
    """
    plot-time.R \
        "${quasar_time.collect()}" \
        "${tensorqtl_cis_nominal_time.collect()}" \
        "${tensorqtl_cis_time.collect()}" \
        "${jaxqtl_cis_nominal_time.collect()}" \
        "${jaxqtl_cis_time.collect()}" \
        "${apex_time.collect()}"
    """
}

process PLOT_FDR {
    label "plot"
    publishDir "output"

    input:
       tuple val(ind), val(cell_type), val(chrs), val(quasar_region), val(quasar_pairs_list), val(quasar_time)
       tuple val(ind), val(cell_type), val(tensorqtl_pairs_list), val(tensorqtl_cis_nominal_time)
       tuple val(ind), val(cell_type), val(chrs), val(jaxqtl_pairs_list), val(jaxqtl_cis_nominal_time)
       tuple val(ind), val(cell_type), val(chrs), val(apex_region_list), val(apex_pairs_list), val(apex_time)
    output: path("plot-other-fdr.pdf")

    script:
    """
    plot-fdr.R \
        "${quasar_pairs_list.collect()}" \
        "${tensorqtl_pairs_list.collect()}" \
        "${jaxqtl_pairs_list.collect()}" \
        "${apex_pairs_list.collect()}"
    """
}
