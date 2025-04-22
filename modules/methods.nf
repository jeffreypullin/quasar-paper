
process RUN_QUASAR {
    label "quasar"

    input: tuple val(chr), val(cell_type), val(pb_type), val(pheno_bed), 
        val(plink_bed), val(covs), val(grm), val(model)
    output: tuple val(chr), val(cell_type),
        path("${chr}-${cell_type}-${model}-cis-region.txt.gz"), 
        path("${chr}-${cell_type}-${model}-cis-variant.txt.gz"),
        path("${chr}-${cell_type}-${model}-time.txt")

    script:
    def prefix = "${plink_bed.getParent().toString() + '/' + plink_bed.getSimpleName()}"
    def grm_flag = (model == "lmm" || model == "glmm") ? "-g ${grm}" : " "
    def apl_flag = (model == "glm-apl") ? "--use-apl" : " "
    def passed_model = (model == "glm-apl") ? "glm" : model
    """
    /usr/bin/time -p -o "${chr}-${cell_type}-${model}-time.txt" \
      /home/jp2045/quasar/build/quasar \
      -p "$prefix" \
      -b "$pheno_bed" \
      -c "$covs" \
      ${grm_flag} \
      -o "${chr}-${cell_type}-${model}" \
      -m $passed_model \
      ${apl_flag} \
      --verbose
    gzip "${chr}-${cell_type}-${model}-cis-variant.txt"
    gzip "${chr}-${cell_type}-${model}-cis-region.txt"
    """
}

process RUN_APEX {
    label "apex"

    input: tuple val(chr), val(cell_type), val(pb_type), val(pheno_bed_gz), val(pheno_bed_gz_tbi), 
        val(vcf), val(vcf_tbi), val(covs), val(sparse_grm)
    output: tuple val(chr), val(cell_type), 
        path("apex-${cell_type}-${chr}.cis_gene_table.txt.gz"),
        path("apex-${cell_type}-${chr}.cis_long_table.txt.gz"),
        path("${cell_type}-${chr}-time.txt")        

    script: 
    """
    /usr/bin/time -p -o "${cell_type}-${chr}-time.txt" \
        /home/jp2045/quasar-paper/apex/bin/apex cis \
        --vcf "$vcf" \
        --bed "$pheno_bed_gz" \
        --cov "$covs" \
        --grm "$sparse_grm" \
        --prefix "apex-${cell_type}-${chr}" \
        --rankNormal \
        --long
    """
}

process RUN_JAXQTL_CIS {
    conda "$projectDir/envs/jaxqtl.yaml"
    label "jaxqtl"

    input: tuple val(chr), val(cell_type), val(pheno), val(chunk), val(chunk_id), val(bed), val(covs), val(chr_num)
    output: tuple val(chr), val(cell_type), 
        path("jaxqtl-${cell_type}-${chr}-${chunk_id}.cis_score.tsv.gz"),
        path("${cell_type}-${chr}-${chunk_id}-time.txt")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "${cell_type}-${chr}-${chunk_id}-time.txt" \
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
    output: tuple val(chr), val(cell_type), 
        path("jaxqtl-${cell_type}-${chr}-${chunk_id}.cis_qtl_pairs.${chr_num}.score.parquet"),
        path("${cell_type}-${chr}-${chunk_id}-time.txt")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "${cell_type}-${chr}-${chunk_id}-time.txt" \
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


process RUN_TENSORQTL_CIS {
    conda "$projectDir/envs/tensorqtl.yaml"
    label "tensorqtl_cis"

    input:  tuple val(cell_type), val(pheno), val(covs), val(bed)
    output: tuple val(cell_type),
        path("onek1k-${cell_type}.cis_qtl.txt.gz"),
        path("${cell_type}-time.txt")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "${cell_type}-time.txt" \
      python3 -m tensorqtl "$prefix" "$pheno" "onek1k-${cell_type}" \
      --covariates "$covs" \
      --mode cis
    """
}

process RUN_TENSORQTL_CIS_NOMINAL {
    conda "$projectDir/envs/tensorqtl.yaml"
    label "tensorqtl"

    input: tuple val(cell_type), val(pheno), val(covs), val(bed)
    output: tuple val(cell_type), 
        path("onek1k-${cell_type}.cis_qtl_pairs.*.parquet"),
        path("${cell_type}-time.txt")

    script:
    def prefix = "${bed.getParent().toString() + '/' + bed.getSimpleName()}"
    """
    /usr/bin/time -p -o "${cell_type}-time.txt" \
       python3 -m tensorqtl "$prefix" "$pheno" "onek1k-${cell_type}" \
      --covariates "$covs" \
      --mode cis_nominal
    """
}