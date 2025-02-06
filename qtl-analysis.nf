#!/usr/bin/env nextflow

params.raw_cov_data = file("randolph-data/individual_meta_data_for_GE_with_scaledCovars_with_geneProps.txt")
params.raw_feat_anno_data = file("randolph-data/GRCh38.92_gene_positions.txt")
params.raw_genotype_data = file("randolph-data/genotypes.txt")
workflow {

    PROCESS_COVARIATE_DATA(params.raw_cov_data)
    PROCESS_FEAT_ANNO_DATA(params.raw_feat_anno_data)
    PROCESS_GENOTYPE_DATA(params.raw_genotype_data)
}

process PROCESS_COVARIATE_DATA {

    input: val raw_cov_data
    output: path("cov-data.tsv")

    script: 
    """
    randolph-process-covariates.R ${raw_cov_data}
    """
}

process PROCESS_FEAT_ANNO_DATA {

    input: val raw_feat_anno_data
    output: path("feat-anno-data.tsv")

    script:
    """
    randolph-process-feat-anno.R ${raw_feat_anno_data}
    """
}

process PROCESS_GENOTYPE_DATA {

    input: val randolph_genotype_data
    output: tuple path("randolph.ped"), path("randolph.map")

    // First convert the file to a ped file than into a bed.
    script:
    """
    module load plink/1.9

    randolph-process-genotype.R ${randolph_genotype_data}
    plink --file randolph --make-bed --out randolph
    """
}

