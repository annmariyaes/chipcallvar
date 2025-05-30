/*
=============================================================================
        MODULE: VEP
=============================================================================
*/

process ENSEMBL_VEP {
    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/annotation/ensembl-vep/${meta.sample}", mode: 'symlink'
    container "${params.VEP_CONTAINER}" 

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.sample}_vep_ann.vcf"), emit: vcf

    script:
    """
    # Run VEP with optimized parameters
    vep \
        -i ${vcf} \
        -o ${meta.sample}_vep_ann.vcf \
        --cache \
        --dir_cache ${params.VEP_CACHE} \
        --fasta ${params.REFERENCE_FASTA} \
        --everything \
        --cell_type list \
        --check_existing \
        --vcf \
        --fork ${task.cpus} \
        --buffer_size 10000 
    """
}
