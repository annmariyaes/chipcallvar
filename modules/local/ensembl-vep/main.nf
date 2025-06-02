/*
=============================================================================
        MODULE: VEP
=============================================================================
*/

process ENSEMBL_VEP {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/annotation/ensembl-vep/${meta.patient}", mode: 'symlink'
    container "${params.VEP_CONTAINER}" 

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.patient}_vep_ann.vcf"), emit: vcf

    script:
    """
    # Run VEP with optimized parameters
    vep \
        -i ${vcf} \
        -o ${meta.patient}_vep_ann.vcf \
        --cache \
        --dir_cache ${params.VEP_CACHE} \
        --fasta ${params.REFERENCE_FASTA} \
        --everything \
        --cell_type list \
        --check_existing \
        --vcf \
        --fork ${task.cpus} \
        --buffer_size 10000 \
        --stats_file ${meta.patient}.variant_effect_output.txt_summary.html
    """
}
