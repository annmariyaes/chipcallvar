/*
=============================================================================
        MODULE: VEP
=============================================================================
*/

process ENSEMBL_VEP {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/annotation/ensembl_vep/${meta.patient}", mode: 'copy'
    container "${params.VEP_CONTAINER}" 

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.patient}_vep_ann.vcf"), emit: vcf
    path("${meta.patient}_vep_summary.html"), emit: vcf_stats

    
    script:
    """
    # Run VEP with optimized parameters
    vep \
        -i ${vcf} \
        -o ${meta.patient}_vep_ann.vcf \
        --cache \
        --dir_cache ${params.VEP_CACHE} \
        --everything \
        --check_existing \
        --vcf \
        --fork ${task.cpus} \
        --buffer_size 10000 \
        --stats_file ${meta.patient}_vep_summary.html  
    """
}
