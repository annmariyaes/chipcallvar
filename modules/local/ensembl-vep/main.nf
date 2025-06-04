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
    tuple val(meta), path("${meta.patient}.macs3.vep.vcf"), emit: vcf
    tuple val(meta), path("${meta.patient}.vep.summary.html"), emit: vep_stats
    
    script:
    """
    # Run VEP with optimized parameters
    vep \
        -i ${vcf} \
        -o ${meta.patient}.macs3.vep.vcf \
        --cache \
        --dir_cache ${params.VEP_CACHE} \
        --everything \
        --check_existing \
        --vcf \
        --fork ${task.cpus} \
        --buffer_size 10000 \
        --stats_file ${meta.patient}.vep.summary.html  
    """
}
