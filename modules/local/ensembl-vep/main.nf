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
        --everything \
        --check_existing \
        --force_overwrite \
        --vcf \
        --fork ${task.cpus} \
        --buffer_size 10000 \
        --stats_html \
        --stats_file ${meta.patient}.vep.summary.html \
        --verbose  
    """
}

