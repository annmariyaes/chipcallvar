/*
=============================================================================
        MODULE: VEP
=============================================================================
*/


process ENSEMBL_VEP {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/annotation/${caller}/${meta.patient}", mode: 'copy'
    container "${params.VEP_CONTAINER}" 

    input:
    tuple val(meta), path(vcf)
    val caller

    output:
    tuple val(meta), path("${meta.patient}.${caller}.vep.vcf"), emit: vcf
    tuple val(meta), path("${meta.patient}.${caller}.vep.summary.html"), emit: vep_stats
    
    script:
    """
    # Run VEP with optimized parameters
    vep \
        -i ${vcf} \
        -o ${meta.patient}.${caller}.vep.vcf \
        --cache \
        --everything \
        --check_existing \
        --force_overwrite \
        --vcf \
        --fork ${task.cpus} \
        --buffer_size 10000 \
        --stats_html \
        --stats_file ${meta.patient}.${caller}.vep.summary.html \
        --verbose  
    """
}

