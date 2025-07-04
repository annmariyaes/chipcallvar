/*
=============================================================================
        MODULE: Freebayes variant calling
=============================================================================
*/

process FREEBAYES {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/variant_calling/freebayes/${meta.patient}", mode: 'copy'
    container "${params.FREEBAYES_CONTAINER}"
    label 'process_high'
    
    input:
    tuple val(meta), path(peaks), path(treat_bams), path(treat_bais), path(reference), path("*"), path(dict), path(intervals)
    
    output:
    tuple val(meta), path("${meta.patient}.${interval}.freebayes.vcf"), emit: vcf
    
    script:
    interval = intervals.baseName
    """
    # create target bed file from peaks
    cut -f1-3 ${peaks} > ${meta.patient}.bed
    # freebayes variant calling
    freebayes \
        -f ${reference} \
        -b ${treat_bams} \
        -t ${meta.patient}.bed \
        > ${meta.patient}.${interval}.freebayes.vcf    
    """
}
