/*
=============================================================================
        MODULE: MACS3 callvar
=============================================================================
*/


process MACS3_CALLVAR {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/variant_calling/macs3/${meta.patient}", mode: 'copy'
    container "${params.MACS3_CONTAINER}"
    
    input:
    tuple val(meta), path(peaks), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais), path(intervals)
    
    output:
    tuple val(meta), path("${meta.patient}_${interval_name}.macs3.vcf"), emit: vcf
    
    script:
    interval_name = intervals.baseName
    def ctrl_flag = ctrl_bams ? "--control $ctrl_bams" : ''
    """
    # create target bed file from peaks
    cut -f1-3 ${peaks} > ${meta.patient}.peaks.bed
    sort -k1,1 -k2,2n ${meta.patient}.peaks.bed > ${meta.patient}.sorted.peaks.bed
    # variant calling
    macs3 callvar \
        --peak ${meta.patient}.sorted.peaks.bed \
        --treatment ${treat_bams} \
        ${ctrl_flag} \
        --multiple-processing ${task.cpus} \
        --outdir . \
        --ofile ${meta.patient}_${interval_name}.macs3.vcf \
    """
}
