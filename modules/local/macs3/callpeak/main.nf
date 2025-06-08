/*
=============================================================================
        MODULE: MACS3 callpeak
=============================================================================
*/


process MACS3_CALLPEAK {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/macs3/peak_calls/${meta.patient}", mode: 'copy'
    container "${params.MACS3_CONTAINER}"
    
    input:
    tuple val(meta), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais)
    
    output:
    tuple val(meta), path("${meta.patient}_peaks.narrowPeak"), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais), emit: peaks
    tuple val(meta), path("${meta.patient}_peaks.xls"), path("${meta.patient}_summits.bed"), emit: peak_stats

    script:
    def format = meta.single_end ? "BAM" : "BAMPE"
    def ctrl_flag = ctrl_bams ? "--control $ctrl_bams" : ''
    """
    macs3 callpeak \
        --treatment ${treat_bams} \
        ${ctrl_flag} \
        --format ${format} \
        --gsize hs \
        --name ${meta.patient} \
        --outdir .
    """
}
