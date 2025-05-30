process MACS3_CALLPEAK {
    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/peak_calls/macs3/${meta.sample}", mode: 'copy'
    container "${params.MACS3_CONTAINER}"
    
    input:
    tuple val(meta), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais)
    
    output:
    tuple val(meta), path("${meta.sample}_peaks.narrowPeak"), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais), emit: peaks
    
    script:
    def format = meta.single_end ? "BAM" : "BAMPE"
    def ctrl_flag = ctrl_bams ? "--control $ctrl_bams" : ''
    """
    macs3 callpeak \\
        -t ${treat_bams} \\
        ${ctrl_flag} \\
        -f ${format} \\
        -g hs \\
        -n ${meta.sample} \\
        --outdir .
    """
}
