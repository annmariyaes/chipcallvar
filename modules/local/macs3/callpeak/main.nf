process MACS3_CALLPEAK {
    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/peak_calls/macs3/${meta.sample}", mode: 'copy'
    container "${params.MACS3_CONTAINER}"

    input:
    tuple val(meta), path(treat_bams), path(ctrl_bams)

    output:
    tuple val(meta), path("${meta.sample}_peaks.narrowPeak"), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais), emit: peaks

    script:
    def format = meta.single_end ? "BAM" : "BAMPE"
    def ctrl_flag = ctrl_bam ? "--control $ctrl_bam" : ''

    """
    macs3 callpeak \
      	-t ${treat_bams} \\
      	${ctrl_flag} \\
      	-f BAM \\
      	-g hs \\
      	-n ${meta.sample} \\
      	--outdir .
    """
}
