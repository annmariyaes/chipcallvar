// chipcallvar/modules/nf-core/macs3/callvar/main.nf

process MACS3_CALLVAR {
    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/variant_calling/macs3/${meta.sample}", mode: 'copy'
    container "${params.MACS3_CONTAINER}" 

    input:
    tuple val(meta), path(treat_bams), path(ctrl_bams) 

    output:
    tuple val(meta), path("${meta.sample}_peaks.vcf"), emit: vcf

    script:
    def ctrl_flag = ctrl_bams ? "-c ${ctrl_bams}" : "" 
   
    """
    macs3 callvar \
        -b ${peaks} \
        -t ${treat_bams} \
        ${ctrl_flag} \
        -m 8 \
        -o ${meta.sample}_peaks.vcf
    """
}
