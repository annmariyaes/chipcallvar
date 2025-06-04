// chipcallvar/modules/nf-core/macs3/callvar/main.nf

process MACS3_CALLVAR {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/macs3/variant_calling/${meta.patient}", mode: 'copy'
    container "${params.MACS3_CONTAINER}" 

    input:
    tuple val(meta), path(peaks), path(treat_bams), path(treat_bais), path(ctrl_bams), path(ctrl_bais)
    
    output:
    tuple val(meta), path("${meta.patient}.macs3.vcf"), emit: vcf

    script:
    def ctrl_flag = ctrl_bams ? "--control $ctrl_bams" : ''   

    """
    macs3 callvar \
        --peak ${peaks} \
        --treatment ${treat_bams} \
        ${ctrl_flag} \
        --multiple-processing ${task.cpus} \
        --outdir . \
        --ofile ${meta.patient}.macs3.vcf 
    """
}
