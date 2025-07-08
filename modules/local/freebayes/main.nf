/*
=============================================================================
        MODULE: Freebayes variant calling
=============================================================================
*/

process FREEBAYES {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/variant_calling/freebayes/${meta.id}", mode: 'copy'
    container "${params.FREEBAYES_CONTAINER}"
    label 'process_high'
    
    input:
    tuple val(meta), path(peaks), path(treat_bams), path(treat_bais), path(reference), path("*"), path(dict)

    output:
    tuple val(meta), path("${meta.id}.freebayes.vcf"), emit: vcf
    
    script:
    
    """
    # create target bed file from peaks
    cut -f1-3 ${peaks} > ${meta.id}.bed
    # freebayes variant calling
    freebayes \
        -f ${reference} \
        -b ${treat_bams} \
        -t ${meta.id}.bed \
        > ${meta.id}.freebayes.vcf    
    """
}
