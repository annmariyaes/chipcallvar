/*
=============================================================================
        MODULE: Freebayes
=============================================================================
*/

process FREEBAYES {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/variant_calling/freebayes/${meta.patient}", mode: 'copy'
    container "${params.FREEBAYES_CONTAINER}"
    label 'process_high'
    
    input:
    tuple val(meta), path(treat_bams), path(treat_bais), path(reference), path("*"), path(dict)
    
    output:
    tuple val(meta), path("${meta.patient}.freebayes.vcf.gz"), emit: vcf
    
    script:
    """
    freebayes \
        -f ${reference} \
        -b ${treat_bams} \
        -v ${meta.patient}.freebayes.vcf
    
    # Compress and index the VCF file
    bgzip ${meta.patient}.freebayes.vcf
    tabix -p vcf ${meta.patient}.freebayes.vcf.gz
    """
}
