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

    when:
    'freebayes' in ${params.TOOLS}

    input:
    tuple val(meta), path(treat_bams), path(treat_bais)

    output:
    tuple val(meta), path("${meta.patient}.mutect2.vcf.gz"), emit: vcf

    script:
    """
    freebayes \
        -f {params.GENOME} \
        -b ${treat_bams} \
        -v \
        -O ${meta.patient}.freebayes.vcf.gz
    """
}
