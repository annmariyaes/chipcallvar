/*
=============================================================================
        MODULE: GATK mutect2
=============================================================================
*/


process GATK_MUTECT2 {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/variant_calling/mutect2/${meta.patient}", mode: 'copy'
    container "${params.GATK_CONTAINER}"
    
    input:
    tuple val(meta), path(treat_bams), path(treat_bais), path(reference), path("*"), path(dict), path(intervals)

    output:
    tuple val(meta), path("${meta.patient}.${interval}.mutect2.vcf.gz"), emit: vcf

    script:
    interval = intervals.baseName
    """
    gatk Mutect2 \\
        -R ${reference} \\
        -I ${treat_bams} \\
        -L ${intervals} \\
        --germline-resource ${params.GNOMAD} \\
        --panel-of-normals ${params.PON} \\
        --sequence-dictionary ${dict} \\
        -O ${meta.patient}.${interval}.mutect2.vcf.gz
    """
}
