/*
=============================================================================
        MODULE: GATK mutect2
=============================================================================
*/


process GATK_MUTECT2 {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/variant_calling/mutect2/${meta.id}", mode: 'copy'
    container "${params.GATK_CONTAINER}"
    
    input:
    tuple val(meta), path(peaks), path(treat_bams), path(treat_bais), path(reference), path("*"), path(dict)

    output:
    tuple val(meta), path("${meta.id}.mutect2.vcf.gz"), emit: vcf

    script:
    """
    # create target bed file from peaks
    cut -f1-3 ${peaks} > ${meta.id}.bed
    # Mutect2 variant calling 
    gatk Mutect2 \
        -R ${reference} \
        -I ${treat_bams} \
        --germline-resource ${params.GNOMAD} \
        --panel-of-normals ${params.PON} \
        --sequence-dictionary ${dict} \
        --intervals ${meta.id}.bed \
        -O ${meta.id}.mutect2.vcf.gz
    """
}
