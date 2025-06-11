/*
=============================================================================
        MODULE: GATK mutect2
=============================================================================
*/


process GATK_MUTECT2 {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/variant_calling/mutect2/${meta.patient}", mode: 'copy'
    container "${params.GATK_CONTAINER}"
    label 'process_high'
    
    input:
    tuple val(meta), path(treat_bams), path(treat_bais), path(reference), path("*")

    output:
    tuple val(meta), path("${meta.patient}.mutect2.vcf.gz"), emit: vcf

    script:
    """
    gatk Mutect2 \
	-R ${reference} \
  	-I ${treat_bams} \
  	--germline-resource ${params.GNOMAD} \
  	--panel-of-normals ${params.PON} \
  	-O ${meta.patient}.mutect2.vcf.gz
    """
}
