/*
=============================================================================
        MODULE: GATK BaseRecalibrator, ApplyBQSR 
Sequencers make systematic errors in assigning base quality scores. 
To correct for these errors, a model is built using covariates encoded in the read groups from all base calls and then applying the adjustments to generate recalibrated base qualities
=============================================================================
*/


process GATK_BQSR {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/preprocessing/recalibrated/${meta.sample}", mode: 'copy'
    container "${params.GATK_CONTAINER}"

    input:
    tuple val(meta), path(treat_bams), path(treat_bais), path(reference), path(fai), path(dict), path(dbsnp), path(dbindel)

    output:
    tuple val(meta), path("${meta.sample}.recal.bam"), path("${treat_bais}"), emit: recalibrated
    path("${meta.sample}.recal.table"), emit: table

 
    script:
    """
    gatk BaseRecalibrator \
        -R "${reference}" \
   	-I "${treat_bams}" \
  	--known-sites "${dbsnp}" \
  	--known-sites "${dbindel}" \
  	-O "${meta.sample}.recal.table"
    
    gatk ApplyBQSR \
  	-R "${reference}" \
  	-I "${treat_bams}" \
  	--bqsr-recal-file "${meta.sample}.recal.table" \
  	-O "${meta.sample}.recal.bam"
    """
}
