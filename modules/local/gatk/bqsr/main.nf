




process GATK_BQSR {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/preprocessing/recalibrated/${meta.sample}", mode: 'copy'
    container "${params.GATK_CONTAINER}"

    input:
    tuple val(meta), path(treat_bams), path(treat_bais), path(dbsnp), path(dbindel)

    output:
    tuple val(meta), path("${meta.sample}.recal.bam"), path("${meta.sample}.recal.bam.bai"), path("${meta.sample}.recal.table"), emit: recalibrated

    script:
    """
    gatk BaseRecalibrator \
  	-I "${treat_bams}" \
  	-R "${reference}" \
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
