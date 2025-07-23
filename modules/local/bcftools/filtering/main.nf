/*
=============================================================================
        MODULE: BCFTOOLS filtering
=============================================================================
*/



process BCFTOOLS {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/annotation/${caller}/${meta.id}", mode: 'copy'
    container "${params.BCFTOOLS_CONTAINER}"
    
    input:
    tuple val(meta), path(vcf)
    val caller
    
    output:
    tuple val(meta), path("${meta.id}.${caller}.filtered.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.${caller}.filtered.vcf.gz.tbi"), emit: vcf_tbi
    
    script:
    def prefix = "${meta.id}.${caller}"
    def filled_vcf = "${prefix}_filled.vcf.gz"
    def filtered_vcf = "${prefix}.filtered.vcf.gz"
    
    """
    # Fill FORMAT/VAF tag
    bcftools +fill-tags ${vcf} -Oz -o ${filled_vcf} -- -t FORMAT/VAF
    
    # Filter using split-vep with improved filter expression
    bcftools +split-vep -i 'FORMAT/DP>=${params.DP} & FORMAT/VAF>=${params.VAF} & (AF<=${params.AF} | AF==".") & (gnomADe_AF<=${params.gnomADe_AF} | gnomADe_AF==".")' ${filled_vcf} -Oz -o ${filtered_vcf}
    
    # Index the filtered VCF
    bcftools index -t ${filtered_vcf}
    
    # Clean up intermediate files
    rm -f ${filled_vcf}
    
    # Validate output
    # bcftools view -h ${filtered_vcf} | tail -1 | grep -q "^#CHROM" || exit 1
    """
}
