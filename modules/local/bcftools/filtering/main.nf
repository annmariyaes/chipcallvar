/*
=============================================================================
        MODULE: BCFTOOLS filtering
=============================================================================
*/



process BCFTOOLS {
    tag "$meta.patient"
    publishDir "${params.OUTDIR}/annotation/${caller}/${meta.patient}", mode: 'copy'
    
    input:
    tuple val(meta), path(vcf)
    val caller
    
    output:
    tuple val(meta), path("${meta.patient}.${caller}.filtered.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.patient}.${caller}.filtered.vcf.gz.tbi"), emit: vcf_tbi
    
    script:
    """
    # Step 1: Fill FORMAT/VAF tag
    bcftools +fill-tags ${vcf} -Oz -o ${meta.patient}_filled.vcf.gz -t FORMAT/VAF
    
    # Step 2: Filter using split-vep
    bcftools +split-vep -i 'INFO/DP>=${params.DP} & FORMAT/VAF>=${params.VAF} & (gnomADe_AF<=${params.gnomADe_AF} | gnomADe_AF==".")' \
                        ${meta.patient}_filled.vcf.gz -Oz -o ${meta.patient}.${caller}.filtered.vcf.gz
    
    # Step 3: Index the filtered VCF
    bcftools index -t ${meta.patient}.${caller}.filtered.vcf.gz
       
    # Step 4: Remove the intermediate filled vcf file
    rm "${meta.patient}_filled.vcf.gz"
    """
}
