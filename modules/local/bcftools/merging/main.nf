/*
=============================================================================
        MODULE: BCFTOOLS merging

https://samtools.github.io/bcftools/bcftools.html#merge
=============================================================================
*/



process BCFTOOLS_MERGING {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/annotation/ensembl_vep/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.id}.filtered.vcf.gz.tbi"), emit: vcf_tbi

    script:
    """
    # Step 1: Fill FORMAT/VAF tag
    bcftools +fill-tags ${vcf} -Oz -o ${meta.id}_filled.vcf.gz -- -t FORMAT/VAF

    # Step 2: Filter using split-vep
    bcftools +split-vep -i 'INFO/DP>=${params.DP} & FORMAT/VAF>=${params.VAF} & (gnomADe_AF<=${params.gnomADe_AF} | gnomADe_AF==".")' \
                        ${meta.id}_filled.vcf.gz -Oz -o ${meta.id}.{meta.caller}.filtered.vcf.gz

    # Step 3: Index the filtered VCF
    bcftools index ${meta.id}.{meta.caller}.filtered.vcf.gz

    # Step 4: Remove the intermediate filled vcf file
    rm "${meta.id}_filled.vcf.gz"
    """
}
