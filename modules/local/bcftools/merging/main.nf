/*
=============================================================================
        MODULE: BCFTOOLS merging

https://samtools.github.io/bcftools/bcftools.html#merge
=============================================================================
*/



process BCFTOOLS_MERGING {
    tag "${meta.sample}"
    publishDir "${params.OUTDIR}/post_processing/merged/${meta.sample}", mode: 'copy'
    container "${params.BCFTOOLS_CONTAINER}"

    input:
    tuple val(meta), path(vcfs), path(vcfs_tbi)

    output:
    tuple val(meta), path("${meta.sample}.merged.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.sample}.merged.vcf.gz.tbi"), emit: vcf_tbi

    script:
    def vcf_list = vcfs.collect { it.getName() }.join(' ')

    """
    bcftools merge \\
        --force-samples \\
        --merge all \\
        --output-type z \\
        --output "${meta.sample}.merged.vcf.gz" \\
        ${vcf_list}
    
    tabix -f -p vcf "${meta.sample}.merged.vcf.gz"
    """
}
