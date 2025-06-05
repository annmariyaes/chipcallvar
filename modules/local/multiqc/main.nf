
process MULTIQC {
    publishDir "${params.OUTDIR}/multiqc", mode: 'copy'
    container "${params.MULTIQC_CONTAINER}"
    
    input:
    path(fastqc_files)
    path(samtools_files)
    path(bcftools_files)
    path(vep_files)
    
    output:
    path "multiqc_report.html", emit: html
    
    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --force \\
        --verbose
    """
}
