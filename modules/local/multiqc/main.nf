
process MULTIQC {
    publishDir "${params.OUTDIR}/multiqc", mode: 'copy'
    container "${params.MULTIQC_CONTAINER}"
    
    input:
    path('fastqc/*')
    path('samtools/*')
    path('bcftools/*') 
    path('vep/*')
    
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
