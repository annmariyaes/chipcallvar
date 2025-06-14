/*
=============================================================================
        MODULE: MULTIQC
=============================================================================
*/


process MULTIQC {
    publishDir "${params.OUTDIR}/multiqc", mode: 'copy'
    container "${params.MULTIQC_CONTAINER}"
    
    input:
    path(fastqc_files)
    path(samtools_files)
    path(mosdepth_files)
    path(bcftools_files)
    path(vep_files)
    path(config_file)
    
    output:
    path "multiqc_report.html", emit: html
    path "multiqc_data/", emit: data, optional: true

    script:
    """
    multiqc . \\
        --filename multiqc_report.html \\
        --force \\
        --config ${config_file} \\
        --verbose
    """
}
