process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    container 'ewels/multiqc'

    input:
    path fastqc_reports
    path samtools_stats
    path bcftools_stats
    path vep_stats

    output:
    path("multiqc_report.html")
    path("multiqc_data")

    script:
    """
    mkdir multiqc_input
    cp $fastqc_reports multiqc_input/
    cp $samtools_stats multiqc_input/
    cp $bcftools_stats multiqc_input/
    cp $vep_stats multiqc_input/

    multiqc multiqc_input
    """
}
