process MULTIQC {
    label 'process_single'

    publishDir "${params.OUTDIR}/multiqc", mode: 'copy'
    container "${params.MULTIQC_CONTAINER}"

    input:
    path(reports)

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data
    path "versions.yml",        emit: versions

    script:
    def args = task.ext.args ?: ''
    def reportFiles = reports.collect { it.toString() }.join(' ')

    """
    multiqc \\
        $args \\
        --filename multiqc_report.html \\
        $reportFiles

    echo \\"${task.process}:\\" > versions.yml
    echo \\"  multiqc: \$(multiqc --version | sed 's/multiqc, version //')\\" >> versions.yml
    """
}
