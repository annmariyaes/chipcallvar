/*
=============================================================================
        MODULE: MAFTOOLS library script
=============================================================================
*/


process MAFTOOLS {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/maftools/${meta.patient}", mode: 'copy'

    input:
    tuple val(meta), path(maf_dir)

    output:
    tuple val(meta), emit: plot

    script:
    """
    Rscript variant_callers.R ${maf_dir}
    """
}
