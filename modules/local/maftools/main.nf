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
    tuple val(meta), path("*.png"), emit: plots, optional: true

    script:
    """
    Rscript ${projectDir}/scripts/variant_callers.R ${maf_dir}
    """
}
