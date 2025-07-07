/*
=============================================================================
        MODULE: MAFTOOLS library script
=============================================================================
*/


process MAFTOOLS {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/mafs_annotated/${meta.caller}/${meta.id}", mode: 'copy'
    container "${params.MAFTOOLS_CONTAINER}"

    input:
    tuple val(meta), path(maf_dir)
    path(r_script)

    output:
    tuple val(meta), path("*.png"), emit: plots
    path "*.png", emit: all_plots

    script:
    """
    set -euo pipefail

    # Run the R script with proper arguments
    Rscript ${r_script} ${maf_dir} ${meta.id} ${meta.caller}
    """
}
