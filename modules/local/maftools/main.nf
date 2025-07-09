/*
=============================================================================
        MODULE: MAFTOOLS library script
=============================================================================
*/


process MAFTOOLS {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/analysis", mode: 'copy'
    container "${params.R_CONTAINER}"

    input:
    tuple val(meta), path(maf_dir)
    path(r_script)

    output:
    tuple val(meta), path("*.png"), emit: plots
    path("*.csv"), emit: csv

    script:
    """
    set -euo pipefail

    # Run the R script with proper arguments
    Rscript --vanilla -e '.libPaths("/usr/local/lib/R/site-library")' ${r_script} ${maf_dir} ${meta.id} ${meta.caller}
    """
}
