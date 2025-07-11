/*
=============================================================================
        MODULE: MAFTOOLS library script
=============================================================================
*/


process MAFTOOLS {
    tag "${meta.id}"
    publishDir "${params.OUTDIR}/analysis", mode: 'copy'
    // container "${params.R_CONTAINER}"

    input:
    tuple val(meta), path(maf_dir)
    path(r_script)

    output:
    tuple val(meta), path("${meta.caller}_recurrent_mutations.csv"), emit: csv
    path("*.png"), emit: csv2
    path("*.png"), emit: plots

    script:
    """
    set -euo pipefail

    # Run the R script with proper arguments
    Rscript ${r_script} ${maf_dir} ${meta.id} ${meta.caller}
    """
}
