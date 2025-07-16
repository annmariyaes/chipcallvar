/*
=============================================================================
        MODULE: MAFTOOLS library script
=============================================================================
*/


process MAFTOOLS {
    tag "${caller}"
    publishDir "${params.OUTDIR}/analysis", mode: 'copy'
    // container "${params.R_CONTAINER}"

    input:
    tuple val(caller), path(maf_files)
    path(r_script)

    output:
    tuple val(caller), path("${caller}_recurrent_mutations.csv"), emit: csv
    path("*.png"), emit: csv2
    path("*.png"), emit: plots

    script:
    """
    set -euo pipefail
    # Create maf directory and copy files there
    mkdir -p maf_dir
    cp ${maf_files} maf_dir/
    ls -l maf_dir/

    # Run the R script with proper arguments
    Rscript ${r_script} maf_dir ${caller}
    """
}
