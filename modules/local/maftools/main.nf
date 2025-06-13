/*
=============================================================================
        MODULE: MAFTOOLS library script
=============================================================================
*/


process MAFTOOLS {
    tag "${meta.patient}"
    publishDir "${params.OUTDIR}/plots", mode: 'copy'
    container "${params.R_CONTAINER}"

    input:
    tuple val(meta), path(maf_files)
    path(r_script)

    output:
    tuple val(meta), path("*.png"), emit: plots
    path "*.png", emit: all_plots

    script:
    """
    # Make R script executable
    chmod +x ${r_script}
    
    # Run the R script with proper arguments
    Rscript ${r_script} \\
        --maf_dir . \\
        --output_dir . \\
        --patient_id ${meta.patient}
    """
}
