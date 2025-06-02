/*
=============================================================================
        MODULE: BWA-MEM2 Indexing
=============================================================================
*/

process BWAMEM2_INDEX {
    tag { new File(params.REFERENCE_GENOME).getName() }
    container "${params.BWAMEM2_CONTAINER}"    
    
    input:
    path reference

    output:
    tuple path(reference), path("*"), emit: index

    script:
    """
    # Check if index files exist in original location
    original_ref="${params.REFERENCE_GENOME}"
    reference_base=\$(basename "${reference}")

    if [ -f "\${original_ref}.0123" ] && [ -f "\${original_ref}.bwt.2bit.64" ]; then
        echo "BWA-MEM2 index files already exist, copying to work directory"
        cp "\${original_ref}".* ./ || true  # Copy all related files, ignore errors
    else
        echo "Creating BWA-MEM2 index for \${reference_base}"
        bwa-mem2 index ${reference}
    fi
    """
}
