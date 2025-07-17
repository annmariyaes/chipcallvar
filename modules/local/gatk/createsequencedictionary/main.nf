/*
=============================================================================
        MODULE: SEQUENCE DICT
=============================================================================
*/

process GATK_CREATE_SEQUENCE_DICTONARY {
    container "${params.GATK_CONTAINER}"

    input:
    path reference

    output:
    path "*.dict", emit: dict

    script:
    """
    # Check if index files exist in original location
    original_ref="${params.REFERENCE_GENOME}"
    reference_base=\$(basename "${reference}")

    if [ -f "\${original_ref}.dict" ]; then
        echo "DICT file already exist, copying to work directory"
        cp "\${original_ref}".* ./ || true  # Copy all related files, ignore errors
    else
        echo "Creating sequence dictionary for \${reference_base}"
        gatk CreateSequenceDictionary -R ${reference}
    fi
    """
}
