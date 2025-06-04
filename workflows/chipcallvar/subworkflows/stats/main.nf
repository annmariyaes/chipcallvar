include { BCFTOOLS_STATS } from '../../../../modules/local/bcftools/stats'


/*
=============================================================================
        SUBWORKFLOW: VCF statistics without vep annotation
=============================================================================
*/


workflow VCF_STATS {
    take:
    ch_vep

    main:
    BCFTOOLS_STATS(ch_vep)

}
