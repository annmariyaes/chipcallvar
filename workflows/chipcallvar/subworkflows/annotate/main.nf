/*
=============================================================================
    SUBWORKFLOW: Variant Annotation
=============================================================================
*/


include { ENSEMBL_VEP } from '../../../../modules/local/ensembl-vep/main.nf'


workflow VARIANT_ANNOTATION {
    take:
    ch_vep // channel: [ meta, vcf ]
    
    main:    
    vep = ENSEMBL_VEP(ch_vep)
    
    emit:
    vcf = vep.vcf  // [ meta, annotated_vcf ]
    vep_stats = vep.vep_stats
}
