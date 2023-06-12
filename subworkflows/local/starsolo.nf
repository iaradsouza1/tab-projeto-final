/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { STAR_ALIGN }                  from '../../modules/nf-core/star/align/main'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { STAR_GENOMEGENERATE }         from '../../modules/nf-core/star/genomegenerate/main'


def multiqc_report    = []

workflow STARSOLO {
    take:
    genome_fasta
    gtf
    star_index
    reads
    save_unmapped
    seq_center
    seq_platform

    main:
    ch_versions = Channel.empty()

    assert star_index || (genome_fasta && gtf):
        "Must provide a genome fasta file ('--fasta') and a gtf file ('--gtf') if no index is given!"

    assert gtf: "Must provide a gtf file ('--gtf') for STARSOLO"

    //
    // Build STAR index if not supplied
    //
    if (!star_index) {
        STAR_GENOMEGENERATE (
        genome_fasta,gtf
    )
    star_index = STAR_GENOMEGENERATE.out.index.collect()
    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }


    //
    // Perform mapping with STAR
    //
    STAR_ALIGN(
        reads,
        star_index,
        gtf,
        save_unmapped,
        true,
        seq_center,
        seq_platform
    )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)


    emit:
    versions              = ch_versions
    index                 = star_index
    result_align          = STAR_ALIGN.out.tab
    result_unmapped       = STAR_ALIGN.out.unmapped
    for_multiqc           = STAR_ALIGN.out.log_final
}
