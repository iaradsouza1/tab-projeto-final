/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
include { STAR_ALIGN }                  from '../../modules/nf-core/star/align/main'

/* --    IMPORT NF-CORE MODULES/SUBWORKFLOWS   -- */
include { STAR_GENOMEGENERATE }         from '../../modules/nf-core/star/genomegenerate/main'


def multiqc_report    = []

workflow STARSOLO {
    take:
    genome_fasta
    gtf
    reads
    save_unmapped
    star_index
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
    }
        ch_idx = STAR_GENOMEGENERATE.out.index.collect()

    //
    // Perform mapping with STAR
    //
    STAR_ALIGN(
        reads,
        ch_idx,
        gtf,
        save_unmapped,
        true,
        seq_center,
        seq_platform
    )
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions, STAR_GENOMEGENERATE.out.versions)


    emit:
    versions              = ch_versions
    index                 = ch_idx
    result_align          = STAR_ALIGN.out.tab
    result_unmapped       = STAR_ALIGN.out.unmapped
    for_multiqc           = STAR_ALIGN.out.log_final
}
