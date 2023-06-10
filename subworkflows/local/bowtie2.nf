include { BOWTIE2_BUILD as BOWTIE2_BUILD_HOST                    } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_HOST                    } from '../../modules/nf-core/bowtie2/align/main'
include { BOWTIE2_BUILD as BOWTIE2_BUILD_ORG                     } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_ORG                     } from '../../modules/nf-core/bowtie2/align/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_BOWTIE2 } from '../../modules/nf-core/subread/featurecounts/main'
include { GATHER_COUNTS as GATHER_COUNTS_BOWTIE2                 } from '../../modules/local/gather_counts'

workflow BOWTIE2_WORKFLOW {

    take:
        fasta_filter
        fasta_align
        gtf_align
        reads
        attribute

    main:

        ch_versions = Channel.empty()

        BOWTIE2_BUILD_HOST (
            [ "index", params.fasta_filter ]
        )

        BOWTIE2_ALIGN_HOST (
            reads,
            BOWTIE2_BUILD_HOST.out.index,
            true,
            true
        )
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_HOST.out.versions)

        BOWTIE2_BUILD_ORG (
            [ "index", params.fasta_align ]
        )
        ch_versions = ch_versions.mix(BOWTIE2_BUILD_ORG.out.versions)

        BOWTIE2_ALIGN_ORG (
            BOWTIE2_ALIGN_HOST.out.fastq,
            BOWTIE2_BUILD_ORG.out.index,
            true,
            true
        )
        ch_versions = ch_versions.mix(BOWTIE2_ALIGN_ORG.out.versions)

        SUBREAD_FEATURECOUNTS_BOWTIE2 (
            BOWTIE2_ALIGN_ORG.out.aligned.map{ [ it[0], it[1], params.gtf_align ] }, attribute
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_BOWTIE2.out.versions)

        GATHER_COUNTS_BOWTIE2(
            SUBREAD_FEATURECOUNTS_BOWTIE2.out.counts.collect{it[1]}
        )

    emit:

        reads_per_gene = SUBREAD_FEATURECOUNTS_BOWTIE2.out.counts
        log_bowtie2 = BOWTIE2_ALIGN_ORG.out.log
        counts_table = GATHER_COUNTS_BOWTIE2.out.count_table_txt
        versions = ch_versions

}
