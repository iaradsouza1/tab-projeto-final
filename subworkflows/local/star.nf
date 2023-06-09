include { STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_HOST     } from '../../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN as STAR_ALIGN_HOST                       } from '../../modules/nf-core/star/align/main'
include { STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_ORG      } from '../../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN as STAR_ALIGN_ORG                        } from '../../modules/nf-core/star/align/main'
include { SUBREAD_FEATURECOUNTS as SUBREAD_FEATURECOUNTS_STAR } from '../../modules/nf-core/subread/featurecounts/main'
include { GATHER_COUNTS as GATHER_COUNTS_STAR                 } from '../../modules/local/gather_counts'

workflow STAR_WORKFLOW {

    take:

        fasta_filter
        gtf_filter
        fasta_align
        gtf_align
        reads
        attribute

    main:

        ch_versions = Channel.empty()

        STAR_GENOMEGENERATE_HOST(
            fasta_filter,
            gtf_filter
        )

        STAR_ALIGN_HOST(
            reads,
            STAR_GENOMEGENERATE_HOST.out.index,
            gtf_filter,
            false,
            '',
            '',
        )
        ch_versions = ch_versions.mix(STAR_ALIGN_HOST.out.versions)

        STAR_GENOMEGENERATE_ORG(
            fasta_align,
            gtf_align
        )

        STAR_ALIGN_ORG(
            STAR_ALIGN_HOST.out.fastq,
            STAR_GENOMEGENERATE_ORG.out.index,
            gtf_align,
            false,
            '',
            '',
        )
        ch_versions = ch_versions.mix(STAR_ALIGN_ORG.out.versions)

        SUBREAD_FEATURECOUNTS_STAR(
            STAR_ALIGN_ORG.out.bam_sorted.map{ [ it[0], it[1], params.gtf_align ] }, attribute
        )
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS_STAR.out.versions)

        GATHER_COUNTS_STAR(
            SUBREAD_FEATURECOUNTS_STAR.out.counts.collect{ it[1] }
        )

    emit:

        reads_per_gene = STAR_ALIGN_ORG.out.read_per_gene_tab
        log_star = STAR_ALIGN_ORG.out.log_final
        counts_table = GATHER_COUNTS_STAR.out.count_table_txt
        versions = ch_versions

}
