// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join

include { GUNZIP as GUNZIP_IDX_HOST     } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_IDX_ORG      } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF          } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GEN          } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GEN1         } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF1         } from '../../modules/nf-core/gunzip/main'


workflow PREPARE_GENOME {

    take:
    fasta_filter
    gtf_filter
    fasta_align
    gtf_align
    // star_org_index
    // star_host_index



    main:
    ch_versions = Channel.empty()

    //
    //Star Index unzip: Optional
    //

    //     if(star_org_index != null) {

    //         if(star_org_index.endsWith('.gz')) {
    //         star_org_index         = GUNZIP_IDX_ORG ( [ [:], star_org_index ] ).gunzip.map { it[1] }
    //         ch_versions        = ch_versions.mix(GUNZIP_IDX_ORG.out.versions)
    //     } else {
    //         star_org_index      = Channel.from( file(star_org_index) )
    //     }
    // }

    // if(star_host_index != null) {
    //     if(star_host_index.endsWith('.gz')) {
    //         star_host_index         = GUNZIP_IDX_HOST ( [ [:], star_host_index ] ).gunzip.map { it[1] }
    //         ch_versions        = ch_versions.mix(GUNZIP_IDX_HOST.out.versions)
    //     } else {
    //         star_host_index      = Channel.from( file(star_host_index) )
    //     }
    // }

    //Host data unzip
    //
    if (gtf_filter.endsWith('.gz')) {
        ch_gtf_filter      = GUNZIP_GTF ( [ [:], gtf_filter ] ).gunzip.map { it[1] }
        ch_versions        = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf_filter = Channel.from( file(gtf_filter) )
    }

    if (fasta_filter.endsWith('.gz')) {
        ch_fasta_filter      = GUNZIP_GEN ( [ [:], fasta_filter ] ).gunzip.map { it[1] }
        ch_versions          = ch_versions.mix(GUNZIP_GEN.out.versions)
    } else {
        ch_fasta_filter      = file(fasta_filter)
    }

    //
    //Target organism data unzip
    //
    if (gtf_align.endsWith('.gz')) {
        ch_gtf_align         = GUNZIP_GTF1 ( [ [:], gtf_align ] ).gunzip.map { it[1] }
        ch_versions          = ch_versions.mix(GUNZIP_GTF1.out.versions)
    } else {
        ch_gtf_align = file(gtf_align)
    }

    if (fasta_align.endsWith('.gz')) {
        ch_fasta_align       = GUNZIP_GEN1 ( [ [:], fasta_align ] ).gunzip.map { it[1] }
        ch_versions          = ch_versions.mix(GUNZIP_GEN1.out.versions)
    } else {
        ch_fasta_align = Channel.from( file(fasta_align) )
    }

    emit:
    fasta_fil         = ch_fasta_filter
    gtf_fil           = ch_gtf_filter
    fasta_al          = ch_fasta_align
    gtf_al            = ch_gtf_align
    versions          = ch_versions
    // star_org_index
    // star_host_index
}
