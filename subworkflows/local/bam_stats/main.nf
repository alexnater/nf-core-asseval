//
// Create depth and coverage statistics for each merged BAM file
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FLAGSTAT                    } from '../../../modules/nf-core/samtools/flagstat'
include { MOSDEPTH                             } from '../../../modules/nf-core/mosdepth'
include { PANDEPTH                             } from '../../../modules/local/pandepth'
include { SUMMARIZE_STATS                      } from '../../../modules/local/summarize_stats'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BAM_STATS {
    take:
    ch_bam_bai     // channel (mandatory): [ val(meta), path(bam), path(bai) ]
    ch_fasta_fai   // channel (mandatory): [ val(meta), path(fasta), path(fai) ]
    bed_file       // file (ptional)

    main:
    ch_versions = Channel.empty()

    // Combine bam files with their reference
    ch_bam_bai
        .map { meta, bam, bai -> [ meta.ref, meta, bam, bai ] }
        .combine(ch_fasta_fai.map { meta, fasta, fai ->
            [ meta.id, meta, fasta, fai ]
            },
            by: 0)
        .multiMap { ref, meta, bam, bai, meta2, fasta, fai ->
            bam_bai: [ meta, bam, bai ]
            fasta:   [ meta2, fasta ]
        }
        .set { ch_mapped }

    // run SAMtools flagstat per merged bam file:
    SAMTOOLS_FLAGSTAT(ch_mapped.bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.first())

    // run mosdepth per merged bam file:
    MOSDEPTH (
        ch_mapped.bam_bai.map { meta, bam, bai -> [ meta, bam, bai, bed_file ] },
        ch_mapped.fasta
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())

    // run pandepth per merged bam file:
    PANDEPTH (
        ch_mapped.bam_bai.map { meta, bam, bai -> [ meta, bam, bai, bed_file ] },
        ch_mapped.fasta,
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(PANDEPTH.out.versions.first())

    // Summarize all bam stats:
    SAMTOOLS_FLAGSTAT.out.flagstat
        .join(MOSDEPTH.out.summary_txt, failOnDuplicate:true, failOnMismatch:true)
        .join(MOSDEPTH.out.global_txt, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, stat, depth, cov ->
            def new_meta = [
                id: 'summary',
                type: meta.type,
                ref: meta.ref
            ]
            [ groupKey(new_meta, meta.samples_per_type), stat, depth, cov, 0, 0 ]
          }
        .groupTuple(sort: true)
        .set { ch_to_summary }

    SUMMARIZE_STATS(ch_to_summary)
    ch_versions = ch_versions.mix(SUMMARIZE_STATS.out.versions)

    emit:
    depth = MOSDEPTH.out.summary_txt     // channel: [ meta, depth ]
    coverage = MOSDEPTH.out.global_txt   // channel: [ meta, coverage ]
    versions = ch_versions               // channel: [ versions.yml ]
}