//
// Map Hi-C reads to reference and produce PretextView file
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP                      } from '../../../modules/nf-core/fastp'
include { BWA_MAP_FILTER             } from '../../../modules/local/bwa/map_filter'
include { GATK4_MARKDUPLICATES       } from '../../../modules/nf-core/gatk4/markduplicates'
include { PRETEXTMAP                 } from '../../../modules/nf-core/pretextmap'
include { PRETEXTSNAPSHOT            } from '../../../modules/nf-core/pretextsnapshot'
include { YAHS                       } from '../../../modules/nf-core/yahs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAP_HIC {

    take:
    ch_reads        // channel: [ meta, fastqs ]
    ch_reference    // channel: [ meta, fasta, fai, index ]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run fastp
    //
    FASTP (
        ch_reads,
        [],
        false,
        false,
        false
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

/*
    // prepare channel for reference genomes:
    ch_fasta_fai
        .combine(ch_reads.first().map { true })     // this prevents building the index if the reads channel is empty.
        .map { meta, fasta, fai, trigger -> [ meta, fasta ] }
        .set { ch_to_index }

    // Create BWA index:
    BWA_INDEX(ch_to_index)
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())
*/

    // Prepare channel that joins reads with reference inidices:
    FASTP.out.reads
        .combine(ch_reference.map { meta, fasta, fai, index -> [ meta, index ] })
        .multiMap { meta, reads, meta2, index ->
            reads: [ meta + [ref: meta2.id], reads ]
            index: [ meta2, index ]
        }
        .set { ch_to_map }

    // Map reads with BWA:
    BWA_MAP_FILTER (
        ch_to_map.reads,
        ch_to_map.index,
        params.mapq_filter
    )
    ch_versions = ch_versions.mix(BWA_MAP_FILTER.out.versions.first())

    // group entries by sample:
    BWA_MAP_FILTER.out.bam
        .map { meta, bam ->
            def new_meta = [
                id: meta.sample,
                sample: meta.sample,
                type: meta.type,
                ref: meta.ref,
                samples_per_type: meta.samples_per_type
            ]
            [ groupKey(new_meta, meta.runs_per_sample), bam ]
        }
        .groupTuple(sort: true)
        .set { ch_bam_bysample }

    // combine bam files with fasta reference and fasta index:
    ch_bam_bysample
        .map { meta, bams -> [ meta.ref, meta.target, bams ] }
        .combine(ch_reference.map { meta, fasta, fai, index ->
            [ meta.id, fasta, fai ]
            },
            by: 0)
        .multiMap { id, meta, bams, fasta, fai ->
            bams:  [ meta, bams ]
            fasta: fasta
            fai:   fai
        }.set { ch_to_dedup }
        
    // merge bam files by sample and mark duplicates:
    GATK4_MARKDUPLICATES (
        ch_to_dedup.bams,
        ch_to_dedup.fasta,
        ch_to_dedup.fai
    )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())

    // join bam files and the corresponding index files: 
    GATK4_MARKDUPLICATES.out.bam
        .join(GATK4_MARKDUPLICATES.out.bai, failOnDuplicate:true, failOnMismatch:true)
        .set { ch_bam_bai }

/*  THIS DOESN'T WORK
    // join all hic reads:
    ch_bam_bai
        .map { meta, bam, bai -> [ [id: meta.ref], bam, bai ] }
        .groupTuple()
        .set { ch_bams_by_ref }
*/

    //
    // MODULE: Run pretextmap
    //
    PRETEXTMAP (
        ch_bam_bai.map { meta, bam, bai -> [ meta, bam ] },
        [ [:], [], [] ]
    )
    ch_versions = ch_versions.mix(PRETEXTMAP.out.versions.first())

    //
    // MODULE: Run pretextsnapshot
    //
    PRETEXTSNAPSHOT (
        PRETEXTMAP.out.pretext
    )
    ch_versions = ch_versions.mix(PRETEXTSNAPSHOT.out.versions.first())

    // combine bam files with fasta reference and fasta index:
    ch_bam_bai
        .map { meta, bam, bai -> [ meta.ref, meta, bam, bai ] }
        .combine(ch_reference.map { meta, fasta, fai, index ->
            [ meta.id, fasta, fai ]
            },
            by: 0)
        .multiMap { id, meta, bam, bai, fasta, fai ->
            bam: [ meta, bam ]
            fasta:   fasta
            fai:     fai
        }.set { ch_to_yahs }

    //
    // MODULE: Run yahs
    //
    YAHS (
        ch_to_yahs.bam,
        ch_to_yahs.fasta,
        ch_to_yahs.fai
    )
    ch_versions = ch_versions.mix(YAHS.out.versions.first())

    emit:
    bam_bai  = ch_bam_bai                 // channel: [ val(meta), path(bam), path(bai) ]
    pretext  = PRETEXTMAP.out.pretext     // channel: [ val(meta), path(pretext) ]
    versions = ch_versions                // channel: [ versions.yml ]
}