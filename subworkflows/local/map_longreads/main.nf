//
// Map  reads to reference and call variants for PacBio HiFi read data
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MINIMAP2_ALIGN       } from '../../../modules/nf-core/minimap2/align'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAP_LONGREADS {

    take:
    ch_reads      // channel: [ meta, fastqs ]
    ch_fasta_fai  // channel: [ meta, fasta, fai ]

    main:

    ch_versions = Channel.empty()

    // Group by sample and combine reads with references
    ch_reads
        .map { meta, fastqs ->
            def new_meta = [
                id: meta.sample,
                sample: meta.sample,
                type: meta.type,
                samples_per_type: meta.samples_per_type
            ]
            [ groupKey(new_meta, meta.runs_per_sample), fastqs[0] ]
        }
        .groupTuple()
        .combine(ch_fasta_fai)
        .multiMap { meta, fastq, meta2, fasta, fai ->
            fastq: [ meta.target + [ref: meta2.id], fastq ]
            ref:   [ meta2, fasta ]
        }
        .set { ch_input }

    //
    // MODULE: Run minimap2
    //
    MINIMAP2_ALIGN (
        ch_input.fastq,
        ch_input.ref,
        true,
        "bai",
        false,
        false        
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.first())

    // join bam files and the corresponding index files: 
    MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index, failOnDuplicate:true, failOnMismatch:true)
        .set { bam_bai }

    emit:
    bam_bai                  // channel: [ val(meta), path(bam), path(bai) ]
    versions = ch_versions   // channel: [ versions.yml ]
}