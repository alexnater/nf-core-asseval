//
// Generate report with depth per site
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_DEPTH                       } from '../../../modules/nf-core/samtools/depth'
include { TABIX_BGZIPTABIX as BGZIPTABIX_DEPTH } from '../../../modules/nf-core/tabix/bgziptabix'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BAM_DEPTH {
    take:
    ch_bam_bai        // channel (mandatory): [ val(meta), path(bam), path(bai) ]

    main:
    ch_versions = Channel.empty()

    // Run SAMtools depth over all merged bam files per reference genome:
    ch_bam_bai
        .map { meta, bam, bai ->
            def new_meta = [
                id: "depth_${meta.ref}",
                type: meta.type,
                ref: meta.ref
            ]
            [ groupKey(new_meta, meta.samples_per_type), bam ]
        }
        .groupTuple(sort: {a, b -> a.name <=> b.name})
        .map { gkey, bams -> [ gkey.target, bams ] }
        .set { ch_to_depth }

    SAMTOOLS_DEPTH (
        ch_to_depth,
        [[], []]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DEPTH.out.versions.first())

    // bgzip and tabix the depth file:
    BGZIPTABIX_DEPTH(SAMTOOLS_DEPTH.out.tsv)
        .gz_tbi
        .set { depth }
    ch_versions = ch_versions.mix(BGZIPTABIX_DEPTH.out.versions.first())

    emit:
    depth                         // channel: [ val(meta), path(depth.tsv.gz), path(depth.tsv.tbi) ]
    versions = ch_versions        // channel: [ versions.yml ]
}