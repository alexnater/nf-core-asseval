//
// Read files in bam folder based on reads and reference input
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow START_FROM_BAM {

    take:
    ch_reads        // channel: [ meta, fastqs ]
    ch_fasta_fai    // channel: [ meta, fasta, fai ]
    bam_folder      // Path to bam files
    ch_bam_folders  // channel: [ meta, bam_dirs ]

    main:

    ch_versions = Channel.empty()

    // Get bam files from bam folder
    Channel.fromFilePairs("${bam_folder}/bwa/*/*_{bam,bai}").view()
    Channel.fromFilePairs("${bam_folder}/minimap2/*/*_{bam,bai}").view()

    // Group by sample and combine reads with references
    ch_reads
        .map { meta, fastqs ->
            [id: meta.sample, sample: meta.sample, type: meta.type]
        }
        .groupTuple()
        .combine(ch_fasta_fai)
        .map { meta, meta2, fasta, fai -> [ meta.sample, meta ] }
        .view()
        .set { ch_input }


    emit:
    bam_bai                  // channel: [ val(meta), path(bam), path(bai) ]
    versions = ch_versions   // channel: [ versions.yml ]
}