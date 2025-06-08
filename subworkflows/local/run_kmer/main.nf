//
// Run k-mer analysis
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MERYL_COUNT                } from '../../../modules/nf-core/meryl/count'
include { MERYL_UNIONSUM             } from '../../../modules/nf-core/meryl/unionsum'
include { MERYL_HISTOGRAM            } from '../../../modules/local/meryl/histogram'
include { MERQURY_MERQURY            } from '../../../modules/nf-core/merqury/merqury'
include { GENOMESCOPE2               } from '../../../modules/nf-core/genomescope2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_KMER {
    take:
    ch_reads      // channel: [ meta, fastq ]
    ch_fasta_fai  // channel: [ meta, fasta, fai ]
    ch_kmers      // channel: [ meta, db ]
    kmer_size     // value: k-mer size

    main:
    ch_versions = Channel.empty()

    // Join reads with pre-existing meryls databases
    ch_reads
        .map { meta, fastq -> [ meta.sample, meta, fastq ] }
        .groupTuple()
        .join (
            ch_kmers
                .filter { meta, db -> meta.type == 'meryl' && meta.kmer_size == kmer_size }
                .map { meta, db -> [ meta.sample, db ] },
            failOnDuplicate: true,
            remainder: true
        )
        .transpose()
        .filter { sample, meta, fastq, db -> !db }
        .map { sample, meta, fastq, db -> [ meta, fastq ] }
        .set { ch_to_count }
    
    //
    // MODULE: Run meryl count
    //
    MERYL_COUNT (
        ch_to_count,
        kmer_size
    )
    ch_versions = ch_versions.mix(MERYL_COUNT.out.versions.first())

    // Group by sample
    MERYL_COUNT.out.meryl_db
        .map { meta, db ->
            def new_meta = [
                id: meta.sample,
                sample: meta.sample,
                type: meta.type,
                samples_per_type: meta.samples_per_type
            ]
            [ groupKey(new_meta, meta.runs_per_sample), db ]
        }
        .groupTuple()
        .map { meta, dbs -> [ meta.target, dbs ] }
        .set { ch_to_union }

    //
    // MODULE: Run meryl unionsum
    //
    MERYL_UNIONSUM (
        ch_to_union,
        kmer_size
    )
    ch_versions = ch_versions.mix(MERYL_UNIONSUM.out.versions.first())

    // Mix newly generated databases with existing databases
    ch_kmers
        .mix(MERYL_UNIONSUM.out.meryl_db.map { meta, db ->
            [ meta + [kmer_size: kmer_size], db ]
        })
        .set { ch_meryl_db }

    //
    // MODULE: Run meryl histogram
    //
    MERYL_HISTOGRAM (
        ch_meryl_db
    )
    ch_versions = ch_versions.mix(MERYL_HISTOGRAM.out.versions.first())

    //
    // MODULE: Run genomescope2
    //
    GENOMESCOPE2 (
        MERYL_HISTOGRAM.out.hist
    )
    ch_versions = ch_versions.mix(GENOMESCOPE2.out.versions.first())

    // Combine merly dbs with reference file
    ch_meryl_db
        .map { meta, db -> [ meta.sample, meta, db ] }
        .combine(ch_fasta_fai.map { meta, fasta, fai -> [ meta.sample, meta, fasta, fai ] }, by: 0)
        .map { sample, meta, db, meta2, fasta, fai -> [ meta + [ref: meta2.id], db, fasta ] }
        .set { ch_to_merqury }

    //
    // MODULE: Run merqury
    //
    MERQURY_MERQURY (
        ch_to_merqury
    )
    ch_versions = ch_versions.mix(MERQURY_MERQURY.out.versions.first())

    emit:
    versions = ch_versions        // channel: [ versions.yml ]
}