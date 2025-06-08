//
// Run k-mer analysis
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTK_FASTK                } from '../../../modules/local/fastk/fastk'
include { MERQURYFK_MERQURYFK        } from '../../../modules/nf-core/merquryfk/merquryfk'
include { GENESCOPEFK                } from '../../../modules/nf-core/genescopefk'
include { SMUDGEPLOT                 } from '../../../modules/local/smudgeplot'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RUN_KMER_FK {
    take:
    ch_reads      // channel: [ meta, fastq ]
    ch_fasta_fai  // channel: [ meta, fasta, fai ]
    ch_kmers      // channel: [ meta, db ]
    kmer_size     // value: k-mer size

    main:
    ch_versions = Channel.empty()

    // Group reads by sample and join with pre-existing fastk databases
    ch_reads
        .map { meta, fastq ->
            def new_meta = [
                id: meta.sample,
                sample: meta.sample,
                type: meta.type,
                samples_per_type: meta.samples_per_type
            ]
            [ meta.sample, new_meta, fastq[0] ]
        }
        .groupTuple(by: [0, 1])
        .join (
            ch_kmers
                .filter { meta, db -> meta.type == 'fastk' && meta.kmer_size == kmer_size }
                .map { meta, db -> [ meta.sample, db ] },
            failOnDuplicate: true,
            remainder: true
        )
        .filter { sample, meta, fastqs, db -> !db }
        .map { sample, meta, fastqs, db -> [ meta, fastqs ] }
        .set { ch_to_count }
    
    //
    // MODULE: Run fastk by sample
    //
    FASTK_FASTK (
        ch_to_count,
        kmer_size
    )
    ch_versions = ch_versions.mix(FASTK_FASTK.out.versions.first())

    // Join outputs of fastk
    FASTK_FASTK.out.ktab
        .join(FASTK_FASTK.out.data, failOnDuplicate: true, failOnMismatch: true)
        .join(FASTK_FASTK.out.hist, failOnDuplicate: true, failOnMismatch: true)
        .join(FASTK_FASTK.out.txt, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, ktab, data, hist, txt ->
            [ meta + [kmer_size: kmer_size], ktab, data, hist, txt ]
        }
        .set { ch_fastk }

    //
    // MODULE: Run genescope.fk
    //
    GENESCOPEFK (
        ch_fastk.map { meta, ktab, data, hist, txt -> [ meta, txt ] }
    )
    ch_versions = ch_versions.mix(GENESCOPEFK.out.versions.first())

    // Combine fastk ktabs with reference file
    ch_fastk
        .map { meta, ktab, data, hist, txt -> [ meta.sample, meta, ktab, data, hist, txt ] }
        .combine(ch_fasta_fai.map { meta, fasta, fai -> [ meta.sample, meta, fasta, fai ] }, by: 0)
        .map { sample, meta, ktab, data, hist, txt, meta2, fasta, fai ->
            [ meta + [ref: meta2.id], hist, ktab, data, fasta, [] ]
        }
        .set { ch_to_merqury }

    //
    // MODULE: Run merqury
    //
    MERQURYFK_MERQURYFK (
        ch_to_merqury,
        [[:], [], []],
        [[:], [], []]
    )
    ch_versions = ch_versions.mix(MERQURYFK_MERQURYFK.out.versions.first())

    //
    // MODULE: Run smudgeplot
    //
    SMUDGEPLOT (
        ch_fastk.map { meta, ktab, data, hist, txt -> [ meta, ktab, data ] },
        4
    )
    ch_versions = ch_versions.mix(SMUDGEPLOT.out.versions.first())

    emit:
    summary = GENESCOPEFK.out.summary       // channel: [ meta, summary ]
    stats = MERQURYFK_MERQURYFK.out.stats   // channel: [ meta, stats ]
    qv = MERQURYFK_MERQURYFK.out.qv         // channel: [ meta, qv ]
    versions = ch_versions                  // channel: [ versions.yml ]
}