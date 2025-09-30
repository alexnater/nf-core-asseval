//
// Generate index files for assemblies
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMTOOLS_FAIDX             } from '../../../modules/nf-core/samtools/faidx'
include { SAMTOOLS_DICT              } from '../../../modules/nf-core/samtools/dict'
include { BWA_INDEX                  } from '../../../modules/nf-core/bwa/index'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_GENOMES {
    take:
    ch_assemblies        // channel (mandatory): [ val(meta), path(fasta), path(gtf), path(yaml) ]
    ch_reads             // channel (mandatory): [ val(meta), path(fastqs) ]

    main:
    ch_versions = Channel.empty()

    // Check if fasta index is already present:
    ch_assemblies
        .branch { meta, fasta, gtf ->
            def fai = file("${fasta}.fai")
            with_fai: fai.exists()
                return [ meta, fasta, fai ]
            no_fai: true
                return [ meta, fasta ]
        }
        .set { ch_to_faidx }

    //
    // MODULE: Run samtools faidx
    //
    SAMTOOLS_FAIDX (
        ch_to_faidx.no_fai,
        [[:], []]
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    // Check if fasta dict is already present:
    ch_to_faidx.no_fai
        .join(SAMTOOLS_FAIDX.out.fai, failOnDuplicate: true, failOnMismatch: true)
        .mix(ch_to_faidx.with_fai)
        .branch { meta, fasta, fai ->
            def dict = file("${fasta.baseName}.dict")
            with_dict: dict.exists()
                return [ meta, fasta, fai, dict ]
            no_dict: true
                return [ meta, fasta, fai ]
        }
        .set { ch_to_dict }

    //
    // MODULE: Run samtools dict
    //
    SAMTOOLS_DICT (
        ch_to_dict.no_dict.map { meta, fasta, fai -> [ meta, fasta ] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions.first())

    ch_to_dict.no_dict
        .join(SAMTOOLS_DICT.out.dict, failOnDuplicate: true, failOnMismatch: true)
        .mix(ch_to_dict.with_dict)
        .join (
            ch_assemblies.map { meta, fasta, gtf -> [ meta, gtf ] },
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .set { ch_genomes }

    // Prepare channel for indexing:
    ch_assemblies
        .combine(
            ch_reads
                .filter { meta, fastqs -> meta.type == 'illumina' || meta.type == 'hic' }
                .first()
                .map { true })     // this prevents building the index if the reads channel is empty.
        .map { meta, fasta, gtf, trigger -> [ meta, fasta ] }
        .set { ch_to_index }

    //
    // MODULE: Run BWA index
    //
    BWA_INDEX (
        ch_to_index
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    ch_genomes
        .join(BWA_INDEX.out.index, failOnDuplicate: true, remainder: true)
        .set { genomes }

    emit:
    genomes                       // channel: [ val(meta), path(fasta), path(fai), path(dict), path(gtf), path(index) ]
    versions = ch_versions        // channel: [ versions.yml ]
}