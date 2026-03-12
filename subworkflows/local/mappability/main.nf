/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GEM2_GEMINDEXER            } from '../../../modules/local/gem2/gemindexer'
include { GEM2_GEMMAPPABILITY        } from '../../../modules/local/gem2/gemmappability'
include { GENMAP_INDEX               } from '../../../modules/local/genmap/index'
include { GENMAP_MAP                 } from '../../../modules/local/genmap/map'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAPPABILITY {

    take:
    ch_fasta_fai   // channel: [ meta, fasta, fai ]
    kmer_sizes     // list

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Run gem2_gemindexer
    //
    GEM2_GEMINDEXER (
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] }
    )
    ch_versions = ch_versions.mix(GEM2_GEMINDEXER.out.versions.first())

    // Combine indices with list of kmer-sizes:
    def ch_to_map = GEM2_GEMINDEXER.out.index
        .combine(kmer_sizes)
        .multiMap { meta, index, kmer ->
            index: [ meta + [id: "${meta.id}_${kmer}"], index ]
            kmer:  kmer
        }

    //
    // MODULE: Run gem2_gemmappability
    //
    GEM2_GEMMAPPABILITY (
        ch_to_map.index,
        ch_to_map.kmer
    )
    ch_versions = ch_versions.mix(GEM2_GEMMAPPABILITY.out.versions.first())

    //
    // MODULE: Run genmap_index
    //
    GENMAP_INDEX (
        ch_fasta_fai
    )
    ch_versions = ch_versions.mix(GENMAP_INDEX.out.versions.first())

    // Join indices with chromsizes and combine with list of kmer-sizes:
    def ch_to_genmap = GENMAP_INDEX.out.index
        .join(GENMAP_INDEX.out.sizes, failOnDuplicate:true, failOnMismatch:true)
        .combine(kmer_sizes)
        .map { meta, index, sizes, kmer ->
            [ meta + [id: "${meta.id}_${kmer}", kmer_size: kmer], index, sizes ]
        }

    //
    // MODULE: Run genmap_map
    //
    GENMAP_MAP (
        ch_to_genmap,
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(GENMAP_MAP.out.versions.first())

    emit:
    versions = ch_versions            // channel: [ path(versions.yml) ]
}
