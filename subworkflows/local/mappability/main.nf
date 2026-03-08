/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GEM2_GEMINDEXER            } from '../../../modules/nf-core/gem2/gemindexer'
include { GEM2_GEMMAPPABILITY        } from '../../../modules/nf-core/gem2/gemmappability'
include { GEM2_GEM2BEDMAPPABILITY    } from '../../../modules/nf-core/gem2/gem2bedmappability'
include { GENMAP_INDEX               } from '../../../modules/nf-core/genmap/index'
include { GENMAP_MAP                 } from '../../../modules/nf-core/genmap/map'

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

    // Combine output with index:
    def ch_to_bed = GEM2_GEMMAPPABILITY.out.map
        .join(ch_to_map.index, failOnDuplicate:true, failOnMismatch:true)
        .multiMap { meta, map, index ->
            map:   [ meta, map ]
            index: [ meta, index ]
        }
    
    //
    // MODULE: Run gem2_gem2bedmappability
    //
    GEM2_GEM2BEDMAPPABILITY (
        ch_to_bed.map,
        ch_to_bed.index
    )
    ch_versions = ch_versions.mix(GEM2_GEM2BEDMAPPABILITY.out.versions.first())

    //
    // MODULE: Run genmap_index
    //
    GENMAP_INDEX (
        ch_fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] }
    )
    ch_versions = ch_versions.mix(GENMAP_INDEX.out.versions.first())

    // Combine indices with list of kmer-sizes:
    def ch_to_genmap = GENMAP_INDEX.out.index
        .combine(kmer_sizes)
        .map { meta, index, kmer ->
            [ meta + [id: "${meta.id}_${kmer}", kmer_size: kmer], index ]
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
