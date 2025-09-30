/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLAIR3                     } from '../../../modules/local/clair3'
include { TABIX_TABIX                } from '../../../modules/nf-core/tabix/tabix'
include { DEEPVARIANT_RUNDEEPVARIANT } from '../../../modules/nf-core/deepvariant/rundeepvariant'
include { GLNEXUS                    } from '../../../modules/local/glnexus'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT_CALLING {

    take:
    ch_bam_bai     // channel: [ meta, bam, bai ]
    ch_fasta_fai   // channel: [ meta, fasta, fai ]
    bed_file       // BED file with genomic intervals
    model_file     // Clair3 model file
    config_file    // GL Nexus config file

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
            bam_bai:   [ meta, bam, bai ]
            fasta_fai: [ meta2, fasta, fai ]
        }
        .set { ch_mapped }

    //
    // MODULE: Run Clair3
    //
    CLAIR3 (
        ch_mapped.bam_bai.map { meta, bam, bai ->
            [ meta, bam, bai, bed_file ]
        },
        ch_mapped.fasta_fai,
        ch_mapped.bam_bai.map { meta, bam, bai ->
            meta.type == 'hifi' ? 'hifi' : meta.type == 'ont' ? 'ont' : 'ilmn'
        },
        model_file
    )
        .gvcf
        .join(ch_mapped.bam_bai, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, bam, bai ->
            def new_meta = [
                id: "joint_${meta.ref}",
                type: meta.type,
                ref: meta.ref,
                caller: 'clair3'
            ]
            [ groupKey(new_meta, meta.samples_per_type), [ meta.sample, gvcf, bam, bai ] ]
        }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples -> 
            def (samples, gvcfs, bams, bais) = tuples.transpose()
            [ meta.target + [samples: tuple(samples)], gvcfs, bams, bais ]
        }
        .set { ch_from_clair3 }
    ch_versions = ch_versions.mix(CLAIR3.out.versions.first())

    //
    // MODULE: tabix_tabix
    //
    TABIX_TABIX (
        CLAIR3.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    CLAIR3.out.vcf
        .join(TABIX_TABIX.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, vcf, tbi ->
            def new_meta = [
                id: "joint_${meta.ref}",
                type: meta.type,
                ref: meta.ref,
                caller: 'clair3'
            ]
            [ groupKey(new_meta, meta.samples_per_type), [ meta.sample, vcf, tbi ] ]
        }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples -> 
            def (samples, vcfs, tbis) = tuples.transpose()
            [ meta.target + [samples: tuple(samples)], vcfs, tbis ]
        }
        .set { ind_vcf_tbi }

/*
    //
    // MODULE: Run DeepVariant
    //
    DEEPVARIANT_RUNDEEPVARIANT (
        ch_mapped.bam_bai.map { meta, bam, bai ->
            [ meta, bam, bai, bed_file ]
        },
        ch_mapped.fasta_fai.map { meta, fasta, fai -> [ meta, fasta ] },
        ch_mapped.fasta_fai.map { meta, fasta, fai -> [ meta, fai ] },
        [ [id:'ref'], [] ],
        [ [:], [] ]
    )
        .gvcf
        .join(ch_mapped.bam_bai, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, bam, bai ->
            def new_meta = [
                id: "joint_${meta.ref}",
                type: meta.type,
                ref: meta.ref,
                caller: 'deepvariant'
            ]
            [ groupKey(new_meta, meta.samples_per_type), [ meta.sample, gvcf, bam, bai ] ]
        }
        .groupTuple(sort: { a, b -> a[0] <=> b[0] })
        .map { meta, tuples -> 
            def (samples, gvcfs, bams, bais) = tuples.transpose()
            [ meta.target + [samples: tuple(samples)], gvcfs, bams, bais ]
        }
        .set { ch_from_dv }
    ch_versions = ch_versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions.first())
*/

    // Mix channels for joint genotyping
    ch_from_clair3
//        .mix(ch_from_dv)
        .multiMap { meta, gvcfs, bams, bais ->
            to_merge: [ meta, gvcfs ]
            bam_bai:  [ meta, bams, bais ]
        }.set { ch_calls }

/*
    //
    // MODULE: Run GL Nexus
    //
    GLNEXUS (
        ch_calls.to_merge,
        [ [id: 'regions'], bed_file ],
        ch_calls.to_merge.map { meta, gvcfs -> meta.caller == 'clair3' ? '' : 'DeepVariant' },
        ch_calls.to_merge.map { meta, gvcfs -> meta.caller == 'clair3' ? config_file : [] }
    )
    ch_versions = ch_versions.mix(GLNEXUS.out.versions.first())
*/

    emit:
//    vcf_tbi  = GLNEXUS.out.vcf_tbi          // channel: [ meta, vcf, tbi ]
    ind_vcf_tbi                             // channel: [ meta, vcf, tbi ]
    bam_bai  = ch_calls.bam_bai             // channel: [ meta, [bam], [bai] ]
    versions = ch_versions                  // channel: [ path(versions.yml) ]
}
