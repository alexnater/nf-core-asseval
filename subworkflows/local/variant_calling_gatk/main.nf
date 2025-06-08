/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_HAPLOTYPECALLER      } from '../../../modules/nf-core/gatk4/haplotypecaller'
include { GATK4_MERGEVCFS            } from '../../../modules/nf-core/gatk4/mergevcfs'
include { GATK4_COMBINEGVCFS         } from '../../../modules/nf-core/gatk4/combinegvcfs'
include { GATK4_GENOTYPEGVCFS        } from '../../../modules/nf-core/gatk4/genotypegvcfs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT_CALLING_GATK {

    take:
    ch_bam_bai     // channel: [ meta, bam, bai ]
    ch_assemblies  // channel: [ meta, fasta, fai, dict ]
    min_length     // int 

    main:

    ch_versions = Channel.empty()

    // Extract list of contigs from each assembly
    ch_assemblies
        .map { meta, fasta, fai, dict ->
            def contigs = WorkflowAssEval.getContigs(fai, min_length)
            [ meta + [ncontigs: contigs.size()], contigs, fasta, fai, dict ]
        }
        .transpose(by: 1)
        .set { ch_regions }

    // Combine bam files with their reference
    ch_bam_bai
        .map { meta, bam, bai -> [ meta.ref, meta, bam, bai ] }
        .combine(ch_regions.map { meta, contig, fasta, fai, dict ->
            [ meta.id, meta + [region: contig], fasta, fai, dict ]
            },
            by: 0)
        .multiMap { ref, meta, bam, bai, meta2, fasta, fai, dict ->
            bam_bai:   [ meta + [region: meta2.region, ncontigs: meta2.ncontigs], bam, bai, [], [] ]
            fasta:     [ meta2, fasta ]
            fai:       [ meta2, fai ]
            dict:      [ meta2, dict ]
        }
        .set { ch_mapped }

    //
    // MODULE: Run GATK HaplotypeCaller
    //
    GATK4_HAPLOTYPECALLER (
        ch_mapped.bam_bai,
        ch_mapped.fasta,
        ch_mapped.fai,
        ch_mapped.dict,
        [ [:], [] ],
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())

    // Group gvcf files by sample
    GATK4_HAPLOTYPECALLER.out.vcf
        .join(GATK4_HAPLOTYPECALLER.out.tbi, failOnDuplicate:true, failOnMismatch:true)
        .map { meta, gvcf, tbi -> 
            [ groupKey(meta + [id: 'joint', sample: 'joint'], meta.samples_per_type), gvcf, tbi ]
        }
        .groupTuple(sort: {a, b -> a.name <=> b.name})
        .map { meta, gvcfs, tbi -> [ meta.ref, meta, gvcfs, tbi ] }
        .combine (
            ch_assemblies
                .map { meta, fasta, fai, dict ->
                    [ meta.id, meta, fasta, fai, dict ]
                },
            by: 0)
        .multiMap { ref, meta, gvcfs, tbi, meta2, fasta, fai, dict ->
            gvcfs: [ meta, gvcfs, tbi, [], [] ]
            fasta: [ meta2, fasta ]
            fai:   [ meta2, fai ]
            dict:  [ meta2, dict ]
        }
        .set { ch_to_genotype }

    //
    // MODULE: Run GATK GenotypeGVCFs
    //
    GATK4_GENOTYPEGVCFS (
        ch_to_genotype.gvcfs,
        ch_to_genotype.fasta,
        ch_to_genotype.fai,
        ch_to_genotype.dict,
        [ [:], [] ],
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions.first())

    // Group vcf files by region
    GATK4_GENOTYPEGVCFS.out.vcf
        .map { meta, vcf -> [ groupKey(meta - meta.subMap(['region']), meta.ncontigs), vcf ] }
        .groupTuple()
        .map { meta, vcfs -> [ meta.ref, meta.target, vcfs ] }
        .combine (
            ch_assemblies
                .map { meta, fasta, fai, dict -> [ meta.id, meta, dict ] },
            by: 0)
        .multiMap { ref, meta, vcfs, meta2, dict ->
            vcfs: [ meta, vcfs ]
            dict: [ meta2, dict ]
        }
        .set { ch_to_concat }

    //
    // MODULE: Run GATK MergeVcfs
    //
    GATK4_MERGEVCFS (
        ch_to_concat.vcfs,
        ch_to_concat.dict
    ).vcf
     .join(GATK4_MERGEVCFS.out.tbi, failOnDuplicate:true, failOnMismatch:true)
     .set { vcf_tbi }
    ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions.first())

    emit:
    vcf_tbi                     // channel: [ meta, vcf, tbi ]
    versions = ch_versions      // channel: [ path(versions.yml) ]
}
