/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { QUAST                  } from '../modules/nf-core/quast'
include { BUSCO_BUSCO            } from '../modules/nf-core/busco/busco'
include { WINDOWS_STATS          } from '../modules/local/windows_stats'
include { PLOT_WINDOWS           } from '../modules/local/plot_windows'
include { GENERATE_EAR           } from '../modules/local/generate_ear'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { PREPARE_GENOMES        } from '../subworkflows/local/prepare_genomes'
include { RUN_KMER_FK            } from '../subworkflows/local/run_kmer_fk'
include { MAP_HIC                } from '../subworkflows/local/map_hic'
include { MAP_LONGREADS          } from '../subworkflows/local/map_longreads'
include { MAP_ILLUMINA           } from '../subworkflows/local/map_illumina'
include { BAM_STATS              } from '../subworkflows/local/bam_stats'
include { BAM_DEPTH              } from '../subworkflows/local/bam_depth'
include { VARIANT_CALLING        } from '../subworkflows/local/variant_calling'
include { VARIANT_CALLING_GATK   } from '../subworkflows/local/variant_calling_gatk'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_asseval_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def steps = params.steps ? params.steps.split(',') : []
if ((steps.contains('stats') || steps.contains('variant_calling')) && !steps.contains('mapping')) {
    steps << "mapping"
}
def busco_lineages_dir = file(params.busco_lineages_path, type: 'dir', checkIfExists: true)
def model_file = params.genotype_model ? file(params.genotype_model, checkIfExists: true) : []
def config_file = file(params.glnexus_config, checkIfExists: true)


workflow ASSEVAL {

    take:
    ch_assemblies // channel: samplesheet read in from --assemblies
    ch_reads      // channel: samplesheet read in from --reads
    ch_kmers      // channel: samplesheet read in from --kmers

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: prepare_genomes
    // 
    PREPARE_GENOMES ( 
        ch_assemblies.map { meta, fasta, gtf, yaml -> [ meta, fasta, gtf ] },
        ch_reads
    )
    ch_versions = ch_versions.mix(PREPARE_GENOMES.out.versions)

    if (steps.contains('qc')) {
        //
        // MODULE: Run FastQC
        //
        FASTQC (
            ch_reads
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    if (steps.contains('evaluation')) {
        
        // Group assemblies by type
        ch_assemblies
            .map { meta, fasta, gtf, yaml ->
                [ [comp: 'by_type', id: meta.type, type: meta.type], [ meta.id, fasta ] ]
            }
            .groupTuple(sort: { a, b -> a[0] <=> b[0] })
            .map { meta, tuples -> [ meta + [labels: tuples.collect { it[0] }], tuples.collect { it[1] } ] }
            .set { ch_by_type }

        // Group assemblies by sample
        ch_assemblies
            .map { meta, fasta, gtf, yaml ->
                [ [comp: 'by_sample', id: meta.sample, sample: meta.sample], [ meta.id, fasta ] ]
            }
            .groupTuple(sort: { a, b -> a[0] <=> b[0] })
            .map { meta, tuples -> [ meta + [labels: tuples.collect { it[0] }], tuples.collect { it[1] } ] }
            .set { ch_by_sample }

        //
        // MODULE: Run quast
        //
        QUAST (
            ch_by_type.mix(ch_by_sample),
            [ [:], [] ],
            [ [:], [] ]
        )
        ch_versions = ch_versions.mix(QUAST.out.versions.first())

        //
        // SUBWORKFLOW: run_kmer_fk
        //    
        RUN_KMER_FK (
            ch_reads.filter { meta, fastq -> meta.type == 'hifi' },
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai ] },
            ch_kmers,
            params.kmer_size
        )
        ch_versions = ch_versions.mix(RUN_KMER_FK.out.versions)

        //
        // MODULE: Run busco
        //
        BUSCO_BUSCO (
            ch_assemblies.map { meta, fasta, gtf, yaml -> [ meta, fasta ] },
            "genome",
            params.busco_lineage,
            busco_lineages_dir,
            []
        )
        ch_versions = ch_versions.mix(BUSCO_BUSCO.out.versions.first())
    }

    if (steps.contains('mapping')) {

        // Branch reads by read type:
        ch_reads.branch { meta, reads ->
            longreads: meta.type == 'hifi' || meta.type == 'ont'
            hic: meta.type == 'hic'
            illumina:  true
        }.set { ch_reads_bytype }

        //
        // SUBWORKFLOW: map_longreads
        // 
        MAP_LONGREADS (
            ch_reads_bytype.longreads,
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai ] }
        )
        ch_versions = ch_versions.mix(MAP_LONGREADS.out.versions)

        //
        // SUBWORKFLOW: map_hic
        //    
        MAP_HIC (
            ch_reads_bytype.hic,
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai, index ] }
        )
        ch_versions = ch_versions.mix(MAP_HIC.out.versions)

        //
        // SUBWORKFLOW: map_illumina
        // 
        MAP_ILLUMINA (
            ch_reads_bytype.illumina,
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai, index ] }
        )
        ch_versions = ch_versions.mix(MAP_ILLUMINA.out.versions)

        // Combine bam files with reference
        MAP_LONGREADS.out.bam_bai
            .mix(MAP_HIC.out.bam_bai)
            .mix(MAP_ILLUMINA.out.bam_bai)
            .set { ch_bam_bai }
    }
    
    if (steps.contains('stats')) {
    
        //
        // SUBWORKFLOW: bam_stats
        //    
        BAM_STATS (
            ch_bam_bai,
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai ] },
            []
        )
        ch_versions = ch_versions.mix(BAM_STATS.out.versions)

        //
        // SUBWORKFLOW: Generate depth per site reports
        //
        BAM_DEPTH (
            ch_bam_bai
        )
        ch_versions = ch_versions.mix(BAM_DEPTH.out.versions)
    }

    if (steps.contains('variant_calling')) {

        //
        // SUBWORKFLOW: variant_calling
        //
        VARIANT_CALLING (
            ch_bam_bai.filter { meta, bam, bai -> meta.type == 'hifi' },
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai ] },
            [],
            model_file,
            config_file
        )
        ch_versions = ch_versions.mix(VARIANT_CALLING.out.versions)

        //
        // SUBWORKFLOW: variant_calling_gatk
        //
        VARIANT_CALLING_GATK (
            ch_bam_bai.filter { meta, bam, bai -> meta.type == 'illumina' },
            PREPARE_GENOMES.out.genomes
                .map { meta, fasta, fai, dict, gtf, index -> [ meta, fasta, fai, dict ] },
            params.min_contig_length
        )
        ch_versions = ch_versions.mix(VARIANT_CALLING_GATK.out.versions)
    }

    if (steps.contains('report')) {

        // Generate window-wise stats
        BAM_STATS.out.depth
            .filter { meta, summary -> meta.type == 'hifi' }
            .map { meta, summary -> [ groupKey(meta.ref, meta.samples_per_type), [ meta.sample, summary ] ] }
            .groupTuple(sort: { a, b -> a[0] <=> b[0] })
            .map { ref, tuples -> [ ref.target, tuples.collect { it[0] }, tuples.collect { it[1] } ] }
            .set { ch_summaries }

        // Join depth and vcf files per reference
        BAM_DEPTH.out.depth
            .filter { meta, depth, tbi -> meta.type == 'hifi' }
            .map { meta, depth, tbi -> [ meta.ref, meta, depth, tbi ] }
            .join(VARIANT_CALLING.out.ind_vcf_tbi
                    .filter { meta, vcf, tbi -> meta.type == 'hifi' && meta.caller == 'clair3' }
                    .map { meta, vcf, tbi -> [ meta.ref, vcf, tbi ] },
                failOnDuplicate: true,
                failOnMismatch: true
            )
            .join(ch_summaries, failOnDuplicate: true, failOnMismatch: true)
            .join(PREPARE_GENOMES.out.genomes
                    .map { meta, fasta, fai, dict, gtf, index -> [ meta.id, meta, fasta, fai ] },
                failOnDuplicate: true,
                failOnMismatch: true
            )
            .multiMap { ref, meta, depth, dtbi, vcf, tbi, samples, summaries, meta2, fasta, fai ->
                input:     [ [id: 'winstats', type: meta.type, ref: ref, samples: tuple(samples)], depth, dtbi, vcf, tbi, summaries ]
                fasta_fai: [ meta2, fasta, fai ]
            }
            .set { ch_to_winstats }

        //
        // MODULE: windows_stats
        //
        WINDOWS_STATS (
            ch_to_winstats.input,
            ch_to_winstats.fasta_fai
        )
        ch_versions = ch_versions.mix(WINDOWS_STATS.out.versions.first())

        PLOT_WINDOWS (
            WINDOWS_STATS.out.bed
        )
        ch_versions = ch_versions.mix(PLOT_WINDOWS.out.versions.first())

/*
        // Stats per sample
        BAM_STATS.out.depth
            .map { meta, depth -> [ meta.sample, [ meta.type, depth ] ] }
            .groupTuple()
            .map { sample, tuples ->
                [ sample, tuples.find(it[0] == 'hifi')[1], tuples.find(it[0] == 'ul'), tuples.find(it[0] == 'hic')[1] ]
            }
            .set { ch_depth_stats }

        RUN_KMER_FK.out.summary
            .join(RUN_KMER_FK.out.stats, failOnDuplicate: true, failOnMismatch: true)
            .join(RUN_KMER_FK.out.qv, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, summary, stats, qv -> [ meta.sample, meta, summay, stats, qv ] }
            .set { ch_kmer_stats }

        ch_assemblies
            .map { meta, fasta, gtf, yaml -> [ meta.sample, meta, yaml ] }
            .join(ch_kmer_stats, failOnDuplicate: true, failOnMismatch: true)
            .join(BUSCO_BUSCO.out.short_summaries_txt
                .map { meta, busco -> [ meta.sample, busco ] },
                failOnDuplicate: true,
                failOnMismatch: true)
            .join(ch_depth_stats, failOnDuplicate: true, failOnMismatch: true)
            .multiMap { sample, meta, yaml, summary, stats, qv, busco, hifi, ul, hic ->
                yaml:    [ meta, yaml ]
                summary: [ meta, summary, [] ]
                depth:   [ meta, hifi, ul, hic ]
                merqury: [ meta, stats, qv ]
                busco:   [ meta, busco ]
            }
            .set { ch_to_ear }

        //
        // MODULE: generate_ear
        //
        GENERATE_EAR (
            ch_to_ear.yaml,
            ch_to_ear.summary,
            ch_to_ear.depth,
            ch_to_ear.merqury,
            ch_to_ear.busco
        )
        ch_versions = ch_versions.mix(GENERATE_EAR.out.versions.first())
*/

    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'asseval_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
