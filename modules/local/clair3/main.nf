process CLAIR3 {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.2.0--py310h779eee5_0' :
        'biocontainers/clair3:1.2.0--py310h779eee5_0' }"

    input:
    tuple val(meta) , path(bam)  , path(bai), path(bed)
    tuple val(meta2), path(fasta), path(fai)
    val(platform)
    path(model)

    output:
    tuple val(meta), path("${prefix}_merge_output.vcf.gz")        , emit: vcf
    tuple val(meta), path("${prefix}_phased_merge_output.gvcf.gz"), emit: phased, optional: true    
    tuple val(meta), path("${prefix}_merge_output.gvcf.gz")       , emit: gvcf  , optional: true
    tuple val(meta), path("${prefix}_pileup.vcf.gz")              , emit: pileup, optional: true
    tuple val(meta), path("${prefix}_full_alignment.vcf.gz")      , emit: full  , optional: true
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bed_arg = bed ? "--bed_fn=${bed}" : ''
    def model_path = model ? "${model}" : "/usr/local/bin/models/${platform}"

    """
    run_clair3.sh \\
        $args \\
        --threads=$task.cpus \\
        --sample_name=${meta.id} \\
        --bam_fn=$bam \\
        --ref_fn=$fasta \\
        --output="." \\
        --platform=$platform \\
        --model_path=$model_path \\
        $bed_arg

    for file in *.{g,}vcf.gz; do
        mv \$file ${prefix}_\${file}
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version | sed 's/Clair3 //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intermediate = args.contain("--remove_intermediate_dir") ? 0 : 1
    def phased = args.contain("--enable_phasing") ? 1 : 0
    def gvcf = args.contain("--gvcf") ? 1 : 0
    """
    if [ $intermediate -eq 1 ]; then
        touch ${prefix}_pileup.vcf.gz ${prefix}_full_alignment.vcf.gz
    fi
    touch ${prefix}_merge_output.vcf.gz
    if [ $phased -eq 1 ]; then
        touch ${prefix}_phased_merge_output.gvcf.gz
    fi
    if [ $gvcf -eq 1 ]; then
        touch ${prefix}_merge_output.gvcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version | sed 's/Clair3 //')
    END_VERSIONS
    """
}