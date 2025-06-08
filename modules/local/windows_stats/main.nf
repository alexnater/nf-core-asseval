process WINDOWS_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.23.0--py39hdd5828d_0' :
        'quay.io/biocontainers/pysam:0.23.0--py39hdd5828d_0' }"

    input:
    tuple val(meta), path(depth), path(dtbi), path(vcf), path(tbi), path(mosdepth)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.summary.tsv"), emit: summary
    tuple val(meta), path("*.bed")        , emit: bed
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    windows_stats.py \\
        --depth $depth \\
        --vcf $vcf \\
        --outprefix $prefix \\
        --samples ${meta.samples.join(' ')} \\
        --summaries $mosdepth \\
        --fasta $fasta \\
        --fai $fai \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}