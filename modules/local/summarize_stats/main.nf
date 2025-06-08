process SUMMARIZE_STATS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(flagstats), path(depths), path(coverages), val(before), val(after)

    output:
    tuple val(meta), path("*_report.tsv"), emit: summary
    path "versions.yml"                  , emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    summarize_stats.py \\
        -f $flagstats \\
        -s $depths \\
        -c $coverages \\
        -b ${before.join(' ')} \\
        -a ${after.join(' ')} \\
        $args \\
        -o ${prefix}_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}