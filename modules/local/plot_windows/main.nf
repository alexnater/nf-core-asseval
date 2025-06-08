process PLOT_WINDOWS {
    tag "$meta.id"
    label 'process_single'

    container "/data/users/anater/singularity_cache/r_container.sif"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.winstats.pdf"), emit: pdf
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    plot_windowstats.r \\
        --bed $bed \\
        --outfile ${prefix}.winstats.pdf \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Rscript: \$(Rscript --version | sed 's/^.*version \\(.*\\) (.*/\\1/')
    END_VERSIONS
    """
}