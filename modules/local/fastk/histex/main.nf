process FASTK_HISTEX {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastk:1.1.0--h71df26d_1':
        'biocontainers/fastk:1.1.0--h71df26d_1' }"

    input:
    tuple val(meta), path(hist)

    output:
    tuple val(meta), path("*_hist.txt")         , emit: hist
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "FASTK_HISTEX module does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def FASTK_VERSION = '1.1.0--h71df26d_1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    Histex \\
        $args \\
        $hist \\
        > ${prefix}_hist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_hist.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
    END_VERSIONS
    """
}
