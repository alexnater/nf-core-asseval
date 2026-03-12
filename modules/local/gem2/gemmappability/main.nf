process GEM2_GEMMAPPABILITY {
    tag "$meta.id"
    label 'process_high'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/23/239a50bbee5154eb584f8d7c46d9c4f6d5a563ac8b3f9b3f883fe308cbc7ef58/data':
        'community.wave.seqera.io/library/gem2_ucsc-bedgraphtobigwig:d215a32896887bf0' }"

    input:
    tuple val(meta), path(index)
    val(read_length)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: ''
    def args2   = task.ext.args2  ?: ''
    def args3   = task.ext.args3  ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20200110' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    gem-mappability \\
        -I ${index} \\
        -o ${prefix} \\
        -l ${read_length} \\
        -T ${task.cpus} \\
        ${args}

    gem-2-bed mappability \\
        --input ${prefix}.mappability \\
        --index ${index} \\
        --output ${prefix} \\
        ${args2}

    bedGraphToBigWig \\
        ${prefix}.bg \\
        ${prefix}.sizes \\
        ${prefix}.bw \\
        ${args3}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem2: $VERSION
        bedGraphToBigWig: \$( echo \$(bedGraphToBigWig 2>&1) | head -n 1 | sed 's/^.*bedGraphToBigWig v //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '20200110' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gem2: $VERSION
        bedGraphToBigWig: \$( echo \$(bedGraphToBigWig 2>&1) | head -n 1 | sed 's/^.*bedGraphToBigWig v //; s/ .*\$//')
    END_VERSIONS
    """
}
