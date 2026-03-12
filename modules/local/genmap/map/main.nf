process GENMAP_MAP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/44/449e9218f70b57e5b473c560ba5f0bf2e8e4a9f66b28d10f75868135868ba46d/data':
        'community.wave.seqera.io/library/genmap_ucsc-bedgraphtobigwig:34c62fd0055b0259' }"

    input:
    tuple val(meta), path(index), path(chromsizes)
    tuple val(meta2), path(regions)

    output:
    tuple val(meta), path("*.bw"), emit: bigwig
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed    = regions ? "--selection ${regions}" : ""

    if ("$index" == "${prefix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    genmap \\
        map \\
        ${args} \\
        ${bed} \\
        --bedgraph \\
        --threads ${task.cpus} \\
        --index ${index} \\
        --output ${prefix}

    bedGraphToBigWig \\
        ${prefix}.bedgraph \\
        $chromsizes \\
        ${prefix}.bw \\
        ${args2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version | sed 's/GenMap version: //; s/SeqAn.*\$//')
        bedGraphToBigWig: \$( echo \$(bedGraphToBigWig 2>&1) | head -n 1 | sed 's/^.*bedGraphToBigWig v //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ("$index" == "${prefix}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    touch ${prefix}.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version | sed 's/GenMap version: //; s/SeqAn.*\$//')
        bedGraphToBigWig: \$( echo \$(bedGraphToBigWig 2>&1) | head -n 1 | sed 's/^.*bedGraphToBigWig v //; s/ .*\$//')
    END_VERSIONS
    """
}
