process GENMAP_INDEX {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/44/449e9218f70b57e5b473c560ba5f0bf2e8e4a9f66b28d10f75868135868ba46d/data':
        'community.wave.seqera.io/library/genmap_ucsc-bedgraphtobigwig:34c62fd0055b0259' }"

    input:
    tuple val(meta), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}") , emit: index
    tuple val(meta), path("*.sizes")   , emit: sizes
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "$meta.id"

    """
    genmap \\
        index \\
        --fasta-file ${fasta} \\
        --index ${prefix} \\
        ${args}

    awk -v OFS='\\t' '{print \$1,\$2}' $fai > ${prefix}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version | sed 's/GenMap version: //; s/SeqAn.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "$meta.id"

    """
    touch ${prefix}
    touch ${prefix}.sizes

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genmap: \$(genmap --version | sed 's/GenMap version: //; s/SeqAn.*\$//')
    END_VERSIONS
    """
}
