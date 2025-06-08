process BWA_MAP_FILTER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0' :
        'biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:1bd8542a8a0b42e0981337910954371d0230828e-0' }"

    input:
    tuple val(meta) , path(reads)
    tuple val(meta2), path(index)
    val(mapq_filter)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    # map forward reads in single-end mode:
    bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        ${reads[0]} \\
        | filter_five_end.pl \\
        | samtools view $args2 --threads $task.cpus -o ${prefix}.R1.bam -

    # map reverse reads in single-end mode:
    bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        ${reads[1]} \\
        | filter_five_end.pl \\
        | samtools view $args2 --threads $task.cpus -o ${prefix}.R2.bam -

    # combine filtered single-end bam files:
    two_read_bam_combiner_mod.pl ${prefix}.R1.bam ${prefix}.R2.bam samtools $mapq_filter \\
    | samtools view $args2 --threads $task.cpus -o ${prefix}.bam -

    rm ${prefix}.R1.bam ${prefix}.R2.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
