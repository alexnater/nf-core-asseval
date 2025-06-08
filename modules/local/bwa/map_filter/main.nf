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
    if (task.cpus < 4) {
        error "Number of CPUs need to be at least 4"
    }
    
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def bwa_cpus = ((task.cpus - 2) / 2) as int

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    mkfifo R1_fifo R2_fifo

    # map forward reads in single-end mode:
    ( bwa mem \\
        $args \\
        -t $bwa_cpus \\
        \$INDEX \\
        ${reads[0]} \\
        | filter_five_end.pl > R1_fifo ) &

    # map reverse reads in single-end mode:
    ( bwa mem \\
        $args \\
        -t $bwa_cpus \\
        \$INDEX \\
        ${reads[1]} \\
        | filter_five_end.pl > R2_fifo ) &

    # combine filtered single-end bam files:
    ( two_read_bam_combiner_fifo.pl R1_fifo R2_fifo $mapq_filter \\
    | samtools view $args2 --threads 1 -o ${prefix}.bam - ) &

    wait

    rm R1_fifo R2_fifo

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
