process AWK_SAMTOOLS_STATS_PARSE_MQ {
    tag "${meta.id}"
    label "process_low"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("${meta.id}_mapq_hist.tsv"), emit: mapq_hist

    script:
    """
    awk '\$1=="MAPQ"{print \$2"\\t"\$3}' ${stats} > ${meta.id}_mapq_hist.tsv
    """
}
