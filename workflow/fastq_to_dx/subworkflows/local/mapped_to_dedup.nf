process GATK_MARKDUPLICATES {

    tag "${meta.id}"

    cpus 4
    memory '12 GB'

    publishDir "${params.outdir ?: 'results'}/dedup_mapped", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}_dedup.bam"), emit: bam
    tuple val(meta), path("${meta.id}_dedup.bam.bai"), emit: bai
    tuple val(meta), path("${meta.id}_dup_metrics.txt"), emit: metrics

    script:
    """
    gatk --java-options "-Xmx10g" MarkDuplicates \
        -I ${bam} \
        -O ${meta.id}_dedup.bam \
        -M ${meta.id}_dup_metrics.txt \
        --REMOVE_DUPLICATES true

    samtools index -@ ${task.cpus} ${meta.id}_dedup.bam
    """
}


workflow MAPPED_TO_DEDUP {

    take:
    ch_mapped_bam

    main:

    ch_input = ch_mapped_bam.map { meta, bam ->

        def dedup_meta = meta + [
            id: meta.original_id ?: meta.id,
            qc_stage: 'dedup_mapped'
        ]

        tuple(
            dedup_meta,
            bam
        )
    }

    GATK_MARKDUPLICATES(
        ch_input
    )

    emit:

    bam = GATK_MARKDUPLICATES.out.bam
    bai = GATK_MARKDUPLICATES.out.bai
    metrics = GATK_MARKDUPLICATES.out.metrics
}
