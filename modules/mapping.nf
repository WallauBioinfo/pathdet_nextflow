process FETCH_REFERENCE_GENOMES {
    tag "$sample_id"
    label 'process_low'

    input:
    tuple val(sample_id), path(lineage_map)
    path ncbi_lineage_tsv
    path ncbi_genome_dir

    output:
    tuple val(sample_id), path("${sample_id}_ref.fasta"), emit: ref_fasta

    script:
    """
    python3 ${projectDir}/bin/fetch_references.py \
        --lineage-map ${lineage_map} \
        --ncbi-lineage ${ncbi_lineage_tsv} \
        --genome-dir ${ncbi_genome_dir} \
        --output ${sample_id}_ref.fasta
    """
}

process ALIGN_AND_VARIANT_CALLING {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.out}/${sample_id}/p05_map", mode: 'copy'

    input:
    tuple val(sample_id), path(ref_fasta), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam

    script:
    """
    bwa index ${ref_fasta}
    bwa mem -t ${task.cpus} ${ref_fasta} ${reads[0]} ${reads[1]} | samtools sort -@ ${task.cpus} -o ${sample_id}.bam
    samtools index ${sample_id}.bam
    """
}