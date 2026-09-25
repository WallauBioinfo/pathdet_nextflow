process FETCH_REFERENCE_GENOMES {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.out}/${sample_id}/p04_references", mode: 'copy'

    input:
    tuple val(sample_id), path(lineage_map)
    path ncbi_lineage_tsv
    path ncbi_genome_dir

    output:
    tuple val(sample_id), path("${sample_id}_ref.fasta"), emit: ref_fasta

    script:
    """
    if [ -s ${lineage_map} ] && [ \$(wc -l < ${lineage_map}) -gt 1 ]; then
        
        python3 ${projectDir}/bin/fetch_genomes.py \\
            ${lineage_map} \\
            ${ncbi_lineage_tsv} \\
            ${ncbi_genome_dir} \\
            ncbi_genome

        if [ -f ncbi_genome/REF.fasta ]; then
            mv ncbi_genome/REF.fasta ${sample_id}_ref.fasta
        else
            touch ${sample_id}_ref.fasta
        fi

    else
        touch ${sample_id}_ref.fasta
    fi
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
    if [ -s ${ref_fasta} ]; then
        bwa index ${ref_fasta}
        bwa mem -t ${task.cpus} ${ref_fasta} ${reads[0]} ${reads[1]} | samtools sort -@ ${task.cpus} -o ${sample_id}.bam
        samtools index ${sample_id}.bam
    else
        touch ${sample_id}.bam
    fi
    """
}
