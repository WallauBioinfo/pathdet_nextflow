process FETCH_REFERENCE_GENOMES {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.out}/${sample_id}/${sample_id}_p04_ref", mode: 'copy'

    input:
    tuple val(sample_id), path(lineage_map)
    path ncbi_lineage_tsv
    path ncbi_genome_dir

    output:
    tuple val(sample_id), path("${sample_id}_ref.fasta"), emit: ref_fasta
    path "${sample_id}_ref_list.txt", optional: true    , emit: ref_list

    script:
    """
    if [ -s ${lineage_map} ] && [ \$(wc -l < ${lineage_map}) -gt 1 ]; then
        python3 ${projectDir}/bin/fetch_genomes.py \\
            ${lineage_map} \\
            ${ncbi_lineage_tsv} \\
            ${ncbi_genome_dir} \\
            .

        if [ -f REF.fasta ]; then
            mv REF.fasta ${sample_id}_ref.fasta
        else
            touch ${sample_id}_ref.fasta
        fi

        if [ -f ref_list.txt ]; then
            mv ref_list.txt ${sample_id}_ref_list.txt
        fi
    else
        touch ${sample_id}_ref.fasta
    fi
    """
}

process ALIGN_AND_VARIANT_CALLING {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.out}/${sample_id}/${sample_id}_p05_map", mode: 'copy'

    input:
    tuple val(sample_id), path(ref_fasta), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.bam")    , emit: bam
    tuple val(sample_id), path("${sample_id}.bam.bai"), emit: bai
    path "*.bw"                                       , optional: true
    path "*.vcf.gz"                                   , optional: true
    path "*.svg"                                      , optional: true

    script:
    """
    if [ -s ${ref_fasta} ]; then
        # 1. Alinhamento com minimap2 (igual ao script Bash original)
        minimap2 -ax sr -t ${task.cpus} ${ref_fasta} ${reads[0]} ${reads[1]} | \\
            samtools sort -@ ${task.cpus} -o ${sample_id}.bam
        samtools index ${sample_id}.bam

        # 2. Chamada de Variantes (bcftools)
        bcftools mpileup -f ${ref_fasta} ${sample_id}.bam | \\
            bcftools call -mv -Oz -o ${sample_id}.vcf.gz
        bcftools index ${sample_id}.vcf.gz

        # 3. Cobertura BigWig
        bamCoverage -b ${sample_id}.bam -o ${sample_id}.bw -p ${task.cpus}

        # 4. Gráfico SVG (spark.py): verifica no caminho original do Conda/git ou no PATH
        if [ -f "\$CONDA_PREFIX/opt/git/SparK/spark.py" ]; then
            python3 "\$CONDA_PREFIX/opt/git/SparK/spark.py" ${sample_id}.bam ${sample_id}
        else
            echo "warning: Spark.py not in \$CONDA_PREFIX/opt/git/SparK/"
        fi
    else
        touch ${sample_id}.bam ${sample_id}.bam.bai
    fi
    """
}
