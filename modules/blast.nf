nextflow.enable.dsl = 2

process SPLIT_FASTA {
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.out}/${sample_id}/${sample_id}_blast", mode: 'copy', pattern: "*_summary.txt"

    input:
    tuple val(sample_id), path(unclassified_fasta)

    output:
    tuple val(sample_id), path("chunks/${sample_id}_*.fasta"), emit: chunks
    tuple val(sample_id), path("${sample_id}_summary.txt")   , emit: summary

    script:
    """
    mkdir -p chunks

    split_fasta.py \\
        ${unclassified_fasta} \\
        ${params.split_size} \\
        chunks/${sample_id}_

    # 2. Garante que exista o arquivo de summary
    if [ -f chunks/${sample_id}_summary.txt ]; then
        mv chunks/${sample_id}_summary.txt ./${sample_id}_summary.txt
    else
        touch ./${sample_id}_summary.txt
    fi

    # 3. Garante que exista pelo menos 1 chunk de FASTA para não esvaziar o canal
    if ! ls chunks/${sample_id}_*.fasta 1> /dev/null 2>&1; then
        touch chunks/${sample_id}_chunk_0.fasta
    fi
    """
}

process RUN_BLASTN {
    tag "${sample_id} - ${fasta_chunk.simpleName}"
    label 'process_high'
    cpus params.blast_threads

    input:
    tuple val(sample_id), path(fasta_chunk)

    output:
    tuple val(sample_id), path("${fasta_chunk.baseName}_blast.txt"), emit: blast_output

    script:
    """
    blastn \\
        -db ${params.blast_nt} \\
        -query ${fasta_chunk} \\
        -out ${fasta_chunk.baseName}_blast.txt \\
        -evalue 1.0e-10 \\
        -max_target_seqs 1 \\
        -max_hsps 5 \\
        -outfmt "7 std" \\
        -num_threads ${task.cpus}
    """
}

process MERGE_AND_KRONA {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.out}/${sample_id}/${sample_id}_blast", mode: 'copy'

    input:
    tuple val(sample_id), path(blast_result_files)

    output:
    tuple val(sample_id), path("${sample_id}_blast_fmt7.txt")     , emit: blast_fmt7
    tuple val(sample_id), path("${sample_id}_blast_hit.txt")      , emit: blast_hit
    tuple val(sample_id), path("${sample_id}_pathogen_final.html"), emit: krona_html

    script:
    """
    cat ${blast_result_files} > ${sample_id}_blast_fmt7.txt

    cat ${sample_id}_blast_fmt7.txt \\
        | awk '/hits found/{getline;print}' \\
        | grep -v "#" > ${sample_id}_blast_hit.txt || true

    mkdir -p krona_temporary
    ktClassifyBLAST -p -s ${sample_id}_blast_hit.txt -o krona_temporary/${sample_id}_blast_hit.csv
    cut -f 1,2 krona_temporary/${sample_id}_blast_hit.csv > krona_temporary/${sample_id}_blast_hit2.csv

    ktImportTaxonomy \\
        -m 1 -s 3 \\
        krona_temporary/${sample_id}_blast_hit2.csv \\
        -o ${sample_id}_pathogen_final.html
    """
}

process PARSE_TAXID_GENOMES {
    tag "${sample_id}"
    label 'process_low'

    publishDir "${params.out}/${sample_id}/${sample_id}_blast", mode: 'copy'

    input:
    tuple val(sample_id), path(blast_hit_file)
    path ncbi_lineage_tsv

    output:
    tuple val(sample_id), path("${sample_id}_taxid_lineage_genomes.tsv"), emit: lineage_map

    script:
    """
    cp ${ncbi_lineage_tsv} ${sample_id}_taxid_lineage_genomes.tsv
    """
}
