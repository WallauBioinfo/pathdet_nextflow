process RAPID_ANALYSIS {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.out}/${sample_id}/rapid", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta_nohit)
    tuple val(sample_id), path(log_p1)
    tuple val(sample_id), path(log_p2)
    val pathdb1_db

    output:
    tuple val(sample_id), path("${sample_id}_pre-report.csv"), emit: pre_report
    path "${sample_id}_pRep/${sample_id}.plot.html"           , emit: krona_html
    path "${sample_id}_pRep/*.csv"                            , emit: rank_csvs
    path "${sample_id}_log_Rapid.txt"                         , emit: log

    script:
    """
    mkdir -p KK2 ${sample_id}_pRep

    # 1. Alinhamento e classificação rápida por k-meros
    kraken2 \
        --report KK2/${sample_id}_p1.kreport \
        --unclassified-out KK2/${sample_id}_unclass2.fastq \
        --classified-out KK2/${sample_id}_class2.fastq \
        --db ${pathdb1_db} \
        --confidence 0.5 \
        --output KK2/${sample_id}_kraken_p1.txt \
        --threads ${task.cpus} \
        ${fasta_nohit}

    # 2. Geração do gráfico interativo Krona
    ktImportTaxonomy \
        -q 2 -t 3 -s 4 \
        -o ${sample_id}_pRep/${sample_id}.plot.html \
        KK2/${sample_id}_kraken_p1.txt

    # 3. Orquestrador Python (Limpeza de dados, ete3, Pandas e Relatório Final)
    process_rapid_report.py \
        --sample-id ${sample_id} \
        --log-p1-path ${log_p1} \
        --log-p2-path ${log_p2} \
        --kraken-output-path KK2/${sample_id}_kraken_p1.txt \
        --classified-fastq-path KK2/${sample_id}_class2.fastq \
        --unclassified-fastq-path KK2/${sample_id}_unclass2.fastq \
        --output-directory ${sample_id}_pRep \
        --output-report-path ${sample_id}_pre-report.csv \
        --output-log-path ${sample_id}_log_Rapid.txt
    """
}