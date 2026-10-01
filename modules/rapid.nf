process RAPID_ANALYSIS {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.out}/${sample_id}/${sample_id}_rapid", mode: 'copy'

    input:
    tuple val(sample_id), path(fasta_nohit)
    tuple val(sample_id), path(log_p1)
    tuple val(sample_id), path(log_p2)
    val pathdb1_db

    output:
    tuple val(sample_id), path("${sample_id}_pre-report.csv"), emit: pre_report
    path "${sample_id}_pRep/${sample_id}.plot.html"           , emit: krona_html
    tuple val(sample_id), path("${sample_id}_pRep")           , emit: rank_csvs
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
   
   # Cria os arquivos vazios APENAS se o Kraken2 não os tiver criado:
    [ -f KK2/${sample_id}_kraken_p1.txt ] || touch KK2/${sample_id}_kraken_p1.txt
    [ -f KK2/${sample_id}_class2.fastq ]   || touch KK2/${sample_id}_class2.fastq
    [ -f KK2/${sample_id}_unclass2.fastq ] || touch KK2/${sample_id}_unclass2.fastq    

    # 2. Geração do gráfico Krona (só roda se o arquivo tiver conteúdo / tamanho > 0)
    if [ -s KK2/${sample_id}_kraken_p1.txt ]; then
        ktImportTaxonomy \
            -q 2 -t 3 -s 4 \
            -o ${sample_id}_pRep/${sample_id}.plot.html \
            KK2/${sample_id}_kraken_p1.txt
    else
        echo "<html><body><h3>Nenhum patógeno detectado (0 reads)</h3></body></html>" > ${sample_id}_pRep/${sample_id}.plot.html
    fi
    
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

process GENERATE_TBL_FOLDER {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.out}/${sample_id}/${sample_id}_tbl", mode: 'copy'

    input:
    tuple val(sample_id), path(pre_report), path(prep_dir)

    output:
    path "${sample_id}_report.csv"         , optional: true
    path "${sample_id}_bacteria_report.csv", optional: true
    path "${sample_id}_bRep/*"             , optional: true

    script:
    """
    mkdir -p ${sample_id}_bRep

    if [ -f ${pre_report} ]; then
        cp ${pre_report} ./${sample_id}_report.csv
        cp ${pre_report} ./${sample_id}_bacteria_report.csv
    fi

    if [ -d ${prep_dir} ]; then
        # Copia mantendo os originais e cria as versões padronizadas Bacteria_f*.csv
        cp ${prep_dir}/*.csv ./${sample_id}_bRep/ 2>/dev/null || true
        
        for file in ./${sample_id}_bRep/${sample_id}_p*.csv; do
            if [ -f "\$file" ]; then
                new_name=\$(echo "\$file" | sed "s/${sample_id}_p/${sample_id}_Bacteria_f/")
                cp "\$file" "\$new_name"
                rm -f "\$file"
            fi
        done
    fi
    """
}
