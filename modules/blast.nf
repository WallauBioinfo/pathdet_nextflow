
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

    publishDir "${params.out}/${sample_id}/${sample_id}_tbl", mode: 'copy'

    input:
    tuple val(sample_id), path(blast_result_files), path(log_p1)

    output:
    tuple val(sample_id), path("${sample_id}_blast_fmt7.txt")                , emit: blast_fmt7
    tuple val(sample_id), path("${sample_id}_blast_hit.txt")                 , emit: blast_hit
    tuple val(sample_id), path("${sample_id}_pathogen_final.html")           , emit: krona_html
    tuple val(sample_id), path("${sample_id}_CSV/${sample_id}_merge_taxid.csv"), emit: merge_taxid
    path "${sample_id}_CSV/*"                                                , emit: csv_dir
    path "${sample_id}_CSV2/*"                                               , emit: csv2_dir
    path "${sample_id}_fRep/*"                                               , emit: frep_dir

    script:
    def taxonkit_db_dir = "${params.datadir}/taxonkit/taxdump"
    """
    #!/bin/bash
    set -e

    export TAXONKIT_DB="${taxonkit_db_dir}"

    mkdir -p ${sample_id}_CSV ${sample_id}_CSV2 ${sample_id}_fRep ${sample_id}_bRep

    cat ${blast_result_files} > ${sample_id}_blast_fmt7.txt

    cat ${sample_id}_blast_fmt7.txt \\
        | awk '/hits found/{getline;print}' \\
        | grep -v "#" > ${sample_id}_blast_hit.txt || true

    ktClassifyBLAST -p -s ${sample_id}_blast_hit.txt -o ${sample_id}_CSV/${sample_id}_blast.csv
    cut -f 1,2 ${sample_id}_CSV/${sample_id}_blast.csv > ${sample_id}_CSV/${sample_id}_merge.csv

    ktImportTaxonomy \\
        -m 1 -s 3 \\
        ${sample_id}_CSV/${sample_id}_merge.csv \\
        -o ${sample_id}_pathogen_final.html

    cat ${sample_id}_CSV/${sample_id}_merge.csv \\
        | sed -e '1d' \\
        | awk '{ tax[\$2] += \$1; } END { for (i in tax) print tax[i] "\t" i; }' \\
        | sort -k 1nr,1 \\
        > ${sample_id}_CSV/${sample_id}_merge_taxid.csv

    # Extrai a contagem de reads de R1 (coluna 5, linha 3)
    TR=\$(awk 'NR==3 {gsub(/,/, "", \$5); print \$5}' ${log_p1})

    # Validação de segurança
    if [ -z "\${TR}" ] || [ "\${TR}" -eq 0 ] 2>/dev/null; then
        echo "[WARNING] Não foi possível ler o número total de reads do log. Usando TR=1."
        TR=1
    fi

    taxonkit lineage -i 2 ${sample_id}_CSV/${sample_id}_merge.csv > ${sample_id}_CSV/${sample_id}_merge_tax.csv

    cat ${sample_id}_CSV/${sample_id}_merge_tax.csv \\
        | awk -v "r=\${TR}" '{print \$1/r*1000000"\t"\$0}' \\
        | cut -f 1,3- \\
        > ${sample_id}_CSV/${sample_id}_merge_tax2.csv

    # Processamento de patentes taxonômicas K, P, F, G, S
    for rank in k p f g s ; do
        # Garante a existência do arquivo temporário
        touch tmp_\${rank}.txt

        # Processa mf (Normalized / RPM)
        taxonkit reformat -j ${task.cpus} -i 3 -f "{\$rank}" ${sample_id}_CSV/${sample_id}_merge_tax2.csv \\
            | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[\$2] += \$1; } END { for (i in tax) print i "\t" tax[i]; }' \\
            | sort -k 2nr,2 | grep -v -e '^\\s' > tmp_\${rank}.txt || true

        # Calcula a Abundância Relativa (RA) se o arquivo não estiver vazio
        if [ -s tmp_\${rank}.txt ]; then
            awk 'NR==FNR { sum += \$2; next } { p = (sum > 0 ? \$2/sum : 0); print \$0 "\t" p }' tmp_\${rank}.txt tmp_\${rank}.txt \\
                > ${sample_id}_CSV2/${sample_id}_mf\${rank^^}.csv
        else
            touch ${sample_id}_CSV2/${sample_id}_mf\${rank^^}.csv
        fi
        rm -f tmp_\${rank}.txt

        # Processa rf (Raw counts)
        taxonkit reformat -j ${task.cpus} -i 3 -f "{\$rank}" ${sample_id}_CSV/${sample_id}_merge_tax.csv \\
            | cut -f 1,4- | sed 's/ /_/g' | awk '{ tax[\$2] += \$1; } END { for (i in tax) print i "\t" tax[i]; }' \\
            | sort -k 2nr,2 | grep -v -e '^\\s' > ${sample_id}_CSV2/${sample_id}_rf\${rank^^}.csv || true

        # Unifica as saídas se existirem dados
        if [ -s ${sample_id}_CSV2/${sample_id}_rf\${rank^^}.csv ] && [ -s ${sample_id}_CSV2/${sample_id}_mf\${rank^^}.csv ]; then
            paste ${sample_id}_CSV2/${sample_id}_rf\${rank^^}.csv ${sample_id}_CSV2/${sample_id}_mf\${rank^^}.csv \\
                | cut -f 1,2,4,5 \\
                > ${sample_id}_fRep/${sample_id}_f\${rank^^}1.csv || true
        else
            touch ${sample_id}_fRep/${sample_id}_f\${rank^^}1.csv
        fi

        # Cabeçalho do relatório final
        echo "Microorganisms,reads,RPM,RA" > ${sample_id}_fRep/${sample_id}_f\${rank^^}.csv
        if [ -s ${sample_id}_fRep/${sample_id}_f\${rank^^}1.csv ]; then
            cat ${sample_id}_fRep/${sample_id}_f\${rank^^}1.csv >> ${sample_id}_fRep/${sample_id}_f\${rank^^}.csv
        fi
        rm -f ${sample_id}_fRep/${sample_id}_f\${rank^^}1.csv
    done
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
