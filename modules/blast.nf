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

    split_fasta.py \
        ${unclassified_fasta} \
        ${params.split_size} \
        chunks/${sample_id}_

    if [ -f chunks/${sample_id}_summary.txt ]; then
        mv chunks/${sample_id}_summary.txt ./${sample_id}_summary.txt
    else
        touch ./${sample_id}_summary.txt
    fi

    if ! ls chunks/${sample_id}_*.fasta 1> /dev/null 2>&1; then
        touch chunks/${sample_id}_01.fasta
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
    blastn \
        -db ${params.blast_nt} \
        -query ${fasta_chunk} \
        -out ${fasta_chunk.baseName}_blast.txt \
        -evalue 1.0e-10 \
        -max_target_seqs 1 \
        -max_hsps 5 \
        -outfmt "7 std" \
        -num_threads ${task.cpus}
    """
}

process MERGE_AND_KRONA {
    tag "${sample_id}"
    label 'process_medium'

    // p03_blast:
    publishDir "${params.out}/${sample_id}/${sample_id}_blast", mode: 'copy',
        saveAs: { fn ->
            fn.startsWith('blast_stage/') ? fn.tokenize('/').last() :
            (fn.endsWith('_blast_fmt7.txt') || fn.endsWith('_blast_hit.txt')) ? fn : null
        }
    // p04_krona:
    publishDir "${params.out}/${sample_id}/${sample_id}_tbl", mode: 'copy',
        saveAs: { fn ->
            (fn.startsWith('blast_stage/') || fn.startsWith("${sample_id}_CSV") ||
             fn.endsWith('_blast_fmt7.txt') || fn.endsWith('_blast_hit.txt')) ? null : fn
        }

    input:
    tuple val(sample_id), path(blast_result_files), path(log_p1)

    output:
    tuple val(sample_id), path("${sample_id}_blast_fmt7.txt")                  , emit: blast_fmt7
    tuple val(sample_id), path("${sample_id}_blast_hit.txt")                   , emit: blast_hit
    path "blast_stage/${sample_id}_pathogen_final.html"                        , emit: blast_stage_html
    tuple val(sample_id), path("${sample_id}_pathogen_final.html")             , emit: krona_html
    tuple val(sample_id), path("${sample_id}_CSV/${sample_id}_merge_taxid.csv"), emit: merge_taxid
    path "${sample_id}_fRep"                                                   , emit: frep_dir
    path "${sample_id}_bRep"                                                   , emit: brep_dir
    path "${sample_id}_report.csv"                                             , emit: report
    path "${sample_id}_bacteria_report.csv"                                    , emit: bacteria_report

    script:
    def taxonkit_db_dir = "${params.datadir}/taxonkit/taxdump"
    """
    set -e

    export TAXONKIT_DB="${taxonkit_db_dir}"
    export LC_ALL=C

    mkdir -p ${sample_id}_BLAST blast_stage ${sample_id}_CSV

    # ---- p03_blast:
    cat *_blast.txt > ${sample_id}_blast_fmt7.txt

    cat ${sample_id}_blast_fmt7.txt \
        | awk '/hits found/{getline;print}' \
        | grep -v "#" > ${sample_id}_blast_hit.txt || true

    ktClassifyBLAST -p -s ${sample_id}_blast_hit.txt -o ${sample_id}_BLAST/${sample_id}_blast_hit.csv
    cut -f 1,2 ${sample_id}_BLAST/${sample_id}_blast_hit.csv > ${sample_id}_BLAST/${sample_id}_blast_hit2.csv

    ktImportTaxonomy \
        -m 1 -s 3 \
        ${sample_id}_BLAST/${sample_id}_blast_hit2.csv \
        -o blast_stage/${sample_id}_pathogen_final.html

    # ---- p04_krona
    ktClassifyBLAST -p -s ${sample_id}_blast_hit.txt -o ${sample_id}_CSV/${sample_id}_blast.csv
    cut -f 1,2 ${sample_id}_CSV/${sample_id}_blast.csv > ${sample_id}_CSV/${sample_id}_merge.csv

    ktImportTaxonomy \
        -m 1 -s 3 \
        ${sample_id}_CSV/${sample_id}_merge.csv \
        -o ${sample_id}_pathogen_final.html

    cat ${sample_id}_CSV/${sample_id}_merge.csv \
        | sed -e '1d' \
        | awk '{ tax[\$2] += \$1; } END { for (i in tax) print tax[i] "\\t" i; }' \
        | sort -k 1nr,1 \
        > ${sample_id}_CSV/${sample_id}_merge_taxid.csv

    TR=\$(cat ${log_p1} | awk 'NR==3 {print \$4}' | tr -d ",")

    if [ -z "\${TR}" ] || [ "\${TR}" -eq 0 ] 2>/dev/null; then
        echo "[WARNING] Could not read total reads count from log. Using TR=1."
        TR=1
    fi

    taxonkit lineage -i 2 ${sample_id}_CSV/${sample_id}_merge.csv > ${sample_id}_CSV/${sample_id}_merge_tax.csv

    cat ${sample_id}_CSV/${sample_id}_merge_tax.csv \
        | awk -v "r=\${TR}" '{print \$1/r*1000000"\\t"\$0}' \
        | cut -f 1,3- \
        > ${sample_id}_CSV/${sample_id}_merge_tax2.csv

    # fRep/, bRep/, report.csv e bacteria_report.csv (Python, mesmos bytes do bash)
    process_blast_report.py \
        --sample-id ${sample_id} \
        --merge-tax ${sample_id}_CSV/${sample_id}_merge_tax.csv \
        --merge-tax2 ${sample_id}_CSV/${sample_id}_merge_tax2.csv \
        --total-reads "\${TR}" \
        --csv2-dir ${sample_id}_CSV2 \
        --frep-dir ${sample_id}_fRep \
        --brep-dir ${sample_id}_bRep \
        --report-path ${sample_id}_report.csv \
        --bacteria-report-path ${sample_id}_bacteria_report.csv \
        --threads ${task.cpus}
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
