process DEHOST_KRAKEN2 {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(reads)
    val kraken2_db

    output:
    tuple val(sample_id), path("${sample_id}_hgs_*.fastq.gz"), emit: reads
    path "${sample_id}.kraken.txt"                           , emit: report

    script:
    """
    kraken2 \
        --db ${kraken2_db} \
        --paired ${reads[0]} ${reads[1]} \
        --unclassified-out ${sample_id}_hgs#.fastq.gz \
        --threads ${task.cpus} \
        > ${sample_id}.kraken.txt
    """
}

process DEHOST_BOWTIE2 {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(reads)
    val bowtie2_index

    output:
    tuple val(sample_id), path("${sample_id}_hgs2.*.fastq"), emit: reads

    script:
    """
    bowtie2 \
        -x ${bowtie2_index} \
        -1 ${reads[0]} -2 ${reads[1]} \
        --un-conc ${sample_id}_hgs2.fastq \
        -S ${sample_id}_hgs2.sam \
        -p ${task.cpus}

    rm -f ${sample_id}_hgs2.sam
    """
}

process DEHOST_HISAT2 {
    tag "$sample_id"
    label 'process_high'

    input:
    tuple val(sample_id), path(reads)
    val hisat2_index

    output:
    tuple val(sample_id), path("${sample_id}_hgs.*.fastq"), emit: reads

    script:
    """
    hisat2 \
        -x ${hisat2_index} \
        -1 ${reads[0]} -2 ${reads[1]} \
        --un-conc ${sample_id}_hgs.fastq \
        -S ${sample_id}_hgs3.sam \
        -p ${task.cpus}

    rm -f ${sample_id}_hgs3.sam
    """
}

process DEHOST_BLASTN {
    tag "$sample_id"
    label 'process_medium'
    publishDir "${params.out}/${sample_id}/${sample_id}_hgs", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(log_p1)
    val blast_host2_db

    output:
    tuple val(sample_id), path("${sample_id}_nohit.fasta"), emit: fasta
    tuple val(sample_id), path("${sample_id}_log_p2.txt") , emit: log
    path "${sample_id}_hgs.txt"                           , emit: blast_out
    path "${sample_id}_nohitID.txt"                       , emit: nohit_id    

    script:
    def read2_param = reads instanceof List && reads.size() > 1 ? "${reads[1]}" : ""
    """
    seqkit fq2fa ${reads[0]} -o ${sample_id}_hgs.1.fasta

    blastn \
        -db ${blast_host2_db} \
        -query ${sample_id}_hgs.1.fasta \
        -evalue 1.0e-10 \
        -max_target_seqs 1 -max_hsps 5 \
        -outfmt "7 std" -num_threads ${task.cpus} > ${sample_id}_hgs.txt

    cat ${sample_id}_hgs.txt \
        | grep -B 2 "# 0 hits found" \
        | grep "Query" \
        | sed -e 's/# Query: //g' \
        | tr -d '\r' \
        > ${sample_id}_nohitID.txt

    if [ -s ${sample_id}_nohitID.txt ]; then
        seqkit grep -n -w0 -f ${sample_id}_nohitID.txt -o ${sample_id}_nohit.fasta ${sample_id}_hgs.1.fasta
    else
        touch ${sample_id}_nohit.fasta
    fi

    START2=\$(date '+%y%m%d%H%M%S')
    echo "\${START2} ${sample_id} trimmed FQ stats" > ${sample_id}_log_p2.txt
    cat ${log_p1} | sed -n '7,9p' >> ${sample_id}_log_p2.txt
    echo "" >> ${sample_id}_log_p2.txt
    echo "\$(date '+%y%m%d%H%M%S') ${sample_id} subtracted FQ/FA stats" >> ${sample_id}_log_p2.txt
    seqkit stats ${reads[0]} ${read2_param} ${sample_id}_nohit.fasta >> ${sample_id}_log_p2.txt
    """
}
