

process QUALITY_CONTROL {
    
    tag "$sample_id"
    publishDir "${params.out}/${sample_id}/${sample_id}_qc", mode: 'copy'

    input:
    
    tuple val(sample_id), path(reads)
    
    output:
    
    tuple val(sample_id), path("${sample_id}_trimmed2_R*.fastq.gz"), emit: sample
    tuple val(sample_id), path("${sample_id}_log_p1.txt"),           emit: log
    tuple val(sample_id), path('*.json'),                            emit: fjson
    tuple val(sample_id), path('*_fp.html'),                         emit: fhtml

    script:

    def bbduk_threads = task.cpus < 2 ? 2 : task.cpus

    """
    START1=\$(date '+%y%m%d%H%M%S')
    echo "\${START1} ${sample_id} raw FQ stats" > ${sample_id}_log_p1.txt
    seqkit stats -j ${task.cpus} ${reads[0]} ${reads[1]} >> ${sample_id}_log_p1.txt

    fastp -i ${reads[0]} \\
          -I ${reads[1]} \\
          -o ${sample_id}_trimmed_R1.fastq.gz \\
          -O ${sample_id}_trimmed_R2.fastq.gz \\
          -h ${sample_id}_fp.html \\
          -j ${sample_id}.json \\
          -3 \\
          -q ${params.fastp_qualified_phred} \\
          -n ${params.fastp_n_base_limit} \\
          -l ${params.fastp_length_required} \\
          -f ${params.fastp_trim_front1} \\
          -t ${params.fastp_trim_tail1} \\
          -T ${params.fastp_trim_tail2} \\
          -w ${task.cpus}

    bbduk.sh \\
	  -Xmx8g \\
          in=${sample_id}_trimmed_R1.fastq.gz \\
          in2=${sample_id}_trimmed_R2.fastq.gz \\
          out=${sample_id}_trimmed2_R1.fastq.gz \\
          out2=${sample_id}_trimmed2_R2.fastq.gz \\
          entropy=${params.bbduk_entropy} \\
          entropywindow=${params.bbduk_entropywindow} \\
          entropyk=${params.bbduk_entropyk} \\
          threads=${bbduk_threads} \\
          unpigz=t
    
    QC=\$(date '+%y%m%d%H%M%S')
    echo -e "\\n\${QC} ${sample_id} trimmed FQ stats" >> ${sample_id}_log_p1.txt
    seqkit stats -j ${task.cpus} ${sample_id}_trimmed2_R1.fastq.gz ${sample_id}_trimmed2_R2.fastq.gz --quiet >> ${sample_id}_log_p1.txt
    """
}
