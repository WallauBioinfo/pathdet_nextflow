process FETCH_REFERENCE_GENOMES {
    tag "$sample_id"

    publishDir "${params.out}/${sample_id}/${sample_id}_map", mode: 'copy'

    input:
    tuple val(sample_id), path(merge_taxid_csv)
    path ncbi_genome_dir

    output:
    tuple val(sample_id), path("MAP/REF.fasta"), optional: true, emit: ref_fasta
    path "MAP/ref_list*.txt"                   , optional: true, emit: ref_lists_map
    path "ref_list*.txt"                       , optional: true, emit: ref_lists
    path "${sample_id}_log_p5.txt"                             , emit: log

    script:
    """
    python3 ${projectDir}/bin/fetch_genomes.py \
        ${merge_taxid_csv} \
        ${ncbi_genome_dir} \
        . \
        ${sample_id}
    """
}

process ALIGN_AND_VARIANT_CALLING {
    tag "$sample_id"

    publishDir "${params.out}/${sample_id}/${sample_id}_map", mode: 'copy'

    input:
    tuple val(sample_id), path(ref_fasta), path(reads)

    output:
    path "${sample_id}_map"  , emit: map_dir
    path "${sample_id}_con"  , emit: con_dir
    path "${sample_id}_bw"   , emit: bw_dir
    path "MAP/REF.fasta.fai" , emit: ref_fai

    script:
    def sort_mem = "${Math.max(256, (task.memory.toMega() / (task.cpus * 2)) as int)}M"
    def heap     = "${Math.max(1, (task.memory.toGiga() * 0.8) as int)}g"
    """
    mkdir -p MAP ${sample_id}_map ${sample_id}_con ${sample_id}_bw

    samtools faidx ${ref_fasta}
    cp ${ref_fasta}.fai MAP/REF.fasta.fai

    # Mapping (p05_map.sh, MODE=short)
    minimap2 -t ${task.cpus} -ax sr ${ref_fasta} ${reads[0]} ${reads[1]} \
        | samtools sort -@ ${task.cpus} -m ${sort_mem} -o ${sample_id}_map/${sample_id}.sort.bam
    samtools index -@ ${task.cpus} ${sample_id}_map/${sample_id}.sort.bam
    pileup.sh -Xmx${heap} in=${sample_id}_map/${sample_id}.sort.bam ref=${ref_fasta} out=${sample_id}_map/${sample_id}_map.csv

    # Variant calling and consensus
    bcftools mpileup -Ou -f ${ref_fasta} ${sample_id}_map/${sample_id}.sort.bam \
        | bcftools call --threads ${task.cpus} -mv -Oz -o ${sample_id}_con/${sample_id}.vcf.gz
    bcftools index --threads ${task.cpus} ${sample_id}_con/${sample_id}.vcf.gz

    bcftools norm --threads ${task.cpus} -f ${ref_fasta} ${sample_id}_con/${sample_id}.vcf.gz \
        -Ob -o ${sample_id}_con/${sample_id}.norm.bcf
    bcftools filter --threads ${task.cpus} --IndelGap 5 ${sample_id}_con/${sample_id}.norm.bcf \
        -Oz -o ${sample_id}_con/${sample_id}.norm.flt-idls.vcf.gz
    bcftools index --threads ${task.cpus} ${sample_id}_con/${sample_id}.norm.flt-idls.vcf.gz

    cat ${ref_fasta} \
        | bcftools consensus ${sample_id}_con/${sample_id}.norm.flt-idls.vcf.gz \
        > ${sample_id}_con/${sample_id}_consensus.fa

    # BigWig and bedGraph
    bamCoverage -b ${sample_id}_map/${sample_id}.sort.bam -o ${sample_id}_bw/${sample_id}.bw -of bigwig \
        --binSize=10 --minMappingQuality 10 --numberOfProcessors=${task.cpus}
    bamCoverage -b ${sample_id}_map/${sample_id}.sort.bam -o ${sample_id}_bw/${sample_id}.bdg -bs 1 -of bedgraph --numberOfProcessors=${task.cpus}
    bgzip -@ ${task.cpus} ${sample_id}_bw/${sample_id}.bdg
    tabix -p bed ${sample_id}_bw/${sample_id}.bdg.gz

    # spark plots (SVG), one by contig
    cp ${ref_fasta} ${sample_id}_bw/${sample_id}_ref.fa
    while read -r contig len _ ; do
        python ${params.spark_script} -pr "\$contig:1-\$len" \
            -cf ${sample_id}_bw/${sample_id}.bdg.gz -o ${sample_id}_bw/spark_\$contig
    done < ${ref_fasta}.fai
    """
}
