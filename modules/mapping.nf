process FETCH_REFERENCE_GENOMES {
    tag "$sample_id"
    label 'process_low'

    publishDir "${params.out}/${sample_id}/${sample_id}_map", mode: 'copy'

    input:
    tuple val(sample_id), path(merge_taxid_csv)
    path ncbi_lineage_tsv
    path ncbi_genome_dir

    output:
    tuple val(sample_id), path("REF.fasta")    , emit: ref_fasta
    tuple val(sample_id), path("ref_list.txt") , emit: ref_list

    script:
    """
    python3 ${projectDir}/bin/fetch_genomes.py \\
        ${merge_taxid_csv} \\
        ${ncbi_lineage_tsv} \\
        ${ncbi_genome_dir} \\
        .
    """
}

process ALIGN_AND_VARIANT_CALLING {
    tag "$sample_id"
    label 'process_high'

    publishDir "${params.out}/${sample_id}/${sample_id}_map", mode: 'copy'

    input:
    tuple val(sample_id), path(ref_fasta), path(reads)

    output:
    path "${sample_id}_map/${sample_id}.sort.bam"               , emit: bam
    path "${sample_id}_map/${sample_id}.sort.bam.bai"           , emit: bai
    path "${sample_id}_map/${sample_id}_map.csv"                , emit: map_csv, optional: true
    path "${sample_id}_con/${sample_id}.vcf.gz"                 , emit: vcf
    path "${sample_id}_con/${sample_id}.vcf.gz.csi"             , emit: vcf_csi
    path "${sample_id}_con/${sample_id}.norm.bcf"               , emit: norm_bcf, optional: true
    path "${sample_id}_con/${sample_id}.norm.flt-idls.vcf.gz"    , emit: flt_vcf, optional: true
    path "${sample_id}_con/${sample_id}.norm.flt-idls.vcf.gz.csi", emit: flt_vcf_csi, optional: true
    path "${sample_id}_con/${sample_id}_consensus.fa"           , emit: consensus
    path "${sample_id}_bw/${sample_id}.bw"                      , emit: bw, optional: true
    path "${sample_id}_bw/${sample_id}.bdg.gz"                  , emit: bdg, optional: true
    path "${sample_id}_bw/${sample_id}.bdg.gz.tbi"              , emit: bdg_tbi, optional: true
    path "${sample_id}_bw/spark_*.svg"                          , emit: spark_svg, optional: true
    path "${sample_id}_bw/${sample_id}_ref.fa"                  , emit: bw_ref, optional: true
    path "MAP/REF.fasta"                                        , emit: ref_fasta
    path "MAP/REF.fasta.fai"                                    , emit: ref_fai
    path "MAP/ref_list*.txt"                                    , emit: ref_lists, optional: true

    script:
    """
    #!/bin/bash
    set -e

    mkdir -p ${sample_id}_map ${sample_id}_con ${sample_id}_bw MAP

    if [ ! -s ${ref_fasta} ]; then
        touch ${sample_id}_map/${sample_id}.sort.bam ${sample_id}_map/${sample_id}.sort.bam.bai
        touch ${sample_id}_con/${sample_id}.vcf.gz ${sample_id}_con/${sample_id}.vcf.gz.csi ${sample_id}_con/${sample_id}_consensus.fa
        touch ${sample_id}_bw/${sample_id}.bw ${sample_id}_bw/${sample_id}.bdg.gz ${sample_id}_bw/${sample_id}.bdg.gz.tbi
        touch MAP/REF.fasta MAP/REF.fasta.fai MAP/ref_list1.txt MAP/ref_list2.txt MAP/ref_list3.txt
        echo "0" > ${sample_id}_map/${sample_id}_map.csv
        exit 0
    fi

    # Cópia e indexação da referência no diretório MAP/
    cp ${ref_fasta} MAP/REF.fasta
    samtools faidx MAP/REF.fasta
    bwa index MAP/REF.fasta

    # Criação dos ref_list1, ref_list2 e ref_list3 conforme o pipeline original
    cut -f1 MAP/REF.fasta.fai > MAP/ref_list1.txt
    cut -f1,2 MAP/REF.fasta.fai > MAP/ref_list2.txt
    cp MAP/ref_list2.txt MAP/ref_list3.txt

    # Alinhamento BWA MEM sem descarte de unmapped no view para manter tamanho idêntico
    bwa mem -t ${task.cpus} MAP/REF.fasta ${reads[0]} ${reads[1]} | \
        samtools view -@ ${task.cpus} -bS - | \
        samtools sort -@ ${task.cpus} -o ${sample_id}_map/${sample_id}.sort.bam -
    
    samtools index ${sample_id}_map/${sample_id}.sort.bam

    # Formato do map.csv (contagem total de reads alinhadas)
    samtools view -c -F 4 ${sample_id}_map/${sample_id}.sort.bam > ${sample_id}_map/${sample_id}_map.csv

    # Cobertura e BigWig
    cp MAP/REF.fasta ${sample_id}_bw/${sample_id}_ref.fa
    bedtools genomecov -ibam ${sample_id}_map/${sample_id}.sort.bam -bg > temp.bdg || true

    if [ -s temp.bdg ]; then
        LC_ALL=C sort -k1,1 -k2,2n temp.bdg > temp.sorted.bdg
        bgzip -c temp.sorted.bdg > ${sample_id}_bw/${sample_id}.bdg.gz
        tabix -p bed ${sample_id}_bw/${sample_id}.bdg.gz || true

        cut -f1,2 MAP/REF.fasta.fai > chrom.sizes
        bedGraphToBigWig temp.sorted.bdg chrom.sizes ${sample_id}_bw/${sample_id}.bw || cp temp.sorted.bdg ${sample_id}_bw/${sample_id}.bw
        rm -f temp.bdg temp.sorted.bdg chrom.sizes
    else
        touch ${sample_id}_bw/${sample_id}.bw ${sample_id}_bw/${sample_id}.bdg.gz ${sample_id}_bw/${sample_id}.bdg.gz.tbi
    fi

    # Chamada de Variantes e Normalização BCFtools
    bcftools mpileup -f MAP/REF.fasta ${sample_id}_map/${sample_id}.sort.bam | \
        bcftools call -mv -Ob -o ${sample_id}_con/${sample_id}.bcf || true

    if [ -s ${sample_id}_con/${sample_id}.bcf ]; then
        bcftools norm -f MAP/REF.fasta ${sample_id}_con/${sample_id}.bcf -Ob -o ${sample_id}_con/${sample_id}.norm.bcf || true
        bcftools filter -s LOWQUAL -e 'QUAL<20 || DP<10' ${sample_id}_con/${sample_id}.norm.bcf -Oz -o ${sample_id}_con/${sample_id}.norm.flt-idls.vcf.gz || true
        
        if [ -s ${sample_id}_con/${sample_id}.norm.flt-idls.vcf.gz ]; then
            bcftools index ${sample_id}_con/${sample_id}.norm.flt-idls.vcf.gz || true
        fi

        bcftools view ${sample_id}_con/${sample_id}.bcf -Oz -o ${sample_id}_con/${sample_id}.vcf.gz || true
        bcftools index ${sample_id}_con/${sample_id}.vcf.gz || true
        bcftools consensus -f MAP/REF.fasta ${sample_id}_con/${sample_id}.vcf.gz > ${sample_id}_con/${sample_id}_consensus.fa || touch ${sample_id}_con/${sample_id}_consensus.fa
    else
        touch ${sample_id}_con/${sample_id}.vcf.gz ${sample_id}_con/${sample_id}.vcf.gz.csi ${sample_id}_con/${sample_id}_consensus.fa
    fi
    """
}
