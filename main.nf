nextflow.enable.dsl = 2

ANSI_GREEN = "\033[1;32m"
ANSI_RED   = "\033[1;31m"
ANSI_RESET = "\033[0m"

include { QUALITY_CONTROL } from './modules/quality_control.nf'
include { DEHOST_KRAKEN2  } from './modules/dehosting.nf'
include { DEHOST_BOWTIE2  } from './modules/dehosting.nf'
include { DEHOST_HISAT2   } from './modules/dehosting.nf'
include { DEHOST_BLASTN   } from './modules/dehosting.nf'
include { RAPID_ANALYSIS; GENERATE_TBL_FOLDER } from './modules/rapid.nf'
include { 
    SPLIT_FASTA; 
    RUN_BLASTN; 
    MERGE_AND_KRONA; 
    PARSE_TAXID_GENOMES 
} from './modules/blast.nf'
include {
    FETCH_REFERENCE_GENOMES;
    ALIGN_AND_VARIANT_CALLING
} from './modules/mapping.nf'

input_csv = Channel.fromPath(params.input_csv, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
        def resolvePath = { path ->
            path.startsWith('assets/') || path.startsWith('./assets/') ? "$launchDir/${path}" : path
        }
        def sample_id = row.sample_id
        def reads_r1 = file(resolvePath(row.reads_r1))
        def reads_r2 = file(resolvePath(row.reads_r2))
        return tuple(sample_id, [reads_r1, reads_r2])
    }

workflow {

    // 1. Controle de Qualidade (p1)
    QUALITY_CONTROL(input_csv)

    // 2. Dehosting (p2)
    DEHOST_KRAKEN2 ( QUALITY_CONTROL.out.sample, params.kraken2_grch38 )
    DEHOST_BOWTIE2 ( DEHOST_KRAKEN2.out.reads,   params.bowtie2_host )
    DEHOST_HISAT2  ( DEHOST_BOWTIE2.out.reads,   params.hisat2_grch38r )
    DEHOST_BLASTN  ( DEHOST_HISAT2.out.reads,    params.blast_host2, QUALITY_CONTROL.out.log )

    // 3. Análise Rápida por k-meros (p03_rapid)
    RAPID_ANALYSIS ( 
        DEHOST_BLASTN.out.fasta, 
        QUALITY_CONTROL.out.log, 
        DEHOST_BLASTN.out.log, 
        params.pathdb1_db 
    )

    GENERATE_TBL_FOLDER (
	RAPID_ANALYSIS.out.pre_report.join(RAPID_ANALYSIS.out.rank_csvs)
    )

    // 4. Módulo BLAST (p03_blast)
    SPLIT_FASTA ( DEHOST_BLASTN.out.fasta )
    ch_blast_chunks = SPLIT_FASTA.out.chunks.transpose()
    RUN_BLASTN ( ch_blast_chunks )
    ch_grouped_blast = RUN_BLASTN.out.blast_output.groupTuple()
    MERGE_AND_KRONA ( ch_grouped_blast )

    PARSE_TAXID_GENOMES (
        MERGE_AND_KRONA.out.blast_hit,
        file(params.ncbi_lineage_tsv)
    )

    // 5. Mapeamento, Cobertura e Chamada de Variantes (p05_map)
    FETCH_REFERENCE_GENOMES (
        PARSE_TAXID_GENOMES.out.lineage_map,
        file(params.ncbi_lineage_tsv),
        file(params.ncbi_genome_dir)
    )

    // Cruzamento das referências geradas com as reads limpas do Quality Control
    ch_align_input = FETCH_REFERENCE_GENOMES.out.ref_fasta.join(DEHOST_HISAT2.out.reads)

    ALIGN_AND_VARIANT_CALLING ( ch_align_input )
}

workflow.onComplete {
    if (workflow.success) {
        log.info """
          ===========================================
          ${ANSI_GREEN}Finished in ${workflow.duration}${ANSI_RESET}
          ${params.cleanup_work_on_success ? "Work directory cleaned up (cleanup_work_on_success=true)." : "Work directory kept (cleanup_work_on_success=false)."}
          ===========================================
        """.stripIndent()
    } else {
        log.info """
          ===========================================
          ${ANSI_RED}Finished with errors!${ANSI_RESET}
          Work directory preserved — you can retry with -resume.
          ===========================================
        """.stripIndent()
    }
}
