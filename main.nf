nextflow.enable.dsl = 2

ANSI_GREEN = "\033[1;32m"
ANSI_RED   = "\033[1;31m"
ANSI_RESET = "\033[0m"

include { QUALITY_CONTROL } from './modules/quality_control.nf'
include { DEHOST_KRAKEN2  } from './modules/dehosting.nf'
include { DEHOST_BOWTIE2  } from './modules/dehosting.nf'
include { DEHOST_HISAT2   } from './modules/dehosting.nf'
include { DEHOST_BLASTN   } from './modules/dehosting.nf'
include { RAPID_ANALYSIS } from './modules/rapid.nf'
include {
    SPLIT_FASTA;
    RUN_BLASTN;
    MERGE_AND_KRONA
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

    // 1. Quality control (p1)
    QUALITY_CONTROL(input_csv)

    // 2. Dehosting (p2)
    DEHOST_KRAKEN2 ( QUALITY_CONTROL.out.sample, params.kraken2_grch38 )
    DEHOST_BOWTIE2 ( DEHOST_KRAKEN2.out.reads,   params.bowtie2_host )
    DEHOST_HISAT2  ( DEHOST_BOWTIE2.out.reads,   params.hisat2_grch38r )

    // Join by sample_id: with multiple samples, channels may not arrive in the same order
    DEHOST_BLASTN  ( DEHOST_HISAT2.out.reads.join(QUALITY_CONTROL.out.log), params.blast_host2 )

    // 3. Rapid Analysis (p03_rapid)
    RAPID_ANALYSIS (
        DEHOST_BLASTN.out.fasta
            .join(QUALITY_CONTROL.out.log)
            .join(DEHOST_BLASTN.out.log),
        params.pathdb1_db
    )

    // 4. Module BLAST (p03_blast) and Krona (p04_krona)
    SPLIT_FASTA ( DEHOST_BLASTN.out.fasta )

    // blastn against nt costs about one pass over the database no matter how few reads
    // there are, so chunks from ALL samples are searched together, params.blast_batch_chunks
    // chunks per blastn run. The result is split back per chunk inside RUN_BLASTN, so
    // everything downstream still receives one <chunk>_blast.txt per chunk, per sample.
    // Sorting by chunk name keeps the batches identical between runs (-resume friendly).
    ch_blast_batches = SPLIT_FASTA.out.chunks
        .transpose()                                         // [sample_id, chunk]
        .toSortedList { a, b -> a[1].name <=> b[1].name }
        .flatMap { it }
        .collate( params.blast_batch_chunks as int )
        .map { group -> tuple( "batch_" + group[0][1].baseName, group.collect { it[1] } ) }

    RUN_BLASTN ( ch_blast_batches )

    // <sample_id>_<NN>_blast.txt  ->  [sample_id, [files...]]
    ch_grouped_blast = RUN_BLASTN.out.blast_files
        .flatten()
        .map { f -> tuple( f.name.replaceFirst(/_\d+_blast\.txt$/, ''), f ) }
        .groupTuple()
    MERGE_AND_KRONA ( ch_grouped_blast.join(QUALITY_CONTROL.out.log) )

    // 5. Mapping, Coverage, and Variant Calling (p05_map)
    FETCH_REFERENCE_GENOMES (
        MERGE_AND_KRONA.out.merge_taxid,
        file(params.ncbi_genome_dir)
    )

    ch_align_input = FETCH_REFERENCE_GENOMES.out.ref_fasta.join(QUALITY_CONTROL.out.sample)
    ALIGN_AND_VARIANT_CALLING ( ch_align_input )
}

workflow.onComplete {
    if (workflow.success) {
        log.info """
          ===========================================
          ${ANSI_GREEN}Finished in ${workflow.duration}${ANSI_RESET}
          ${params.cleanup_work_on_success ? "Work directory cleaned up (cleanup_work_on_success=true)." : "Work directory kept."}
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
