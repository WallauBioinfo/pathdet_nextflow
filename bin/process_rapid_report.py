#!/usr/bin/env python3
"""Script orquestrador principal da etapa Rapid."""

import argparse
import sys
from pathlib import Path

# Resolve o caminho para 'bin/utils' explicitamente antes de importar os utilitários
sys.path.append(str(Path(__file__).resolve().parent / "utils"))

from diversity import append_prompt_report_section
from log_parser import (
    calculate_pipeline_metrics,
    count_fastq_reads,
)
from taxonomy import (
    aggregate_counts_by_rank,
    parse_kraken_output,
    write_rank_csv_reports,
)


def _parse_arguments() -> argparse.Namespace:
    """Configura e lê os argumentos de linha de comando."""
    parser = argparse.ArgumentParser(
        description="Processa relatórios e métricas de taxonomia para a etapa Rapid."
    )
    parser.add_argument(
        "--sample-id", required=True, help="Identificador único da amostra"
    )
    parser.add_argument(
        "--log-p1-path",
        required=True,
        type=Path,
        help="Caminho para o log da etapa p1 (QC)",
    )
    parser.add_argument(
        "--log-p2-path",
        required=True,
        type=Path,
        help="Caminho para o log da etapa p2 (Dehosting)",
    )
    parser.add_argument(
        "--kraken-output-path",
        required=True,
        type=Path,
        help="Caminho para a saída kraken_p1.txt do Kraken2",
    )
    parser.add_argument(
        "--classified-fastq-path",
        required=True,
        type=Path,
        help="Caminho do FASTQ de leituras classificadas pelo Kraken2",
    )
    parser.add_argument(
        "--unclassified-fastq-path",
        required=True,
        type=Path,
        help="Caminho do FASTQ de leituras não-classificadas pelo Kraken2",
    )
    parser.add_argument(
        "--output-directory",
        required=True,
        type=Path,
        help="Diretório onde os CSVs de patente taxonômica serão salvos",
    )
    parser.add_argument(
        "--output-report-path",
        required=True,
        type=Path,
        help="Caminho de saída do pre-report consolidado",
    )
    parser.add_argument(
        "--output-log-path",
        required=True,
        type=Path,
        help="Caminho de saída do log final da etapa Rapid",
    )

    return parser.parse_args()


def _write_initial_pre_report(
    report_path: Path, sample_id: str, pipeline_metrics: dict
):
    """Escreve o cabeçalho e os dados iniciais do pipeline no relatório."""
    with open(report_path, "w", encoding="utf-8") as report_file:
        report_file.write(f"Sample_id,{sample_id}\n\n")
        report_file.write("[FASTQ_STATUS]\n")
        report_file.write(
            f"Total_raw_reads\t{pipeline_metrics['total_raw_reads']}\n"
        )
        report_file.write(
            "Total_quality_filtered_reads\t"
            f"{pipeline_metrics['total_quality_filtered_reads']}\n"
        )
        report_file.write(
            "Total_non_human_reads\t"
            f"{pipeline_metrics['total_non_human_reads']}\n"
        )
        report_file.write(
            "Reads_per_million_non_human\t"
            f"{pipeline_metrics['reads_per_million_non_human']}\n\n"
        )


def _write_kraken_status_section(
    report_path: Path, classified_reads: int, unclassified_reads: int
):
    """Anexa a seção de contagens do Kraken2 ao relatório consolidado."""
    with open(report_path, "a", encoding="utf-8") as report_file:
        report_file.write("[KRAKEN2_STATUS]\n")
        report_file.write(f"Classified_reads\t{classified_reads}\n")
        report_file.write(f"Unclassified_reads\t{unclassified_reads}\n\n")


def _write_execution_log(
    log_path: Path,
    sample_id: str,
    pipeline_metrics: dict,
    classified_reads: int,
    unclassified_reads: int,
):
    """Gera o arquivo de log legível de finalização do módulo Rapid."""
    with open(log_path, "w", encoding="utf-8") as log_file:
        log_file.write(f"=== Log de Execução Rapid - Amostra: {sample_id} ===\n\n")
        log_file.write(
            f"Reads Brutas: {pipeline_metrics['total_raw_reads']}\n"
        )
        log_file.write(
            "Reads Filtradas por Qualidade: "
            f"{pipeline_metrics['total_quality_filtered_reads']}\n"
        )
        log_file.write(
            "Reads Não-Humanas: "
            f"{pipeline_metrics['total_non_human_reads']}\n"
        )
        log_file.write(
            "Classificadas no Kraken2: "
            f"{classified_reads}\n"
        )
        log_file.write(
            "Não Classificadas no Kraken2: "
            f"{unclassified_reads}\n"
        )


def main():
    args = _parse_arguments()

    # Garantir existência do diretório de saída
    args.output_directory.mkdir(parents=True, exist_ok=True)

    # 1. Leitura e cálculo de métricas básicas a partir dos logs do p1 e p2
    pipeline_metrics = calculate_pipeline_metrics(
        args.log_p1_path, args.log_p2_path
    )

    # 2. Contagem das reads FASTQ classificadas/não-classificadas
    classified_reads = count_fastq_reads(args.classified_fastq_path)
    unclassified_reads = count_fastq_reads(args.unclassified_fastq_path)

    # 3. Geração inicial do relatório consolidado
    _write_initial_pre_report(
        args.output_report_path, args.sample_id, pipeline_metrics
    )
    _write_kraken_status_section(
        args.output_report_path, classified_reads, unclassified_reads
    )

    # 4. Agregação taxonômica (TaxIDs -> Níveis de linhagem NCBI)
    taxid_counts = parse_kraken_output(args.kraken_output_path)
    aggregated_rank_data = aggregate_counts_by_rank(taxid_counts)

    # 5. Exportação dos relatórios CSV por rank via Pandas (RPM e RA)
    write_rank_csv_reports(
        aggregated_rank_data,
        pipeline_metrics["total_quality_filtered_reads"],
        args.output_directory,
        args.sample_id,
    )

    # 6. Cálculo dos índices de diversidade e Top 3 patógenos no relatório
    append_prompt_report_section(
        args.output_report_path, args.output_directory, args.sample_id
    )

    # 7. Escrita do arquivo final de log
    _write_execution_log(
        args.output_log_path,
        args.sample_id,
        pipeline_metrics,
        classified_reads,
        unclassified_reads,
    )


if __name__ == "__main__":
    main()