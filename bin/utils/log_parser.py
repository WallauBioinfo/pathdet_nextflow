#!/usr/bin/env python3
"""Métricas lidas dos logs de QC/HGS, com a mesma semântica do bash original."""

from pathlib import Path
from typing import Dict

from awk_compat import bc_rpm, count_newlines


def extract_read_count_from_log(log_file_path: Path, line_index: int) -> int:
    """awk 'NR==N {print $4}' | tr -d ',' (line_index começa em 0)."""
    with open(log_file_path, "r", encoding="utf-8") as file_handle:
        log_lines = file_handle.readlines()
    return int(log_lines[line_index].split()[3].replace(",", ""))


def count_fastq_reads(fastq_file_path: Path) -> int:
    """cat ficheiro | wc -l; expr $n / 4 (igual ao bash, mesmo com FASTA de 2 linhas)."""
    return count_newlines(fastq_file_path) // 4


def calculate_pipeline_metrics(log_p1_path: Path, log_p2_path: Path) -> Dict[str, object]:
    """WCA1, RMDA2, HG2 e PRstd do rapid() do bash.

    * Total e Available: nº de reads do R1 vezes 2 (linhas 3 e 8 do log_p1);
    * Non-human: linha 8 do log_p2, SEM multiplicar por 2 (como no bash);
    * RPM: bc com scale=0 (trunca).
    """
    total_raw_reads = extract_read_count_from_log(log_p1_path, 2) * 2
    total_quality_filtered_reads = extract_read_count_from_log(log_p1_path, 7) * 2
    total_non_human_reads = extract_read_count_from_log(log_p2_path, 7)
    return {
        "total_raw_reads": total_raw_reads,
        "total_quality_filtered_reads": total_quality_filtered_reads,
        "total_non_human_reads": total_non_human_reads,
        "reads_per_million_non_human": bc_rpm(
            total_non_human_reads, total_quality_filtered_reads
        ),
    }
