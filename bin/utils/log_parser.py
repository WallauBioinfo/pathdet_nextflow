#!/usr/bin/env python3

from pathlib import Path
from typing import Dict


def extract_read_count_from_log(log_file_path: Path, line_index: int) -> int:
    with open(log_file_path, "r", encoding="utf-8") as file_handle:
        log_lines = file_handle.readlines()
        raw_read_value = log_lines[line_index].split()[3].replace(",", "")
        return int(raw_read_value)


def count_fastq_reads(fastq_file_path: Path) -> int:
    if not fastq_file_path.exists():
        return 0
    with open(fastq_file_path, "r", encoding="utf-8") as file_handle:
        total_lines = sum(1 for _ in file_handle)
        return total_lines // 4


def calculate_pipeline_metrics(log_p1_path: Path, log_p2_path: Path) -> Dict[str, int]:
    # Leitura das contagens dos logs anteriores
    total_raw_reads = extract_read_count_from_log(log_p1_path, line_index=2) * 2
    total_quality_filtered_reads = (
        extract_read_count_from_log(log_p1_path, line_index=7) * 2
    )
    total_non_human_reads = (
        extract_read_count_from_log(log_p2_path, line_index=7) * 2
    )

    # Cálculo da métrica de Reads Per Million (RPM) de sequências não-humanas
    reads_per_million_non_human = (
        round((total_non_human_reads * 1_000_000) / total_quality_filtered_reads)
        if total_quality_filtered_reads > 0
        else 0
    )

    return {
        "total_raw_reads": total_raw_reads,
        "total_quality_filtered_reads": total_quality_filtered_reads,
        "total_non_human_reads": total_non_human_reads,
        "reads_per_million_non_human": reads_per_million_non_human,
    }
