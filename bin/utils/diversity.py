#!/usr/bin/env python3
"""Módulo para cálculo de estatísticas de diversidade e relatórios usando Pandas e NumPy."""

from pathlib import Path
from typing import Tuple, List, Dict, Any
import numpy as np
import pandas as pd


def calculate_diversity_indices(
    csv_file_path: Path,
) -> Tuple[float, float, List[str]]:
    """Calcula o Índice de Simpson (D), Índice de Shannon (H) e extrai o Top 3 patógenos com Pandas."""
    if not csv_file_path.exists():
        return 0.0, 0.0, []

    dataframe = pd.read_csv(csv_file_path)

    if dataframe.empty or "RA" not in dataframe.columns:
        return 0.0, 0.0, []

    relative_abundances = dataframe["RA"].values

    # Cálculo vetorial com NumPy: Simpson Index (D = 1 - sum(p_i^2))
    simpson_index = round(float(1.0 - np.sum(relative_abundances**2)), 6)

    # Shannon Index (H = -sum(p_i * ln(p_i))), ignorando zeros
    positive_abundances = relative_abundances[relative_abundances > 0]
    shannon_index = (
        round(float(-np.sum(positive_abundances * np.log(positive_abundances))), 6)
        if len(positive_abundances) > 0
        else 0.0
    )

    # Extração das 3 primeiras linhas formatadas como string CSV
    top_three_dataframe = dataframe.head(3)
    top_three_lines = [
        f"{row['Microorganisms']},{int(row['reads'])},{row['RPM']:.2f},{row['RA']:.6f}"
        for _, row in top_three_dataframe.iterrows()
    ]

    return simpson_index, shannon_index, top_three_lines


def _collect_diversity_data_by_rank(
    output_directory: Path,
    sample_id: str,
    ranks_to_process: List[Tuple[str, str]],
) -> Dict[str, Dict[str, Any]]:
    """Função auxiliar: Itera pelos arquivos CSV dos ranks e calcula as métricas de diversidade."""
    diversity_results = {}

    for rank_label, rank_code in ranks_to_process:
        csv_file_path = output_directory / f"{sample_id}_p{rank_code}.csv"
        simpson_d, shannon_h, top_three_list = calculate_diversity_indices(
            csv_file_path
        )
        diversity_results[rank_code] = {
            "label": rank_label,
            "simpson": simpson_d,
            "shannon": shannon_h,
            "top_three": "\n".join(top_three_list),
        }

    return diversity_results


def _format_diversity_metrics_block(
    diversity_results: Dict[str, Dict[str, Any]],
    ranks_to_process: List[Tuple[str, str]],
    metric_key: str,
    metric_display_name: str,
) -> str:
    """Função auxiliar: Formata o bloco de texto para um índice específico (Simpson ou Shannon)."""
    formatted_lines = []
    for _, rank_code in ranks_to_process:
        label = diversity_results[rank_code]["label"]
        value = diversity_results[rank_code][metric_key]
        formatted_lines.append(f"{metric_display_name}({label})\t{value}")
    return "\n".join(formatted_lines)


def _format_top_pathogens_block(
    diversity_results: Dict[str, Dict[str, Any]],
    ranks_to_process: List[Tuple[str, str]],
) -> str:
    """Função auxiliar: Formata os blocos dos 3 principais organismos por rank."""
    formatted_blocks = []
    for _, rank_code in ranks_to_process:
        label = diversity_results[rank_code]["label"]
        top_three_content = diversity_results[rank_code]["top_three"]
        formatted_blocks.append(
            f"\n>Top_3hit_pathogen({label}) Name,reads,RPM,RA\n{top_three_content}"
        )
    return "\n".join(formatted_blocks)


def append_prompt_report_section(
    report_file_path: Path,
    output_directory: Path,
    sample_id: str,
):
    """Anexa os índices de diversidade e o Top 3 ao arquivo pre-report."""
    ranks_to_process = [
        ("family", "F"),
        ("genus", "G"),
        ("species", "S"),
    ]

    diversity_results = _collect_diversity_data_by_rank(
        output_directory, sample_id, ranks_to_process
    )

    simpson_block = _format_diversity_metrics_block(
        diversity_results, ranks_to_process, "simpson", "Simpson_Index"
    )
    shannon_block = _format_diversity_metrics_block(
        diversity_results, ranks_to_process, "shannon", "Shannon_Index"
    )
    top_pathogens_block = _format_top_pathogens_block(
        diversity_results, ranks_to_process
    )

    with open(report_file_path, "a", encoding="utf-8") as report_file:
        report_file.write("[PROMPT_REPORT]\n>Diversity_Index\n")
        report_file.write(f"{simpson_block}\n\n")
        report_file.write(f"{shannon_block}\n")
        report_file.write(f"{top_pathogens_block}\n")
