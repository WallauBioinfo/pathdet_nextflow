#!/usr/bin/env python3
"""Módulo para processamento taxonômico usando a biblioteca ete3 e exportação via Pandas."""

from pathlib import Path
from typing import Dict, Tuple
import pandas as pd
from ete3 import NCBITaxa


def parse_kraken_output(kraken_output_path: Path) -> Dict[int, int]:
    """Mapeia os TaxIDs classificados pelo Kraken2 para suas respectivas frequências.

    Lê a saída kraken_p1.txt e filtra apenas as sequências classificadas ('C').
    """
    taxid_read_counts = {}
    with open(kraken_output_path, "r", encoding="utf-8") as file_handle:
        for line in file_handle:
            columns = line.strip().split("\t")
            if len(columns) >= 3 and columns[0] == "C":
                taxonomy_id = int(columns[2])
                if taxonomy_id > 0:
                    taxid_read_counts[taxonomy_id] = (
                        taxid_read_counts.get(taxonomy_id, 0) + 1
                    )
    return taxid_read_counts


def _map_taxonomy_id_to_ancestors(
    taxonomy_id: int,
    read_count: int,
    ncbi_database: NCBITaxa,
    rank_code_mapping: Dict[str, str],
    aggregated_rank_data: Dict[str, Dict[str, int]],
):
    """Função auxiliar: Mapeia os ancestrais de um TaxID e incrementa as contagens acumuladas."""
    try:
        lineage_ids = ncbi_database.get_lineage(taxonomy_id)
        rank_mapping = ncbi_database.get_rank(lineage_ids)
        translated_names = ncbi_database.get_taxid_translator(lineage_ids)

        for ancestor_id in lineage_ids:
            rank_name = rank_mapping.get(ancestor_id)
            if rank_name in rank_code_mapping:
                short_rank_code = rank_code_mapping[rank_name]
                organism_name = translated_names.get(
                    ancestor_id, f"taxid_{ancestor_id}"
                ).replace(" ", "_")

                current_count = aggregated_rank_data[short_rank_code].get(
                    organism_name, 0
                )
                aggregated_rank_data[short_rank_code][organism_name] = (
                    current_count + read_count
                )
    except Exception:
        # Ignora TaxIDs não encontrados no banco do ete3
        pass


def aggregate_counts_by_rank(
    taxid_counts: Dict[int, int]
) -> Dict[str, Dict[str, int]]:
    """Agrupa a contagem de reads por patente taxonômica (Reino, Filo, Família, Gênero, Espécie)."""
    ncbi_database = NCBITaxa()

    rank_code_mapping = {
        "superkingdom": "K",
        "phylum": "P",
        "family": "F",
        "genus": "G",
        "species": "S",
    }
    aggregated_rank_data = {rank_code: {} for rank_code in rank_code_mapping.values()}

    for taxonomy_id, read_count in taxid_counts.items():
        _map_taxonomy_id_to_ancestors(
            taxonomy_id,
            read_count,
            ncbi_database,
            rank_code_mapping,
            aggregated_rank_data,
        )

    return aggregated_rank_data


def _build_rank_dataframe(
    organism_counts: Dict[str, int], total_quality_filtered_reads: int
) -> pd.DataFrame:
    """Função auxiliar: Cria e preenche um DataFrame Pandas para uma única patente taxonômica."""
    if not organism_counts:
        return pd.DataFrame(columns=["Microorganisms", "reads", "RPM", "RA"])

    dataframe = pd.DataFrame(
        list(organism_counts.items()),
        columns=["Microorganisms", "reads"],
    )

    dataframe = dataframe.sort_values(by="reads", ascending=False)
    total_rank_reads = dataframe["reads"].sum()

    dataframe["RPM"] = (
        (dataframe["reads"] / total_quality_filtered_reads * 1_000_000).round(2)
        if total_quality_filtered_reads > 0
        else 0.0
    )

    dataframe["RA"] = (
        (dataframe["reads"] / total_rank_reads).round(6)
        if total_rank_reads > 0
        else 0.0
    )

    return dataframe


def write_rank_csv_reports(
    aggregated_rank_data: Dict[str, Dict[str, int]],
    total_quality_filtered_reads: int,
    output_directory: Path,
    sample_id: str,
):
    """Gera os relatórios CSV para cada patente taxonômica."""
    for rank_code, organism_counts in aggregated_rank_data.items():
        rank_dataframe = _build_rank_dataframe(
            organism_counts, total_quality_filtered_reads
        )
        csv_output_path = output_directory / f"{sample_id}_p{rank_code}.csv"
        rank_dataframe.to_csv(csv_output_path, index=False)
