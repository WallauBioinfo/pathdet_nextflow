#!/usr/bin/env python3
"""
fetch_genomes.py
Busca os genomas de referência para os 3 principais TaxIDs identificados na etapa Rapid/BLAST
utilizando a tabela pré-gerada refseq_lineage_taxid_fasta.tsv.
"""

import sys
import os
import gzip
from pathlib import Path
from typing import List, TextIO
import pandas as pd

# Garante a resolução explícita para módulos contidos na pasta bin/utils
sys.path.append(str(Path(__file__).resolve().parent / "utils"))


def find_fasta_paths_for_single_taxid(
    target_taxid: str, 
    reference_lineage_dataframe: pd.DataFrame, 
    database_base_directory: Path
) -> List[Path]:
    """Localiza os caminhos dos arquivos FASTA associados a um TaxID específico."""
    matching_fasta_paths: List[Path] = []
    
    # Busca exata por TaxID na coluna tax_id
    matching_rows = reference_lineage_dataframe[
        reference_lineage_dataframe['tax_id'].astype(str) == target_taxid
    ]
    
    # Se não encontrar, busca se o TaxID está presente na linhagem taxonômica
    if matching_rows.empty:
        matching_rows = reference_lineage_dataframe[
            reference_lineage_dataframe['lineage'].astype(str).str.contains(rf'(^|,){target_taxid}(,|$)')
        ]

    if not matching_rows.empty:
        for raw_fasta_path in matching_rows['fasta_path']:
            sanitized_path = str(raw_fasta_path).replace("db/ncbi_genome/", "")
            full_genome_path = database_base_directory / sanitized_path
            
            if full_genome_path.exists():
                matching_fasta_paths.append(full_genome_path)
            elif Path(raw_fasta_path).exists():
                matching_fasta_paths.append(Path(raw_fasta_path))
                
    return matching_fasta_paths


def get_ordinal_description(rank_position: int) -> str:
    """Retorna a representação ordinal para uma posição numérica."""
    ordinal_mapping = {1: "1st", 2: "2nd", 3: "3rd"}
    return ordinal_mapping.get(rank_position, f"{rank_position}th")


def process_single_taxid_reference(
    rank_position: int,
    target_taxid: str,
    reference_lineage_dataframe: pd.DataFrame,
    database_base_directory: Path,
    log_file_handle: TextIO
) -> List[Path]:
    """Processa um único TaxID, escrevendo no log e retornando os caminhos FASTA encontrados."""
    ordinal_description = get_ordinal_description(rank_position)
    log_file_handle.write(f">{ordinal_description}_Hit microbe reference genome ({target_taxid})\n")
    
    found_fasta_paths = find_fasta_paths_for_single_taxid(
        target_taxid, 
        reference_lineage_dataframe, 
        database_base_directory
    )
    
    if found_fasta_paths:
        for fasta_file_path in found_fasta_paths:
            log_file_handle.write(f"{fasta_file_path.name}\n")
        return found_fasta_paths
    
    log_file_handle.write("Not applicable\n")
    log_file_handle.write(f"Automatic acquisition of {ordinal_description} microbe complete genome failed.\n")
    return []


def collect_fasta_paths_for_top_taxids(
    top_taxid_list: List[str],
    reference_lineage_dataframe: pd.DataFrame,
    database_base_directory: Path,
    log_file_path: Path
) -> List[Path]:
    """Itera sobre a lista de principais TaxIDs registrando logs e acumulando os FASTAs."""
    accumulated_fasta_paths: List[Path] = []
    
    with open(log_file_path, "w", encoding="utf-8") as log_file_handle:
        for rank_position, target_taxid in enumerate(top_taxid_list, start=1):
            fasta_paths_for_taxid = process_single_taxid_reference(
                rank_position,
                target_taxid,
                reference_lineage_dataframe,
                database_base_directory,
                log_file_handle
            )
            accumulated_fasta_paths.extend(fasta_paths_for_taxid)
            
    return accumulated_fasta_paths


def concatenate_compressed_fastas(
    source_fasta_paths: List[Path], 
    output_fasta_path: Path
) -> None:
    """Concatena arquivos .fna.gz diretamente via fluxo de bytes."""
    with open(output_fasta_path, "wb") as output_file_stream:
        for fasta_file_path in source_fasta_paths:
            with gzip.open(fasta_file_path, "rb") as gzip_input_stream:
                output_file_stream.write(gzip_input_stream.read())


def main() -> None:
    if len(sys.argv) < 5:
        print(
            "Uso: python3 fetch_genomes.py <merge_taxid_csv> <ncbi_lineage_tsv> <database_base_directory> <output_directory>", 
            file=sys.stderr
        )
        sys.exit(1)

    taxid_csv_path = Path(sys.argv[1])
    lineage_tsv_path = Path(sys.argv[2])
    database_base_directory = Path(sys.argv[3])
    output_directory = Path(sys.argv[4])
    
    output_directory.mkdir(parents=True, exist_ok=True)
    destination_fasta_path = output_directory / "REF.fasta"
    reference_log_path = output_directory / "ref_list.txt"

    if not taxid_csv_path.exists() or os.stat(taxid_csv_path).st_size == 0:
        print(f"Aviso: Arquivo {taxid_csv_path} está vazio ou não foi encontrado.", file=sys.stderr)
        destination_fasta_path.touch()
        sys.exit(0)

    # 1. Carrega os Top 3 TaxIDs do arquivo de resultado do BLAST/Krona
    # Usa sep=None com o motor do python para inferir automaticamente se é delimitado por TSV, CSV ou Espaço
    top_hits_dataframe = pd.read_csv(
        taxid_csv_path, 
        sep=None, 
        engine='python', 
        header=None,
        on_bad_lines='skip'
    )
    
    top_taxid_list = (
        top_hits_dataframe.iloc[:3, 1].astype(str).tolist() 
        if not top_hits_dataframe.empty and top_hits_dataframe.shape[1] > 1 
        else []
    )

    # 2. Carrega a tabela TSV pré-processada de genomas
    reference_lineage_dataframe = pd.read_csv(
        lineage_tsv_path, 
        sep='\t', 
        header=None, 
        names=['lineage', 'tax_id', 'fasta_path']
    )

    # 3. Processa a busca de genomas e gera o log
    collected_fasta_paths = collect_fasta_paths_for_top_taxids(
        top_taxid_list,
        reference_lineage_dataframe,
        database_base_directory,
        reference_log_path
    )

    # 4. Concatenação direta de bytes (I/O puro)
    if collected_fasta_paths:
        concatenate_compressed_fastas(collected_fasta_paths, destination_fasta_path)
    else:
        destination_fasta_path.touch()


if __name__ == '__main__':
    main()
