#!/usr/bin/env python3
"""
taxid_parser.py

Maps local reference genomes (assembly report + fna.gz file)
based on the taxonomy tree (nodes.dmp) provided by NCBI.
"""

import sys
import re
from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd


ORGANISM_GROUPS = ['bacteria', 'fungi', 'viral']


def load_ncbi_parent_map(nodes_file_path: str) -> Dict[int, int]:
    dataframe_nodes = pd.read_csv(
        nodes_file_path,
        sep=r'\s*\|\s*',
        engine='python',
        header=None,
        usecols=[0, 2],
        names=['child_tax_id', 'parent_tax_id']
    )
    return dict(zip(dataframe_nodes['child_tax_id'], dataframe_nodes['parent_tax_id']))


def build_taxonomy_lineage(child_tax_id: int, parent_taxid_map: Dict[int, int]) -> str:
    lineage_ids: List[str] = [str(child_tax_id)]
    current_tax_id = child_tax_id

    while current_tax_id in parent_taxid_map:
        current_tax_id = parent_taxid_map[current_tax_id]
        lineage_ids.append(str(current_tax_id))
        if current_tax_id == 1:
            break

    return ",".join(lineage_ids)


def extract_taxid_from_assembly_report(assembly_report_path: Path) -> Optional[int]:
    if not assembly_report_path.is_file():
        return None

    with open(assembly_report_path, 'r', encoding='utf-8', errors='ignore') as report_file:
        for line in report_file:
            cleaned_line = line.strip()
            if cleaned_line.startswith("# Taxid:"):
                match = re.search(r"^# Taxid:\s+(\d+)$", cleaned_line)
                if match:
                    return int(match.group(1))

    return None


def parse_accession_directory(accession_directory: Path, group_name: str) -> Optional[dict]:
    accession_number = accession_directory.name
    report_files = list(accession_directory.glob(f"{accession_number}_*_assembly_report.txt"))

    if not report_files:
        return None

    assembly_report_path = report_files[0]
    assembly_name_match = re.search(
        rf"^{re.escape(accession_number)}_(.+)_assembly_report\.txt$", 
        assembly_report_path.name
    )
    
    if not assembly_name_match:
        return None

    assembly_name = assembly_name_match.group(1)
    fasta_file_path = accession_directory / f"{accession_number}_{assembly_name}_genomic.fna.gz"

    tax_id = extract_taxid_from_assembly_report(assembly_report_path)
    if tax_id is None:
        return None

    return {
        "tax_id": tax_id,
        "accession": accession_number,
        "group": group_name,
        "fasta_path": str(fasta_file_path)
    }


def collect_genomes_dataframe(genomes_base_directory: str) -> pd.DataFrame:
    base_path = Path(genomes_base_directory)
    records: List[dict] = []

    for group_name in ORGANISM_GROUPS:
        group_directory = base_path / group_name
        if not group_directory.is_dir():
            continue

        for accession_directory in group_directory.iterdir():
            if accession_directory.is_dir():
                metadata = parse_accession_directory(accession_directory, group_name)
                if metadata:
                    records.append(metadata)

    if not records:
        return pd.DataFrame(columns=["tax_id", "accession", "group", "fasta_path"])

    return pd.DataFrame(records)


def process_taxonomy_and_format_output(
    dataframe_genomes: pd.DataFrame, 
    parent_taxid_map: Dict[int, int]
) -> pd.DataFrame:
    if dataframe_genomes.empty:
        return dataframe_genomes

    dataframe_genomes = dataframe_genomes.drop_duplicates(subset=["tax_id", "accession"]).copy()
    dataframe_genomes["lineage"] = dataframe_genomes["tax_id"].apply(
        lambda tax_id: build_taxonomy_lineage(tax_id, parent_taxid_map)
    )
    dataframe_genomes = dataframe_genomes.sort_values(by="tax_id")
    return dataframe_genomes[["lineage", "tax_id", "fasta_path"]]


def main() -> None:
    if len(sys.argv) < 3:
        print("Use: python3 taxid_parser.py <genomes_directory> <nodes.dmp>", file=sys.stderr)
        sys.exit(1)

    genomes_base_directory = sys.argv[1]
    ncbi_nodes_file_path = sys.argv[2]
    parent_taxid_map = load_ncbi_parent_map(ncbi_nodes_file_path)
    dataframe_genomes = collect_genomes_dataframe(genomes_base_directory)
    dataframe_output = process_taxonomy_and_format_output(dataframe_genomes, parent_taxid_map)

    if not dataframe_output.empty:
        dataframe_output.to_csv(sys.stdout, sep="\t", header=False, index=False)


if __name__ == '__main__':
    main()