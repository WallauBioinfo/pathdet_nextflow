#!/usr/bin/env python3
"""
fetch_genomes.py
Replica a busca de genomas de referência do p05_map.sh (funções fetch/get_fasta):
para cada um dos 3 principais TaxIDs (<id>_merge_taxid.csv) procura primeiro na
tabela `reference/` e, se não houver, na tabela `all/` (refseq_lineage_taxid_fasta.tsv).

Saídas em <output_directory>:
  MAP/REF.fasta e MAP/ref_list{1,2,3}.txt   (só se algum genoma foi encontrado)
  ref_list{1,2,3}.txt na raiz                (se nenhum genoma foi encontrado, como o bash)
  <sample_id>_log_p5.txt

Uso: fetch_genomes.py <merge_taxid_csv> <ncbi_genome_dir> <output_directory> <sample_id>
"""

import gzip
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Set

# Resolve 'bin/utils' explicitamente antes de importar os utilitários
sys.path.append(str(Path(__file__).resolve().parent / "utils"))

from awk_compat import cut_field, read_lines  # noqa: E402


class LineageTable:
    """Tabela lineage<TAB>taxid<TAB>fasta_path, lida como o bash (grep sobre linhas)."""

    def __init__(self, table_path: Path):
        self.lines = read_lines(table_path)
        body = self.lines[1:]  # o bash faz `sed 1d` ao montar as listas de membros
        self.taxids: Set[str] = {cut_field(line, 2) for line in body}  # id_*1
        self.lineage_ids: Set[str] = {  # id_*2
            token for line in body for token in cut_field(line, 1).split(",")
        }

    def fetch(self, pattern: str) -> List[str]:
        """grep -E "pattern" tabela | cut -f 3."""
        regex = re.compile(pattern)
        return [cut_field(line, 3) for line in self.lines if regex.search(line)]

    def by_taxid(self, taxid: str) -> List[str]:
        return self.fetch(r"\s%s\s" % re.escape(taxid))

    def by_lineage(self, taxid: str) -> List[str]:
        return self.fetch(r"(^|,)%s(,|\s)" % re.escape(taxid))


def ordinal(number: int) -> str:
    return {1: "1st", 2: "2nd"}.get(number, "3rd")


def search_genomes(taxid: str, number: int, reference: LineageTable, registered: LineageTable):
    """Devolve (linhas_do_log, caminhos_brutos) para um dos 3 principais TaxIDs."""
    name = ordinal(number)
    log = [f">{name}_Hit microbe reference genome"]
    if taxid in reference.taxids:
        return log, reference.by_taxid(taxid)
    if taxid in reference.lineage_ids:
        return log, reference.by_lineage(taxid)
    log += ["Not applicable", f">{name}_Hit microbe registered genome"]
    if taxid in registered.taxids:
        return log, registered.by_taxid(taxid)
    if taxid in registered.lineage_ids:
        return log, registered.by_lineage(taxid)
    log += [
        "Not applicable",
        f"Automatic acquisition of {name} microbe complete genome failed.",
    ]
    return log, []


def main() -> None:
    if len(sys.argv) < 5:
        print(
            "Uso: fetch_genomes.py <merge_taxid_csv> <ncbi_genome_dir> "
            "<output_directory> <sample_id>",
            file=sys.stderr,
        )
        sys.exit(1)

    merge_taxid_path = Path(sys.argv[1])
    genome_dir = sys.argv[2].rstrip("/")
    output_dir = Path(sys.argv[3])
    sample_id = sys.argv[4]
    output_dir.mkdir(parents=True, exist_ok=True)

    reference = LineageTable(Path(genome_dir) / "reference" / "refseq_lineage_taxid_fasta.tsv")
    registered = LineageTable(Path(genome_dir) / "all" / "refseq_lineage_taxid_fasta.tsv")
    ranked = read_lines(merge_taxid_path)

    ref_lists, fasta_paths = [], []
    for number in (1, 2, 3):
        # sed -n ${num}p | cut -f 2
        taxid = cut_field(ranked[number - 1], 2) if len(ranked) >= number else ""
        log_lines, raw_paths = search_genomes(taxid, number, reference, registered)
        for raw in raw_paths:
            full = f"{genome_dir}/{raw.replace('db/ncbi_genome/', '')}"
            fasta_paths.append(full)
            log_lines.append(Path(full).name)
        ref_lists.append(log_lines)

    # REF.fasta (zcat >> REF.fasta) só se houve algum caminho
    map_dir = output_dir / "MAP"
    if fasta_paths:
        map_dir.mkdir(parents=True, exist_ok=True)
        with open(map_dir / "REF.fasta", "ab") as destination:
            for path in fasta_paths:
                try:
                    with gzip.open(path, "rb") as source:
                        shutil.copyfileobj(source, destination)
                except OSError as error:
                    print(f"zcat: {path}: {error}", file=sys.stderr)

    # <id>_log_p5.txt: carimbo + conteúdo dos três ref_list
    stamp = datetime.now().strftime("%y%m%d%H%M%S")
    log_text = f"{stamp} reference genome search finished.\n"
    for number, lines in enumerate(ref_lists, start=1):
        text = "".join(f"{line}\n" for line in lines)
        log_text += text
        # mv ref_list*.txt MAP/ só funciona se MAP existir
        target = map_dir if fasta_paths else output_dir
        (target / f"ref_list{number}.txt").write_text(text, encoding="utf-8")
    (output_dir / f"{sample_id}_log_p5.txt").write_text(log_text, encoding="utf-8")


if __name__ == "__main__":
    main()
