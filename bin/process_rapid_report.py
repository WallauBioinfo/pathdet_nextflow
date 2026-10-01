#!/usr/bin/env python3
"""Etapa Rapid: tudo o que o rapid() do bash fazia depois do kraken2 e do Krona.

Gera, com os mesmos nomes e conteúdos do bash:
  KK2/<id>_p2.kreport, <id>_tax.kreport, <id>_tax2.kreport
  CSV/<id>_mp{K,P,F,G,S}.csv e <id>_rp{K,P,F,G,S}.csv
  <id>_pRep/<id>_p{K,P,F,G,S}.csv
  <id>_pre-report.csv e <id>_log_Rapid.txt
"""

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

# Resolve 'bin/utils' explicitamente antes de importar os utilitários
sys.path.append(str(Path(__file__).resolve().parent / "utils"))

from awk_compat import bc_percent, read_lines, write_lines  # noqa: E402
from diversity import head3, shannon_index, simpson_index  # noqa: E402
from log_parser import calculate_pipeline_metrics, count_fastq_reads  # noqa: E402
from taxonomy import (  # noqa: E402
    build_rank_table,
    kreport_to_p2,
    p2_to_tax_lines,
    tax_to_rpm_lines,
    taxonkit_lineage,
    with_header,
)

RANKS = "kpfgs"
_START = time.time()


def _timestamp() -> str:
    return datetime.now().strftime("%y%m%d%H%M%S")


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Relatórios da etapa Rapid.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--log-p1-path", required=True, type=Path)
    parser.add_argument("--log-p2-path", required=True, type=Path)
    parser.add_argument("--kreport-path", required=True, type=Path,
                        help="<id>_p1.kreport gerado pelo kraken2 (--report)")
    parser.add_argument("--classified-fastq-path", required=True, type=Path)
    parser.add_argument("--unclassified-fastq-path", required=True, type=Path)
    parser.add_argument("--kk2-dir", required=True, type=Path)
    parser.add_argument("--csv-dir", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path,
                        help="<id>_pRep")
    parser.add_argument("--output-report-path", required=True, type=Path)
    parser.add_argument("--output-log-path", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=1)
    return parser.parse_args()


def _pre_report_text(sample_id, metrics, diversity, top3) -> str:
    return (
        f"Sample_ID\t{sample_id}\n"
        f"Total_read(reads)\t{metrics['total_raw_reads']}\n"
        f"Available_read(reads)\t{metrics['total_quality_filtered_reads']}\n"
        f"Non-human_read(reads)\t{metrics['total_non_human_reads']}\n"
        f"Non-human_read(RPM)\t{metrics['reads_per_million_non_human']}\n"
        "\n[PROMPT_REPORT]\n"
        ">Diversity_Index\n"
        f"Simpson_Index(family)\t{diversity['simpson']['f']}\n"
        f"Simpson_Index(genus)\t{diversity['simpson']['g']}\n"
        f"Simpson_Index(species)\t{diversity['simpson']['s']}\n"
        f"\nShannon_Index(family)\t{diversity['shannon']['f']}\n"
        f"Shannon_Index(genus)\t{diversity['shannon']['g']}\n"
        f"Shannon_Index(species)\t{diversity['shannon']['s']}\n"
        f"\n>Top_3hit_pathogen(family) Name,reads,RPM,RA\n{top3['f']}\n"
        f"\n>Top_3hit_pathogen(genus) Name,reads,RPM,RA\n{top3['g']}\n"
        f"\n>Top_3hit_pathogen(species) Name,reads,RPM,RA\n{top3['s']}\n"
    )


def _elapsed() -> str:
    seconds = int(time.time() - _START)
    return f"{seconds // 3600}:{(seconds % 3600) // 60}:{seconds % 60}"


def main() -> None:
    args = _parse_arguments()
    sid = args.sample_id
    for directory in (args.kk2_dir, args.csv_dir, args.output_directory):
        directory.mkdir(parents=True, exist_ok=True)

    metrics = calculate_pipeline_metrics(args.log_p1_path, args.log_p2_path)
    wca1 = metrics["total_raw_reads"]
    rmda2 = metrics["total_quality_filtered_reads"]
    hg2 = metrics["total_non_human_reads"]

    # --- contagens do kraken2 e log_Rapid (primeiras 2 linhas)
    kraken_stamp = _timestamp()
    classified = count_fastq_reads(args.classified_fastq_path)
    unclassified = count_fastq_reads(args.unclassified_fastq_path)
    log_text = (
        f"{kraken_stamp} kraken2 finished: {classified}reads "
        f"({bc_percent(classified, hg2)}% of karaken2 target)\n"
        f"\t* un-classified read: {unclassified}reads "
        f"({bc_percent(unclassified, wca1)}%)\n"
    )

    # --- KK2/<id>_p2.kreport, _tax.kreport, _tax2.kreport
    p2_lines = kreport_to_p2(read_lines(args.kreport_path))
    write_lines(args.kk2_dir / f"{sid}_p2.kreport", p2_lines)
    tax_lines = taxonkit_lineage(p2_to_tax_lines(p2_lines), 2)
    write_lines(args.kk2_dir / f"{sid}_tax.kreport", tax_lines)
    tax2_lines = tax_to_rpm_lines(tax_lines, rmda2)
    write_lines(args.kk2_dir / f"{sid}_tax2.kreport", tax2_lines)

    # --- CSV/ e pRep/ por rank (o kingdom não filtra nomes vazios, como no bash)
    mp_by_rank, final_by_rank = {}, {}
    for rank in RANKS:
        rp, mp, final = build_rank_table(
            tax_lines, tax2_lines, rank, args.threads, drop_empty_names=(rank != "k")
        )
        upper = rank.upper()
        write_lines(args.csv_dir / f"{sid}_mp{upper}.csv", mp)
        write_lines(args.csv_dir / f"{sid}_rp{upper}.csv", rp)
        write_lines(args.output_directory / f"{sid}_p{upper}.csv", with_header(final))
        mp_by_rank[rank], final_by_rank[rank] = mp, final

    # --- pre-report
    diversity = {
        "simpson": {r: simpson_index(mp_by_rank[r]) for r in "fgs"},
        "shannon": {r: shannon_index(mp_by_rank[r]) for r in "fgs"},
    }
    top3 = {r: head3(final_by_rank[r]) for r in "fgs"}
    args.output_report_path.write_text(
        _pre_report_text(sid, metrics, diversity, top3), encoding="utf-8"
    )

    # --- log_Rapid (3.ª linha)
    log_text += f"{_timestamp()} PROMPT REPORT finished ({_elapsed()})\n"
    args.output_log_path.write_text(log_text, encoding="utf-8")


if __name__ == "__main__":
    main()
