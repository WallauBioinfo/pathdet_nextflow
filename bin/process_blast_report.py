#!/usr/bin/env python3
"""
Step p04 (BLAST tables): reproduces p04_krona.sh from
<id>_merge_tax.csv (reads + lineage) and <id>_merge_tax2.csv (RPM + lineage).

Generates, with the exact same names and contents as the bash script:
  <id>_CSV2/<id>_{mf,rf}{K,P,F,G,S}.csv, _mBacteria*/_rBacteria*, _mVirus/_rVirus, _mFungi/_rFungi
  <id>_fRep/<id>_f{K,P,F,G,S}.csv, _Virus.csv, _Fungi.csv
  <id>_bRep/<id>_Bacteria_f{F,G,S}.csv
  <id>_report.csv and <id>_bacteria_report.csv
"""

import argparse
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent / "utils"))

from awk_compat import awk_sum_second_column, read_lines, write_lines  # noqa: E402
from diversity import head3, shannon_index, simpson_index  # noqa: E402
from taxonomy import build_rank_table, with_header  # noqa: E402


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="p04 reports (BLAST).")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--merge-tax", required=True, type=Path,
                        help="<id>_merge_tax.csv (taxonkit lineage on merge.csv)")
    parser.add_argument("--merge-tax2", required=True, type=Path,
                        help="<id>_merge_tax2.csv (RPM, taxid, lineage)")
    parser.add_argument("--total-reads", required=True,
                        help="TR: number of R1 reads (log_p1)")
    parser.add_argument("--csv2-dir", required=True, type=Path)
    parser.add_argument("--frep-dir", required=True, type=Path)
    parser.add_argument("--brep-dir", required=True, type=Path)
    parser.add_argument("--report-path", required=True, type=Path)
    parser.add_argument("--bacteria-report-path", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    return parser.parse_args()


def _top_or_message(final_lines, message: str) -> str:
    return message if len(final_lines) == 0 else head3(final_lines)


def _report_text(sid, tr, pr2, prstd2, div, top, bac, vir, fun) -> str:
    return (
        f"[{sid}_REPORT]\n"
        f"Total_read(reads)\t{tr}\n"
        f"Pathogen_read(reads)\t{pr2}\n"
        f"Pathogen_read(RPM)\t{prstd2}\n"
        "\n>Diversity_Index\n"
        f"Simpson_Index(family)\t{div['simpson']['f']}\n"
        f"Simpson_Index(genus)\t{div['simpson']['g']}\n"
        f"Simpson_Index(species)\t{div['simpson']['s']}\n"
        f"\nShannon_Index(family)\t{div['shannon']['f']}\n"
        f"Shannon_Index(genus)\t{div['shannon']['g']}\n"
        f"Shannon_Index(species)\t{div['shannon']['s']}\n"
        f"\n>Top_3hit_pathogen(family) Name,reads,RPM,RA\n{top['f']}\n"
        f"\n>Top_3hit_pathogen(genus) Name,reads,RPM,RA\n{top['g']}\n"
        f"\n>Top_3hit_pathogen(species) Name,reads,RPM,RA\n{top['s']}\n"
        "\n>Top_3hit_Bacteria(species) Name,reads,RPM,RA\n"
        f"{_top_or_message(bac['s'][1], 'No Bacteria could be found.')}\n"
        "\n>Top_3hit_Virus(species) Name,reads,RPM,RA\n"
        f"{_top_or_message(vir[1], 'No Virus could be found.')}\n"
        "\n>Top_3hit_Fungi(species) Name,reads,RPM,RA\n"
        f"{_top_or_message(fun[1], 'No Fungi could be found.')}\n"
    )


def _bacteria_report_text(sid, tr, pr3, prstd3, div, top) -> str:
    return (
        f"[{sid}_BACTERIA_REPORT]\n"
        f"Total_read(reads)\t{tr}\n"
        f"Bacterial_read(reads)\t{pr3}\n"
        f"Bacterial_read(RPM)\t{prstd3}\n"
        "\n>Diversity_Index\n"
        f"Simpson_Index(family)\t{div['simpson']['f']}\n"
        f"Simpson_Index(genus)\t{div['simpson']['g']}\n"
        f"Simpson_Index(species)\t{div['simpson']['s']}\n"
        f"\nShannon_Index(family)\t{div['shannon']['f']}\n"
        f"Shannon_Index(genus)\t{div['shannon']['g']}\n"
        f"Shannon_Index(species)\t{div['shannon']['s']}\n"
        f"\n>Top_3hit_bacteria(family) Name,reads,RPM,RA\n{top['f']}\n"
        f"\n>Top_3hit_bacteria(genus) Name,reads,RPM,RA\n{top['g']}\n"
        f"\n>Top_3hit_bacteria(species) Name,reads,RPM,RA\n{top['s']}\n"
    )


def main() -> None:
    args = _parse_arguments()
    sid, threads = args.sample_id, args.threads
    for directory in (args.csv2_dir, args.frep_dir, args.brep_dir):
        directory.mkdir(parents=True, exist_ok=True)

    raw = read_lines(args.merge_tax)      # bash: ${OUT}_merge_tax.csv
    rpm = read_lines(args.merge_tax2)     # bash: ${OUT}_merge_tax2.csv

    # --- all organisms: K (without empty name filtering), then P, F, G, S
    total = {}
    for rank in "kpfgs":
        rp, mp, final = build_rank_table(raw, rpm, rank, threads, rank != "k")
        upper = rank.upper()
        write_lines(args.csv2_dir / f"{sid}_mf{upper}.csv", mp)
        write_lines(args.csv2_dir / f"{sid}_rf{upper}.csv", rp)
        write_lines(args.frep_dir / f"{sid}_f{upper}.csv", with_header(final))
        total[rank] = (mp, final)

    # --- subsets via grep (substring in any column, same as bash)
    def subset(lines, word):
        return [line for line in lines if word in line]

    # Bacterias: species, genus, family
    bacteria = {}
    raw_b, rpm_b = subset(raw, "Bacteria"), subset(rpm, "Bacteria")
    for rank in "sgf":
        rp, mp, final = build_rank_table(raw_b, rpm_b, rank, threads, True)
        upper = rank.upper()
        write_lines(args.csv2_dir / f"{sid}_mBacteria{upper}.csv", mp)
        write_lines(args.csv2_dir / f"{sid}_rBacteria{upper}.csv", rp)
        write_lines(args.brep_dir / f"{sid}_Bacteria_f{upper}.csv", with_header(final))
        bacteria[rank] = (mp, final)

    # Vírus and fungi: only species
    extras = {}
    for word, name in (("Virus", "Virus"), ("Fungi", "Fungi")):
        rp, mp, final = build_rank_table(
            subset(raw, word), subset(rpm, word), "s", threads, True
        )
        write_lines(args.csv2_dir / f"{sid}_m{name}.csv", mp)
        write_lines(args.csv2_dir / f"{sid}_r{name}.csv", rp)
        write_lines(args.frep_dir / f"{sid}_{name}.csv", with_header(final))
        extras[name] = (mp, final)

    # --- parameters
    pr2 = awk_sum_second_column(total["f"][1])        # sum of reads (fF1)
    pr3 = awk_sum_second_column(bacteria["f"][1])     # sum of reads (BacteriaF1)
    prstd2 = awk_sum_second_column(total["k"][0])     # sum of RPM (mfK)
    prstd3 = awk_sum_second_column(bacteria["f"][0])  # sum of RPM (mBacteriaF)

    div_total = {
        "simpson": {r: simpson_index(total[r][0]) for r in "fgs"},
        "shannon": {r: shannon_index(total[r][0]) for r in "fgs"},
    }
    top_total = {r: head3(total[r][1]) for r in "fgs"}
    div_bac = {
        "simpson": {r: simpson_index(bacteria[r][0]) for r in "fgs"},
        "shannon": {r: shannon_index(bacteria[r][0]) for r in "fgs"},
    }
    top_bac = {r: head3(bacteria[r][1]) for r in "fgs"}

    args.report_path.write_text(
        _report_text(sid, args.total_reads, pr2, prstd2, div_total, top_total,
                     bacteria, extras["Virus"], extras["Fungi"]),
        encoding="utf-8",
    )
    args.bacteria_report_path.write_text(
        _bacteria_report_text(sid, args.total_reads, pr3, prstd3, div_bac, top_bac),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
