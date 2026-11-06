#!/usr/bin/env python3
"""Processamento taxonómico idêntico ao do bash original.

Em vez do ete3 (banco próprio, taxdump possivelmente diferente) chamamos o
`taxonkit` — o mesmo binário e o mesmo taxdump ($TAXONKIT_DB) do bash — e
replicamos em Python as pipelines `awk | sort | paste | cut`.
"""

import subprocess
from typing import List, Sequence, Tuple

from awk_compat import (
    add_fraction_column,
    awk_float,
    awk_fields,
    awk_greater_than_one,
    awk_num,
    awk_sum_by_key,
    cut_fields,
    cut_from_fourth,
    cut_from_third,
    drop_leading_whitespace,
    paste_cut,
    sort_numeric_desc,
)

_ENC = {"encoding": "utf-8", "errors": "surrogateescape"}


# ----------------------------------------------------------------------- taxonkit
def _run_taxonkit(args: Sequence[str], lines: Sequence[str]) -> List[str]:
    if not lines:
        return []
    payload = "".join(f"{line}\n" for line in lines).encode(**_ENC)
    result = subprocess.run(
        ["taxonkit", *args], input=payload, capture_output=True, check=False
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"taxonkit {' '.join(args)} falhou: {result.stderr.decode(**_ENC)}"
        )
    out = result.stdout.decode(**_ENC).split("\n")
    if out and out[-1] == "":
        out.pop()
    return out


def taxonkit_lineage(lines: Sequence[str], taxid_column: int = 2) -> List[str]:
    """taxonkit lineage -i N (acrescenta a linhagem como última coluna)."""
    return _run_taxonkit(["lineage", "-i", str(taxid_column)], lines)


def taxonkit_reformat(
    lines: Sequence[str], lineage_column: int, rank: str, threads: int
) -> List[str]:
    """taxonkit reformat -j T -i N -f {rank}."""
    return _run_taxonkit(
        ["reformat", "-j", str(threads), "-i", str(lineage_column), "-f", "{%s}" % rank],
        lines,
    )


# -------------------------------------------------- etapas comuns (Rapid e BLAST)
def build_rank_table(
    raw_lines: Sequence[str],
    rpm_lines: Sequence[str],
    rank: str,
    threads: int,
    drop_empty_names: bool,
) -> Tuple[List[str], List[str], List[str]]:
    """Devolve (rp, mp, final) de um rank, como no bash.

    rp    = reads por nome        (`..._rp<RANK>.csv` / `..._rf<RANK>.csv`)
    mp    = RPM por nome + RA     (`..._mp<RANK>.csv` / `..._mf<RANK>.csv`)
    final = paste rp mp | cut -f 1,2,4,5   (corpo de `_p<RANK>.csv` / `_f<RANK>.csv`)

    `drop_empty_names=False` só no kingdom (o bash não aplica o grep -v '^\\s' aí).
    """
    mp = _sorted_rank_lines(rpm_lines, rank, threads, drop_empty_names)
    mp = add_fraction_column(mp)
    rp = _sorted_rank_lines(raw_lines, rank, threads, drop_empty_names)
    return rp, mp, paste_cut(rp, mp)


def _sorted_rank_lines(
    lines: Sequence[str], rank: str, threads: int, drop_empty_names: bool
) -> List[str]:
    reformatted = taxonkit_reformat(lines, 3, rank, threads)
    names = [cut_from_fourth(line).replace(" ", "_") for line in reformatted]
    totals = awk_sum_by_key(names)
    rows = [f"{name}\t{awk_num(total)}" for name, total in totals.items()]
    rows = sort_numeric_desc(rows, 2)
    return drop_leading_whitespace(rows) if drop_empty_names else rows


def with_header(body: Sequence[str]) -> List[str]:
    return ["Microorganisms,reads,RPM,RA", *body]


# ----------------------------------------------------------------- só o Rapid
def kreport_to_p2(kreport_lines: Sequence[str]) -> List[str]:
    """grep -v 'unclassified' | cut -f 3,5 | awk -F, '$1 > 1'."""
    out = []
    for line in kreport_lines:
        if "unclassified" in line:
            continue
        selected = cut_fields(line, (3, 5))
        if awk_greater_than_one(selected.split(",")[0]):
            out.append(selected)
    return out


def p2_to_tax_lines(p2_lines: Sequence[str]) -> List[str]:
    """awk soma por taxid | sort -k 1nr,1  ->  'reads<TAB>taxid'."""
    totals = awk_sum_by_key(p2_lines)
    rows = [f"{awk_num(total)}\t{taxid}" for taxid, total in totals.items()]
    return sort_numeric_desc(rows, 1)


def tax_to_rpm_lines(tax_lines: Sequence[str], total_reads: int) -> List[str]:
    """awk '{print $1/r*1000000"\\t"$0}' | cut -f 1,3-  ->  'RPM<TAB>taxid<TAB>linhagem'."""
    out = []
    for line in tax_lines:
        fields = awk_fields(line)
        reads = awk_float(fields[0]) if fields else 0.0
        rpm = awk_num(reads / total_reads * 1000000)
        out.append(cut_from_third(f"{rpm}\t{line}"))
    return out
