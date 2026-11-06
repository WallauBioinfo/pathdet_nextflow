#!/usr/bin/env python3
"""Índices de diversidade e blocos de relatório com a formatação do bash original."""

import math
from typing import List, Sequence

from awk_compat import awk_float, cut_field, perl_num


def _relative_abundances(mp_lines: Sequence[str]) -> List[float]:
    """cut -f3 (coluna RA do `mp`), numerificada como o perl."""
    return [awk_float(cut_field(line, 3)) for line in mp_lines]


def simpson_index(mp_lines: Sequence[str]) -> str:
    """perl -nle '$D += $_*$_; print 1-$D if(eof);' (sem linhas não imprime nada)."""
    values = _relative_abundances(mp_lines)
    if not values:
        return ""
    return perl_num(1 - sum_in_order(v * v for v in values))


def shannon_index(mp_lines: Sequence[str]) -> str:
    """perl -nle '$H += $_*log($_); print $H*(-1) if(eof);' (log(0) aborta o perl)."""
    values = _relative_abundances(mp_lines)
    if not values or any(v <= 0 for v in values):
        return ""
    return perl_num(-sum_in_order(v * math.log(v) for v in values))


def sum_in_order(values) -> float:
    total = 0.0
    for value in values:
        total += value
    return total


def head3(lines: Sequence[str]) -> str:
    """head -n3 guardado numa variável do shell (sem '\\n' final)."""
    return "\n".join(lines[:3])
