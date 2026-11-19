#!/usr/bin/env python3
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

_NUM_PREFIX = re.compile(r"[ \t\n]*[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")
_FIELD = re.compile(r"[^ \t\n]+")
_GNU_NUM = re.compile(r"[ \t]*(-?)(\d*)(?:\.(\d*))?")
_ENC = {"encoding": "utf-8", "errors": "surrogateescape"}


def read_lines(path) -> List[str]:
    p = Path(path)
    if not p.exists():
        return []
    text = p.read_text(**_ENC)
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    return lines


def write_lines(path, lines: Sequence[str]) -> None:
    Path(path).write_text("".join(f"{line}\n" for line in lines), **_ENC)


def count_newlines(path) -> int:
    p = Path(path)
    if not p.exists():
        return 0
    total = 0
    with open(p, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            total += chunk.count(b"\n")
    return total


def awk_float(text: str) -> float:
    match = _NUM_PREFIX.match(text)
    return float(match.group(0)) if match else 0.0


def awk_num(value: float) -> str:
    if math.isfinite(value) and value == int(value) and abs(value) < 1e16:
        return str(int(value))
    return "%.6g" % value


def perl_num(value: float) -> str:
    text = "%.15g" % value
    return "0" if text == "-0" else text


def awk_fields(line: str) -> List[str]:
    return _FIELD.findall(line)


def awk_sum_by_key(lines: Iterable[str]) -> Dict[str, float]:
    acc: Dict[str, float] = {}
    for line in lines:
        fields = awk_fields(line)
        key = fields[1] if len(fields) > 1 else ""
        acc[key] = acc.get(key, 0.0) + (awk_float(fields[0]) if fields else 0.0)
    return acc


def awk_sum_second_column(lines: Sequence[str]) -> str:
    if not lines:
        return ""
    total = 0.0
    for line in lines:
        fields = awk_fields(line)
        total += awk_float(fields[1]) if len(fields) > 1 else 0.0
    return awk_num(total)


def awk_greater_than_one(text: str) -> bool:
    stripped = text.strip(" \t\n")
    if stripped and _NUM_PREFIX.fullmatch(stripped):
        return float(stripped) > 1
    return text > "1"


def _skip_fields(line: str, count: int) -> int:
    pos, size = 0, len(line)
    for _ in range(count):
        while pos < size and line[pos] in " \t":
            pos += 1
        while pos < size and line[pos] not in " \t":
            pos += 1
    return pos


def gnu_numeric_key(line: str, field: int) -> float:
    match = _GNU_NUM.match(line, _skip_fields(line, field - 1))
    if not match or not (match.group(2) or match.group(3)):
        return 0.0
    return float(f"{match.group(1)}{match.group(2) or '0'}.{match.group(3) or '0'}")


def sort_numeric_desc(lines: Iterable[str], field: int) -> List[str]:
    return sorted(
        lines,
        key=lambda line: (-gnu_numeric_key(line, field), line.encode(**_ENC)),
    )


def drop_leading_whitespace(lines: Iterable[str]) -> List[str]:
    return [line for line in lines if not (line and line[0] in " \t\n\r\f\v")]


def cut_fields(line: str, fields: Sequence[int]) -> str:
    if "\t" not in line:
        return line
    cols = line.split("\t")
    return "\t".join(cols[i - 1] for i in fields if i - 1 < len(cols))


def cut_field(line: str, field: int) -> str:
    return cut_fields(line, [field])


def cut_from_fourth(line: str) -> str:
    if "\t" not in line:
        return line
    cols = line.split("\t")
    return "\t".join([cols[0]] + cols[3:])


def cut_from_third(line: str) -> str:
    if "\t" not in line:
        return line
    cols = line.split("\t")
    return "\t".join([cols[0]] + cols[2:])


def paste_cut(left: Sequence[str], right: Sequence[str]) -> List[str]:
    out = []
    for i in range(max(len(left), len(right))):
        a = left[i] if i < len(left) else ""
        b = right[i] if i < len(right) else ""
        out.append(cut_fields(f"{a}\t{b}", (1, 2, 4, 5)))
    return out


def add_fraction_column(lines: Sequence[str]) -> List[str]:
    if not lines:
        return []
    total = 0.0
    for line in lines:
        fields = awk_fields(line)
        total += awk_float(fields[1]) if len(fields) > 1 else 0.0
    if total == 0:
        return []
    out = []
    for line in lines:
        fields = awk_fields(line)
        value = awk_float(fields[1]) if len(fields) > 1 else 0.0
        out.append(f"{line}\t{awk_num(value / total)}")
    return out


def bc_percent(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    quotient = (numerator * 10000) // denominator
    scaled = quotient * 100
    if scaled == 0:
        return "0"
    integer, fraction = divmod(scaled, 10000)
    return f"{integer}.{fraction:04d}" if integer else f".{fraction:04d}"


def bc_rpm(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return ""
    return str((numerator * 1000000) // denominator)
