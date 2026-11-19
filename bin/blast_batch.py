#!/usr/bin/env python3
import argparse
import re
import sys
from contextlib import ExitStack
from pathlib import Path
from typing import Dict, List, Tuple

SEPARATOR = "__"
DESCRIPTION = 'Run one blastn for many FASTA chunks and split the result back per chunk.'
_QUERY_PREFIX = "# Query: "
_BLOCK_START = re.compile(r"^# [A-Z]*BLAST[A-Z]* \d")
_FOOTER = re.compile(r"^# BLAST processed \d+ queries")
_ENC = {"encoding": "utf-8", "errors": "surrogateescape"}


def _chunk_name(path: str) -> str:
    name = Path(path).name
    return name[: -len(".fasta")] if name.endswith(".fasta") else name


def _prefix(index: int) -> str:
    return f"c{index:06d}{SEPARATOR}"


def concat(chunks: List[str], fasta_out: Path, manifest_out: Path) -> int:
    reads = 0
    with open(fasta_out, "w", **_ENC) as fasta, open(manifest_out, "w", **_ENC) as manifest:
        for index, chunk in enumerate(chunks, start=1):
            prefix = _prefix(index)
            manifest.write(f"{prefix}\t{_chunk_name(chunk)}\n")
            with open(chunk, "r", **_ENC) as handle:
                for line in handle:
                    if line.startswith(">"):
                        reads += 1
                        fasta.write(">" + prefix + line[1:])
                    else:
                        fasta.write(line)
                    if not line.endswith("\n"):
                        fasta.write("\n")
    return reads


def read_manifest(path: Path) -> Tuple[Dict[str, str], List[str]]:
    chunk_by_prefix: Dict[str, str] = {}
    order: List[str] = []
    with open(path, "r", **_ENC) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            prefix, chunk = line.split("\t", 1)
            if prefix in chunk_by_prefix:
                raise ValueError(f"duplicated prefix in manifest: {prefix}")
            chunk_by_prefix[prefix] = chunk
            order.append(prefix)
    return chunk_by_prefix, order


def split(blast_path: Path, manifest_path: Path, outdir: Path) -> Dict[str, int]:
    chunk_by_prefix, order = read_manifest(manifest_path)
    outdir.mkdir(parents=True, exist_ok=True)
    counts = {chunk_by_prefix[p]: 0 for p in order}
    footer_template = None

    with ExitStack() as stack:
        outputs = {
            prefix: stack.enter_context(
                open(outdir / f"{chunk_by_prefix[prefix]}_blast.txt", "w", **_ENC)
            )
            for prefix in order
        }
        if not blast_path.exists() or blast_path.stat().st_size == 0:
            return counts

        pending: List[str] = []
        current = None
        with open(blast_path, "r", **_ENC) as handle:
            for line in handle:
                if _FOOTER.match(line):
                    footer_template = line
                    current = None
                    continue
                if _BLOCK_START.match(line):
                    pending, current = [line], None
                    continue
                if line.startswith(_QUERY_PREFIX):
                    query = line[len(_QUERY_PREFIX):]
                    cut = query.find(SEPARATOR)
                    prefix = query[: cut + len(SEPARATOR)] if cut >= 0 else ""
                    if prefix not in outputs:
                        raise ValueError(f"query without a known prefix: {line.rstrip()}")
                    current = prefix
                    counts[chunk_by_prefix[prefix]] += 1
                    out = outputs[prefix]
                    out.writelines(pending)
                    pending = []
                    out.write(_QUERY_PREFIX + query[len(prefix):])
                    continue
                if current is None:
                    pending.append(line)
                    continue
                if line.startswith("#"):
                    outputs[current].write(line)
                elif line.startswith(current):
                    outputs[current].write(line[len(current):])
                elif line.strip():
                    raise ValueError(f"hit line does not start with {current}: {line.rstrip()}")
                else:
                    outputs[current].write(line)

        if footer_template is not None:
            for prefix in order:
                count = counts[chunk_by_prefix[prefix]]
                if count:
                    outputs[prefix].write(re.sub(r"\d+", str(count), footer_template, count=1))
    return counts


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    sub = parser.add_subparsers(dest="command", required=True)

    p_concat = sub.add_parser("concat", help="concatenate chunks with prefixed read IDs")
    p_concat.add_argument("--fasta", required=True, type=Path)
    p_concat.add_argument("--manifest", required=True, type=Path)
    p_concat.add_argument("chunks", nargs="+")

    p_split = sub.add_parser("split", help="split the BLAST output back per chunk")
    p_split.add_argument("--blast", required=True, type=Path)
    p_split.add_argument("--manifest", required=True, type=Path)
    p_split.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = _parse_arguments()
    if args.command == "concat":
        reads = concat(args.chunks, args.fasta, args.manifest)
        print(f"concat: {len(args.chunks)} chunks, {reads} reads", file=sys.stderr)
    else:
        counts = split(args.blast, args.manifest, args.outdir)
        print(f"split: {len(counts)} chunks, {sum(counts.values())} queries", file=sys.stderr)


if __name__ == "__main__":
    main()
