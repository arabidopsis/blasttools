from __future__ import annotations

from itertools import chain
from itertools import islice
from os.path import splitext
from pathlib import Path
from typing import Iterator
from typing import TypeVar

import click
from Bio.SeqRecord import SeqRecord

from .blastapi import read_fasta
from .blastapi import write_fasta_iter

sentinel = object()

T = TypeVar("T")


def batched(iterable: Iterator[T], n: int) -> Iterator[Iterator[T]]:
    it = iter(iterable)
    while True:
        batch = islice(it, n)
        b = next(batch, sentinel)
        if b is sentinel:
            break
        yield chain([b], batch)  # type: ignore


def split_fasta(
    fastafile: str | Path,
    batch: int = 10000,
    target_dir: Path | str | None = None,
    fmt: str | None = None,
) -> list[Path]:
    fastafile = Path(fastafile)
    fname, _ = splitext(fastafile.name)

    fname = f"{fname}-{{num:03d}}.fa.gz" if fmt is None else fmt

    ret = []
    if target_dir is None:
        target_dir = fastafile.parent
    else:
        click.secho(f'creating directory: "{target_dir}"', fg="yellow", err=True)
        target_dir = Path(target_dir)
        target_dir.mkdir(parents=True, exist_ok=True)
    for num, it in enumerate(batched(read_fasta(fastafile), batch), start=1):
        sname = target_dir / fname.format(num=num)
        ret.append(sname)
        write_fasta_iter(sname, it)

    return ret


def fasta_filter(fastafile: str | Path, ids: set[str]) -> Iterator[SeqRecord]:
    for rec in read_fasta(fastafile):
        if rec.id in ids:
            ids.remove(rec.id)
            yield rec


def fasta_filter_out(fastafile: str | Path, outfile: str | Path, ids: set[str]) -> None:
    write_fasta_iter(Path(outfile), fasta_filter(fastafile, ids))


def is_rna(s: str) -> bool:
    return set(s).issubset({"A", "C", "G", "U", "T"})


def translate(fasta: str | Path) -> Iterator[SeqRecord]:
    for rec in read_fasta(fasta):
        if rec.seq is not None and is_rna(str(rec.seq)):
            rec.seq = rec.seq.transcribe().translate()
            # print('translated', rec.id,rec.seq, file=sys.stderr)
        yield rec
