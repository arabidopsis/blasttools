from __future__ import annotations

from itertools import chain
from itertools import islice
from os.path import splitext
from pathlib import Path
from typing import Iterator
from typing import TypeVar

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

    fname = f"{fname}-{{num}}.fa.gz" if fmt is None else fmt

    ret = []
    if target_dir is None:
        target_dir = fastafile.parent
    else:
        target_dir = Path(target_dir)
        target_dir.mkdir(parents=True, exist_ok=True)
    for num, it in enumerate(batched(read_fasta(fastafile), batch), start=1):
        sname = target_dir / fname.format(num=num)
        ret.append(sname)
        write_fasta_iter(sname, it)

    return ret
