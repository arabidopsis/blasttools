from __future__ import annotations

import re
import subprocess
from collections.abc import Sequence
from dataclasses import asdict
from dataclasses import dataclass
from dataclasses import fields
from pathlib import Path
from typing import Any
from typing import Iterator
from uuid import uuid4

import click
import pandas as pd
from Bio.Blast.NCBIXML import parse  # type: ignore
from Bio.Blast.Record import Alignment  # type: ignore
from Bio.Blast.Record import Blast
from Bio.Blast.Record import HSP

from .blastapi import Blast6
from .blastapi import remove_files


class BlastXML(Blast6):
    def __init__(
        self,
        header: Sequence[str] | None = None,
        num_threads: int = 1,
        blastp: bool = True,
    ):
        if header is None:
            header = HEADER
        super().__init__(header, num_threads, blastp)
        self._h = set(header)

    def runner(self, queryfasta: str | Path, blastdb: str | Path) -> Iterator[Blast]:
        outfmt = "5"
        key = uuid4()
        out = f"{key}.xml"
        blast = self.get_blast()
        try:
            r = subprocess.run(
                [
                    blast,
                    "-outfmt",
                    outfmt,
                    "-query",
                    queryfasta,
                    "-db",
                    blastdb,
                    "-out",
                    out,
                    "-num_threads",
                    str(self.num_threads),
                ],
                check=False,
            )
            if r.returncode:
                raise click.ClickException(f"Can't blast {queryfasta}")
            with open(out, encoding="utf-8") as fp:
                yield from parse(fp)
        finally:
            remove_files([out])

    def todict(self, hit: Hit) -> dict[str, Any]:
        return {k: v for k, v in asdict(hit).items() if k in self._h}

    def run(self, queryfasta: str | Path, blastdb: str | Path) -> pd.DataFrame:
        todict = self.todict
        return pd.DataFrame(
            [todict(hit) for hit in hits(self.runner(queryfasta, blastdb))],
        )


def out5_to_df(xmlfile: str) -> pd.DataFrame:
    def run() -> Iterator[Blast]:
        with open(xmlfile, encoding="utf-8") as fp:
            yield from parse(fp)

    return pd.DataFrame([asdict(hit) for hit in hits(run())])


GAPS = re.compile("[-]+").finditer


@dataclass
class Hit:
    qaccver: str
    qlen: int
    saccver: str
    slen: int
    length: int
    bitscore: float
    score: float
    evalue: float
    nident: int
    positive: int
    gaps: int
    match: str  # not in 6
    query: str  # not in 6
    qstart: int
    qend: int
    sbjct: str  # not in 6
    sstart: int
    send: int
    gapopen: int
    mismatch: int
    pident: float


HEADER = [f.name for f in fields(Hit)]


def unwind(xml: Iterator[Blast]) -> Iterator[tuple[Blast, Alignment, HSP]]:
    b: Blast
    a: Alignment
    h: HSP
    for b in xml:
        for a in b.alignments:
            for h in a.hsps:
                yield b, a, h


def mismatch(hsp: HSP) -> int:
    return hsp.align_length - hsp.identities - hsp.gaps


def pident(hsp: HSP) -> float:
    return hsp.identities * 100.0 / hsp.align_length


def gapopen(hsp: HSP) -> int:
    # sbjct.str.count("[-]+") + query.str.count("[-]+")
    return sum(1 for _ in GAPS(hsp.sbjct)) + sum(1 for _ in GAPS(hsp.query))


def hits(xml: Iterator[Blast], full: bool = False) -> Iterator[Hit]:
    for b, a, h in unwind(xml):
        # b.query is the full line in the query fasta
        # actually <query-def>
        queryid = b.query.split(None, 1)[0] if not full else b.query
        yield Hit(
            qaccver=queryid,
            qlen=b.query_length,
            saccver=a.accession,
            slen=a.length,
            length=h.align_length,  # alignment length
            bitscore=h.bits,
            score=h.score,
            evalue=h.expect,
            nident=h.identities,
            gaps=h.gaps,
            positive=h.positives,
            match=h.match,  # not in --format=6
            query=h.query,  # not in --format=6
            qstart=h.query_start,
            qend=h.query_end,
            sbjct=h.sbjct,  # not in --format=6
            sstart=h.sbjct_start,
            send=h.sbjct_end,
            gapopen=gapopen(h),
            mismatch=mismatch(h),
            pident=pident(h),
        )


# this will work with df.apply(hsp_match, axis=1)
def hsp_match(hsp: Hit, width: int = 50, right: int = 0) -> str:
    lines = [
        f"Score {hsp.score:.0f} ({hsp.bitscore:.0f} bits), expectation {hsp.evalue:.1e},"
        f" alignment length {hsp.length}",
    ]
    if width <= 0:
        width = hsp.length
    if hsp.length <= width:
        lines.append(
            f"Query:{str(hsp.qstart).rjust(8)} {hsp.query} {hsp.qend}",
        )
        lines.append(f"               {hsp.match}")
        lines.append(
            f"Sbjct:{str(hsp.sstart).rjust(8)} {hsp.sbjct} {hsp.send}",
        )
    elif right <= 0:
        query_end = hsp.qstart
        sbjct_end = hsp.sstart
        for q in range(0, hsp.length, width):
            query = hsp.query[q : q + width]
            sbjct = hsp.sbjct[q : q + width]

            s = " " * (width - len(query))

            query_start = query_end
            sbjct_start = sbjct_end
            query_end += len(query) - query.count("-")
            sbjct_end += len(sbjct) - sbjct.count("-")
            lines.append(
                f"Query:{str(query_start).rjust(8)} {query}{s} {query_end - 1}",
            )
            lines.append(f"{' '*15}{hsp.match[q:q+width]}")
            lines.append(
                f"Sbjct:{str(sbjct_start).rjust(8)} {sbjct}{s} {sbjct_end - 1}",
            )
            lines.append("")
        del lines[-1]
    else:
        left = width - right - 3 + 1

        lines.append(
            f"Query:{str(hsp.qstart).rjust(8)} {hsp.query[:left]}...{hsp.query[-right:]} {hsp.qend}",
        )
        lines.append(f"               {hsp.match[:left]}...{hsp.match[-right:]}")
        lines.append(
            f"Sbjct:{str(hsp.sstart).rjust(8)} {hsp.sbjct[:left]}...{hsp.sbjct[-right:]} {hsp.send}",
        )
    return "\n".join(lines)


def blastxml_to_df(
    queryfasta: str | Path,
    blastdb: str | Path,
    num_threads: int = 1,
    blastp: bool = True,
) -> pd.DataFrame:
    bs = BlastXML(HEADER, num_threads=num_threads, blastp=blastp)
    return bs.run(queryfasta, blastdb)
