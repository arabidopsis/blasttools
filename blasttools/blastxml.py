from __future__ import annotations

import re
import subprocess
import sys
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
from .blastapi import get_cat
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

    def get_output(self) -> str:
        return f"{uuid4()}.xml"

    def runner(
        self,
        queryfasta: str | Path,
        blastdb: str | Path,
        *,
        needs_translation: bool = False,
    ) -> Iterator[Blast]:
        outfmt = "5"  # also 16...
        out = self.get_output()
        blast = self.get_blast()

        queryfasta = Path(queryfasta)

        if not needs_translation:
            cat_cmd = get_cat(queryfasta)
            cat_cmd.append(queryfasta.name)
        else:
            cat_cmd = [sys.executable, "-m", "blasttools.translate", queryfasta.name]
        try:
            with subprocess.Popen(
                cat_cmd,
                stdout=subprocess.PIPE,
                cwd=str(queryfasta.parent),
            ) as p1:
                with subprocess.Popen(
                    [
                        blast,
                        "-outfmt",
                        outfmt,
                        "-query",
                        "-",
                        "-db",
                        str(blastdb),
                        "-out",
                        out,
                        "-num_threads",
                        str(self.num_threads),
                        *self.get_blast_args(),
                    ],
                    stdin=p1.stdout,
                ) as p2:
                    if p1.stdout:
                        p1.stdout.close()
                    p2.wait()
                    p1.wait()

                    if p2.returncode:
                        b = "blastp" if self.blastp else "blastn"
                        raise click.ClickException(f"Can't run {b} using {queryfasta}")
                with open(out, encoding="utf-8") as fp:
                    yield from parse(fp)
        finally:
            remove_files([out])

    def todict(self, hit: Hit) -> dict[str, Any]:
        return {k: v for k, v in asdict(hit).items() if k in self._h}

    def run(
        self,
        queryfasta: str | Path,
        blastdb: str | Path,
        *,
        needs_translation: bool = False,
    ) -> pd.DataFrame:
        todict = self.todict
        return pd.DataFrame(
            [
                todict(hit)
                for hit in hits(
                    self.runner(
                        queryfasta,
                        blastdb,
                        needs_translation=needs_translation,
                    ),
                )
            ],
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


def gaps(hsp: HSP) -> int:
    # same as hsp.gaps
    return hsp.sbjct.count("-") + hsp.query.count("-")


def hits(xml: Iterator[Blast], full: bool = False) -> Iterator[Hit]:
    for b, a, h in unwind(xml):
        # b.query is the full line in the query fasta
        # actually <query-def>
        queryid = b.query.split(None, maxsplit=1)[0] if not full else b.query
        # assert h.gaps == gaps(h)
        # a.hit_id == sseqid == a.accession
        # undocumented? (see set_hit_accession in Bio/Blast/NCBIXML.py)
        # assert a.hit_id == a.accession, (a.hit_id, a.accession)
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
