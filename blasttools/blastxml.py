from __future__ import annotations

import subprocess
from typing import Iterator

from dataclasses import dataclass
from uuid import uuid4
from Bio.Blast.NCBIXML import parse
from Bio.Blast.Record import Blast, Alignment, HSP
from .blast import safe_which, blast_dir, remove_files


def doblastxml(queryfasta: str, plant: str, release: int) -> Iterator[Blast]:
    blastdir = blast_dir(release)
    outfmt = "5"
    out = f"{uuid4()}.xml"
    blastp = safe_which("blastp")
    try:
        r = subprocess.run(
            [
                blastp,
                "-outfmt",
                outfmt,
                "-query",
                queryfasta,
                "-db",
                blastdir / plant,
                "-out",
                out,
            ],
            check=False,
        )
        if r.returncode:
            raise RuntimeError("can't blast")
        with open(out, "rt", encoding="utf-8") as fp:
            yield from parse(fp)
    finally:
        remove_files([out])




@dataclass
class Hit:
    query: str
    query_length: int
    accession: str
    accession_length: int  # accession length
    hsp_align_length: int
    hsp_bits: float
    hsp_expect: float
    hsp_identities: int
    hsp_score: float

    hsp_match: str
    hsp_query: str
    hsp_query_end: int
    hsp_query_start: int
    hsp_sbjct: str
    hsp_sbjct_end: int
    hsp_sbjct_start: int


def unwind(xml: Iterator[Blast]) -> Iterator[tuple[Blast, Alignment, HSP]]:
    b: Blast
    a: Alignment
    h: HSP
    for b in xml:
        for a in b.alignments:
            for h in a.hsps:
                yield b, a, h


def hits(xml: Iterator[Blast]) -> Iterator[Hit]:
    for b, a, h in unwind(xml):
        yield Hit(
            query=b.query,
            query_length=b.query_length,
            accession=a.accession,
            accession_length=a.length,
            hsp_align_length=h.align_length,
            hsp_bits=h.bits,
            hsp_expect=h.expect,
            hsp_identities=h.identities,
            hsp_score=h.score,
            hsp_match=h.match,
            hsp_query=h.query,
            hsp_query_start=h.query_start,
            hsp_query_end=h.query_end,
            hsp_sbjct=h.sbjct,
            hsp_sbjct_start=h.sbjct_start,
            hsp_sbjct_end=h.sbjct_end,
        )


def hsp_match(hsp: HSP, width: int = 50, right: int = 0) -> str:
    lines = [
        f"Score {hsp.score:.0f} ({hsp.bits:.0f} bits), expectation {hsp.expect:.1e},"
        f" alignment length {hsp.align_length}"
    ]
    if width <= 0:
        width = hsp.align_length
    if hsp.align_length <= width:
        lines.append(
            f"Query:{str(hsp.query_start).rjust(8)} {hsp.query} {hsp.query_end}"
        )
        lines.append(f"               {hsp.match}")
        lines.append(
            f"Sbjct:{str(hsp.sbjct_start).rjust(8)} {hsp.sbjct} {hsp.sbjct_end}"
        )
    elif right <= 0:
        query_end = hsp.query_start
        sbjct_end = hsp.sbjct_start
        for q in range(0, hsp.align_length, width):
            query = hsp.query[q : q + width]
            sbjct = hsp.sbjct[q : q + width]

            s = " " * (width - len(query))

            query_start = query_end
            sbjct_start = sbjct_end
            query_end += len(query) - query.count("-")
            sbjct_end += len(sbjct) - sbjct.count("-")
            lines.append(
                f"Query:{str(query_start).rjust(8)} {query}{s} {query_end - 1}"
            )
            lines.append(f"{' '*15}{hsp.match[q:q+width]}")
            lines.append(
                f"Sbjct:{str(sbjct_start).rjust(8)} {sbjct}{s} {sbjct_end - 1}"
            )
            lines.append("")
        del lines[-1]
    else:
        left = width - right - 3 + 1

        lines.append(
            f"Query:{str(hsp.query_start).rjust(8)} {hsp.query[:left]}...{hsp.query[-right:]} {hsp.query_end}"
        )
        lines.append(f"               {hsp.match[:left]}...{hsp.match[-right:]}")
        lines.append(
            f"Sbjct:{str(hsp.sbjct_start).rjust(8)} {hsp.sbjct[:left]}...{hsp.sbjct[-right:]} {hsp.sbjct_end}"
        )
    return "\n".join(lines)
