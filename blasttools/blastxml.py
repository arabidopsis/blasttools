from __future__ import annotations

from pathlib import Path
import subprocess
from typing import Iterator
from collections.abc import Sequence

from dataclasses import dataclass, asdict
from uuid import uuid4
import pandas as pd
import click
from Bio.Blast.NCBIXML import parse  # type: ignore
from Bio.Blast.Record import Blast, Alignment, HSP  # type: ignore
from .blastapi import safe_which, remove_files, fasta_to_df, find_best, fetch_seq_df


class BlastXML:
    def runner(self, queryfasta: str, blastdb: str) -> Iterator[Blast]:
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
                    blastdb,
                    "-out",
                    out,
                ],
                check=False,
            )
            if r.returncode:
                raise click.ClickException("can't blast")
            with open(out, "rt", encoding="utf-8") as fp:
                yield from parse(fp)
        finally:
            remove_files([out])

    def run(self, queryfasta: str, blastdb: str) -> pd.DataFrame:
        return pd.DataFrame(
            [asdict(hit) for hit in hits(self.runner(queryfasta, blastdb))]
        )


def out5_to_df(xmlfile: str) -> pd.DataFrame:
    def run() -> Iterator[Blast]:
        with open(xmlfile, "rt", encoding="utf-8") as fp:
            yield from parse(fp)

    return pd.DataFrame([asdict(hit) for hit in hits(run())])


@dataclass
class Hit:
    query: str  # full string from fasta description line
    query_length: int
    accession: str
    accession_length: int  # accession length
    hsp_align_length: int
    hsp_bits: float
    hsp_score: float
    hsp_expect: float
    hsp_identities: int

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


def hits(xml: Iterator[Blast], full: bool = False) -> Iterator[Hit]:
    for b, a, h in unwind(xml):
        # b.query is the full line in the query fasta
        # actually <query-def>
        query = b.query.split(None, 1)[0] if not full else b.query
        yield Hit(
            query=query,
            query_length=b.query_length,
            accession=a.accession,  # saccver
            accession_length=a.length,
            hsp_align_length=h.align_length,  # alignment length
            hsp_bits=h.bits,  # bitscore
            hsp_score=h.score,  # bitscore?
            hsp_expect=h.expect,  # evalue
            hsp_identities=h.identities,
            hsp_match=h.match,
            hsp_query=h.query,
            hsp_query_start=h.query_start,  # qstart
            hsp_query_end=h.query_end,  # qend
            hsp_sbjct=h.sbjct,
            hsp_sbjct_start=h.sbjct_start,  # sstart
            hsp_sbjct_end=h.sbjct_end,  # send
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


def blastxml_to_df(queryfasta: str, blastdb: str) -> pd.DataFrame:
    bs = BlastXML()
    return pd.DataFrame([asdict(hit) for hit in hits(bs.runner(queryfasta, blastdb))])


def blastall(
    query: str, blastdbs: Sequence[str], best: int, with_seq: bool
) -> pd.DataFrame:
    df = fasta_to_df(query)
    if not df["id"].is_unique:
        raise click.ClickException(
            f'sequences IDs are not unique for query file "{query}"'
        )
    res = []

    b5 = BlastXML()
    for blastdb in blastdbs:
        rdf = b5.run(query, blastdb)

        print(rdf.columns, df.columns)

        if with_seq and "accession" in rdf.columns:
            saccver = list(rdf["accession"])
            sdf = fetch_seq_df(saccver, blastdb)
            sdf.rename(columns={"saccver": "accession"}, inplace=True)
            rdf = pd.merge(rdf, sdf, left_on="accession", right_on="accession")

        myrdf = find_best(
            rdf, df, nevalues=best, evalue_col="hsp_expect", query_col="query"
        )
        print(len(myrdf))
        myrdf["blastdb"] = Path(blastdb).name
        res.append(myrdf)
    ddf = pd.concat(res, axis=0)

    return ddf
