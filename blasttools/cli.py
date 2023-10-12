from __future__ import annotations
from collections.abc import Sequence
import click

from .blastapi import merge_fasta, buildall, blastall, save_df, toblastdb


@click.group()
@click.version_option()
def blast():
    pass


@blast.command(name="merge-fasta")
@click.argument("fasta1", type=click.Path(exists=True, dir_okay=False))
@click.argument("fasta2", type=click.Path(exists=True, dir_okay=False))
@click.argument("outfasta", type=click.Path(exists=False, dir_okay=False))
def merge_fasta_cmd(fasta1: str, fasta2: str, outfasta: str) -> None:
    """merge 2 fasta files based on sequence identity"""
    merge_fasta(fasta1, fasta2, outfasta)


@blast.command(name="build")
@click.option(
    "-b", "--build", "builddir", type=click.Path(file_okay=False, dir_okay=True)
)
@click.argument("fastas", nargs=-1)
def build_cmd(fastas: Sequence[str], builddir: str | None):
    """build blast databases from fasta files"""

    buildall(fastas, builddir=builddir)


@blast.command(name="blast")
@click.option("--out", help="output file", type=click.Path(dir_okay=False))
@click.option("--best", default=2, help="best (lowest) evalues [=0 take all]")
@click.option("--with_seq", is_flag=True, help="added sequence data to output")
@click.argument("query", type=click.Path(exists=True, dir_okay=False))
@click.argument("blastdbs", nargs=-1)
def blast_cmd(
    query: str,
    blastdbs: Sequence[str],
    best: int,
    with_seq: bool,
    out: str | None,
):
    """blast databases"""

    df = blastall(query, toblastdb(blastdbs), with_seq=with_seq, best=best)
    if out is None:
        out = query + ".csv"

    click.secho(f"writing {out}", fg="green")
    save_df(df, out, index=False)
