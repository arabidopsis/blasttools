from __future__ import annotations
from collections.abc import Sequence
import click

from .blastapi import merge_fasta, buildall, blastall, save_df, toblastdb


@click.group(
    epilog=click.style("Commands to run blasts\n", fg="magenta"),
)
@click.version_option()
def blast() -> None:
    pass


@blast.command()
def update() -> None:
    """Update this package"""
    import subprocess
    import sys

    from .config import REPO

    ret = subprocess.call([sys.executable, "-m", "pip", "install", "-U", REPO])
    if ret:
        click.secho(f"can't install {REPO}", fg="red", err=True)
        raise click.Abort()


@blast.command(name="fasta-merge")
@click.argument("fasta1", type=click.Path(exists=True, dir_okay=False))
@click.argument("fasta2", type=click.Path(exists=True, dir_okay=False))
@click.argument("outfasta", type=click.Path(exists=False, dir_okay=False))
def merge_fasta_cmd(fasta1: str, fasta2: str, outfasta: str) -> None:
    """merge 2 fasta files based on sequence identity"""
    merge_fasta(fasta1, fasta2, outfasta)


@blast.command(name="build")
@click.option(
    "-b",
    "--build",
    "builddir",
    type=click.Path(file_okay=False, dir_okay=True),
    help="build all databases in this directory",
)
@click.option("-m", "--merge", help="merge all fastafiles into one")
@click.argument(
    "fastas", nargs=-1, type=click.Path(exists=True, file_okay=True, dir_okay=False)
)
def build_cmd(fastas: Sequence[str], builddir: str | None, merge: str | None) -> None:
    """build blast databases from fasta files"""
    buildall(fastas, builddir=builddir, merge=merge)


@blast.command(name="blast")
@click.option("--out", help="output filename (default is to write <query>.csv)", type=click.Path(dir_okay=False))
@click.option("--best", default=0, help="best (lowest) evalues [=0 take all]")
@click.option("--xml", is_flag=True, help="run with xml output")
@click.option("--with-seq", is_flag=True, help="add sequence data to output")
@click.option(
    "-c",
    "--columns",
    help="space separated list of columns (see columns cmd for a list of valid columns)",
)
@click.argument("query", type=click.Path(exists=True, dir_okay=False))
@click.argument("blastdbs", nargs=-1)
def blast_cmd(
    query: str,
    blastdbs: Sequence[str],
    best: int,
    with_seq: bool,
    out: str | None,
    xml: bool,
    columns: str | None,
) -> None:
    """blast databases"""
    from .blastapi import mkheader, has_pdatabase
    from .blastxml import blastall as blastall5

    if len(blastdbs) == 0:
        return
    blastdbs = toblastdb(blastdbs)
    missing = {b for b in blastdbs if not has_pdatabase(b)}
    if missing:
        m = ", ".join(missing)
        raise click.BadParameter(f"missing databases {m}", param_hint="blastdbs")

    myheader = None
    if columns is not None:
        if xml:
            raise click.BadParameter("can't have header with xml", param_hint="xml")
        myheader = mkheader(columns)

    if xml:
        df = blastall5(query, blastdbs, with_seq=with_seq, best=best)
    else:
        df = blastall(query, blastdbs, with_seq=with_seq, best=best, header=myheader)
    if out is None:
        out = query + ".csv"

    click.secho(f"writing {out}", fg="green")
    save_df(df, out, index=False)


@blast.command(name="columns")
def columns_cmd() -> None:
    """Show possible output columns for blast"""
    from .columns import VALID
    from .blastapi import HEADER

    mx = len(max(VALID, key=len))

    for k, v in VALID.items():
        s = "*" if k in HEADER else " "
        click.echo(f"{k:<{mx}}{s}: {v}")
    click.echo()
    click.secho(
        "'*' means default. See `blastp -help` form more information.", fg="yellow"
    )
