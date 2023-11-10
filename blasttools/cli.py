from __future__ import annotations

from collections.abc import Sequence

import click

from .blastapi import blast_options
from .blastapi import blastall
from .blastapi import BlastConfig
from .blastapi import buildall
from .blastapi import check_ext
from .blastapi import merge_fasta
from .blastapi import read_df
from .blastapi import save_df
from .blastapi import test_save
from .blastapi import toblastdb


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
        click.secho(f"Can't install {REPO}", fg="red", err=True)
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
@click.option("-n", "--nucl", is_flag=True, help="nucleotide blastn")
@click.option(
    "-m",
    "--merge",
    help="merge all fastafiles into one (and create one blast database)",
)
@click.argument(
    "fastas",
    nargs=-1,
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
def build_cmd(
    fastas: Sequence[str],
    builddir: str | None,
    merge: str | None,
    nucl: bool,
) -> None:
    """Build blast databases from fasta files"""
    buildall(fastas, builddir=builddir, merge=merge, blastp=not nucl)


@blast.command(name="blast")
@click.option(
    "--out",
    help="output filename (default is to write <query>.csv)",
    type=click.Path(dir_okay=False),
)
@click.option("--xml", is_flag=True, help="run with xml output")
@click.option(
    "-c",
    "--columns",
    help="space separated list of columns (see columns cmd for a list of valid columns)",
)
@click.option("-n", "--nucl", is_flag=True, help="nucleotide blastn")
@blast_options
@click.argument("query", type=click.Path(exists=True, dir_okay=False))
@click.argument("blastdbs", nargs=-1)
def blast_cmd(
    query: str,
    blastdbs: Sequence[str],
    out: str | None,
    xml: bool,
    columns: str | None,
    nucl: bool,
    # blast options
    best: int,
    with_seq: bool,
    num_threads: int,
    with_description: bool,
    expr: str,
    without_query_seq: bool,
) -> None:
    """Run a blast query over specified databases"""
    from .blastapi import mkheader, has_pdatabase
    from .blastxml import blastall as blastall5

    if len(blastdbs) == 0:
        return

    if out is not None:
        check_ext(out)  # fail early
        test_save(out)

    blastdbs = toblastdb(blastdbs)
    missing = {b for b in blastdbs if not has_pdatabase(b)}
    if missing:
        m = ", ".join(missing)
        raise click.BadParameter(f"missing databases {m}", param_hint="blastdbs")

    myheader = None
    if columns is not None:
        if xml:
            raise click.BadParameter(
                'Can\'t have "--columns" with "--xml"',
                param_hint="xml",
            )
        myheader = mkheader(columns)

    config = BlastConfig(
        best=best,
        with_seq=with_seq,
        header=myheader,
        num_threads=num_threads,
        with_description=with_description,
        expr=expr,
        blastp=not nucl,
        without_query_seq=without_query_seq,
    )
    if xml:
        df = blastall5(query, blastdbs, config=config)
    else:
        df = blastall(query, blastdbs, config=config)
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
        "'*' means default. See `blastp -help` form more information.",
        fg="yellow",
    )


@blast.command(name="concat")
@click.option(
    "--out",
    help="output filename. If not specified, write CSV to stdout",
    type=click.Path(dir_okay=False),
)
@click.argument("dataframes", nargs=-1, type=click.Path(dir_okay=False, exists=True))
def concat_cmd(dataframes: Sequence[str], out: str | None) -> None:
    """Concatenate multiple saved DataFrames"""
    import sys
    import pandas as pd

    if out is not None:
        check_ext(out)
        test_save(out)

    dfs = []
    for df in dataframes:
        res = read_df(df)
        if res is None:
            click.secho(f"Can't read {df}", err=True, bold=True, fg="red")
            continue
        dfs.append(res)

    odf = pd.concat(dfs, axis=0)
    if out is not None:
        save_df(odf, out)
    else:
        odf.to_csv(sys.stdout, index=False)


@blast.command(name="fasta-split")
@click.option("-d", "--directory", type=click.Path(dir_okay=True, file_okay=False))
@click.option("--fmt", help="format of split filenames")
@click.argument("fastafile", type=click.Path(dir_okay=False, exists=True))
@click.argument("batch", type=int)
def fasta_split_cmd(
    fastafile: str,
    batch: int,
    fmt: str | None,
    directory: str | None,
) -> None:
    """Split a fasta file into batches"""
    from .utils import split_fasta

    split_fasta(fastafile, batch, target_dir=directory, fmt=fmt)
