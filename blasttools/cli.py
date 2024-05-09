from __future__ import annotations

from collections.abc import Sequence

import click

from .options import blast_options


@click.group(
    epilog=click.style(
        "Commands for turning blast queries into pandas dataframes.\n",
        fg="magenta",
    ),
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
def fasta_merge_cmd(fasta1: str, fasta2: str, outfasta: str) -> None:
    """merge 2 fasta files based on sequence identity"""
    from .blastapi import fasta_merge

    fasta_merge(fasta1, fasta2, outfasta)


@blast.command(name="fasta-xref")
@click.option(
    "--out",
    help="output filename (default is to write CSV to stdout)",
    type=click.Path(dir_okay=False),
)
@click.argument("fasta1", type=click.Path(exists=True, dir_okay=False))
@click.argument("fasta2", type=click.Path(exists=True, dir_okay=False))
def fasta_xref_cmd(fasta1: str, fasta2: str, out: str | None) -> None:
    """match IDs based on sequence identity"""
    import sys
    from .blastapi import fasta_xref, save_df, test_save

    if out is not None:
        test_save(out)
    df = fasta_xref(fasta1, fasta2)

    if out is not None:
        save_df(df, out)
    else:
        df.to_csv(sys.stdout, index=False)


@blast.command(name="build")
@click.option(
    "-b",
    "--build",
    "builddir",
    type=click.Path(file_okay=False, dir_okay=True),
    help="build all databases in this directory — otherwise in same directory as fasta file(s)",
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
    from .blastapi import buildall

    buildall(fastas, builddir=builddir, merge=merge, blastp=not nucl)


@blast.command(name="blast")
@click.option(
    "--out",
    help="output filename (default is to write <query>.csv)",
    type=click.Path(dir_okay=False),
)
@click.option(
    "-c",
    "--columns",
    help="space or comma separated list of columns (see columns cmd for a list of valid columns)",
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
    with_subject_seq: bool,
    num_threads: int,
    with_description: bool,
    expr: str,
    without_query_seq: bool,
    needs_translation: bool,
) -> None:
    """Run a blast query over specified databases"""
    from .blastapi import (
        mkheader,
        has_pdatabase,
        blastall,
        save_df,
        test_save,
        check_ext,
        toblastdb,
        BlastConfig,
    )

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
        with_subject_seq=with_subject_seq,
        header=myheader,
        num_threads=num_threads,
        with_description=with_description,
        expr=expr,
        blastp=not nucl,
        without_query_seq=without_query_seq,
        xml=xml,
        needs_translation=needs_translation,
    )

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
    from .blastapi import read_df, save_df, test_save, check_ext

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

    odf = pd.concat(dfs, axis=0, ignore_index=True)
    if out is not None:
        save_df(odf, out)
    else:
        odf.to_csv(sys.stdout, index=False)


@blast.command(name="fasta-split")
@click.option(
    "-d",
    "--directory",
    type=click.Path(dir_okay=True, file_okay=False),
    help="directory to write split fastas",
)
@click.option(
    "-Z",
    "--Z",
    "use_null",
    is_flag=True,
    help="separate filenames by '\\0' in output",
)
@click.option("--fmt", help="format of split filenames. Must have a {num} key.")
@click.argument("fastafile", type=click.Path(dir_okay=False, exists=True))
@click.argument("batch", type=int)
def fasta_split_cmd(
    fastafile: str,
    batch: int,
    fmt: str | None,
    directory: str | None,
    use_null: bool,
) -> None:
    """Split a fasta file into batches"""
    from .utils import split_fasta
    from .blastapi import list_out

    ret = split_fasta(fastafile, batch, target_dir=directory, fmt=fmt)
    if ret is None:
        raise click.ClickException(f"Can't split fasta file {fastafile}")
    # So you can do say:
    # fastas=$(blasttools fasta-split fastafile.fa.gz 20000)
    # parallel blasttools blast --out=my{}.pkl ::: $fastas ::: blastdb
    # blasttools concat --out=final.csv my*.pkl && rm my*.pkl
    list_out((p.resolve() for p in ret), use_null=use_null)
