from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import click
from .config import RELEASE
from .blastapi import EVALUE, BlastConfig

from .plants import (
    build,
    blastall,
    fetch_fastas,
    find_fasta_names,
    find_species,
)

from .cli import blast
from .blastapi import save_df, mkheader


@dataclass
class Config:
    release: int = RELEASE


pass_config = click.make_pass_decorator(Config, ensure=True)


@blast.group(help=click.style("blast commands that understand Ensembl", fg="magenta"))
@click.option(
    "-r", "--release", default=RELEASE, show_default=True, help="release number"
)
@click.pass_context
def plants(ctx: click.Context, release: int) -> None:
    """Run blast on ensembl plant genomes"""
    ctx.obj = Config(release=release)


BUILD = """
Blast databases will be named after the species name and
placed in the directory 'ensembleblast-{release}' (which will be created)
(e.g. 'ensemblblast-57/zea_mays.p*')
"""


@plants.command(name="build", epilog=click.style(BUILD, fg="magenta"))
@click.option(
    "-b",
    "--build",
    "builddir",
    type=click.Path(file_okay=False, dir_okay=True),
    help="build all 'ensembleblast-{release}' directories in this directory",
)
@click.argument("species", nargs=-1)
@pass_config
def build_cmd(cfg: Config, species: Sequence[str], builddir: str | None) -> None:
    """Download and build blast databases"""
    build(species, release=cfg.release, path=builddir)


@plants.command(name="blast")
@click.option(
    "--out",
    help="output filename (default is to write <query>.csv)",
    type=click.Path(dir_okay=False),
)
@click.option(
    "--best", default=0, help="best (lowest) evalues [=0 take all]  (see also --expr)"
)
@click.option("--with-seq", is_flag=True, help="add sequence data to output")
@click.option(
    "-n", "--names", is_flag=True, help="use descriptive column names in output"
)
@click.option(
    "-d", "--with-description", is_flag=True, help="include query description"
)
@click.option(
    "--expr", help="expression to minimize", default=EVALUE, show_default=True
)
@click.option(
    "-c",
    "--columns",
    help="space separated list of columns (see columns cmd for a list of valid columns)",
)
@click.option(
    "-a",
    "--all",
    "all_db",
    is_flag=True,
    help="try all available databases",
)
@click.option(
    "-b",
    "--build",
    "builddir",
    type=click.Path(file_okay=False, dir_okay=True),
    help="find all 'ensembleblast-{release}' directories in this directory",
)
@click.option("-t", "--num-threads", help="number of threads to use", default=1)
@click.argument("query", type=click.Path(exists=True, dir_okay=False))
@click.argument("species", nargs=-1)
@pass_config
def blast_cmd(
    cfg: Config,
    query: str,
    species: Sequence[str],
    best: int,
    with_seq: bool,
    out: str | None,
    names: bool,
    columns: str | None,
    num_threads: int,
    builddir: str | None,
    all_db: bool,
    with_description: bool,
    expr: str,
) -> None:
    """Run blast on query fasta file"""
    from .blastapi import check_ext
    from .plants import available_species
    from .columns import VALID

    if len(species) == 0 and not all_db:
        return
    if out is not None:
        check_ext(out)

    myheader = None
    if columns is not None:
        myheader = mkheader(columns)

    if all_db:
        species = available_species(cfg.release)
        if len(species) == 0:
            return
    config = BlastConfig(
        best=best,
        with_seq=with_seq,
        header=myheader,
        num_threads=num_threads,
        with_description=with_description,
        expr=expr,
    )
    df = blastall(query, species, release=cfg.release, path=builddir, config=config)
    if out is None:
        out = query + ".csv"
    if names:
        df.rename(columns=VALID, inplace=True)
    click.secho(f"writing {out}", fg="green")
    save_df(df, out, index=False)


@plants.command(name="fasta-fetch")
@click.argument("species", nargs=-1)
@pass_config
def fasta_fetch_cmd(cfg: Config, species: Sequence[str]) -> None:
    """Download fasta files from FTP site"""
    fetch_fastas(species, release=cfg.release)


@plants.command()
@click.option("-f", "--full", is_flag=True, help="show full URL to file")
@click.argument("species", nargs=-1)
@pass_config
def fasta_filenames(cfg: Config, species: Sequence[str], full: bool) -> None:
    """Find fasta filenames for plant species"""
    from .plants import ENSEMBL

    for info in find_fasta_names(species, release=cfg.release):
        if info.fasta is None:
            click.secho(f"{info.species}: no fasta!", fg="red", bold=True)
        else:
            if full:
                fasta = ENSEMBL.format(
                    release=cfg.release, plant=info.species, file=info.fasta
                )
                click.secho(f"{info.species}: {fasta}")
            else:
                click.secho(f"{info.species}: {info.fasta}")


@plants.command(name="species")
@pass_config
def species_cmd(cfg: Config) -> None:
    """Available species at Ensembl"""
    sl = find_species(cfg.release)
    for s in sorted(sl):
        click.echo(s)
