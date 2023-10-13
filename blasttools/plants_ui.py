from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import click
from .config import RELEASE

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


@plants.command(name="build")
@click.argument("species", nargs=-1)
@pass_config
def build_cmd(cfg: Config, species: Sequence[str]) -> None:
    """Download and build blast databases"""
    build(species, release=cfg.release)


@plants.command(name="blast")
@click.option("--out", help="output file", type=click.Path(dir_okay=False))
@click.option("--best", default=0, help="best (lowest) evalues [=0 take all]")
@click.option("--with-seq", is_flag=True, help="add sequence data to output")
@click.option(
    "-n", "--names", is_flag=True, help="use descriptive column names in output"
)
@click.option("-h", "--header", help="headers")
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
    header: str | None,
) -> None:
    """Run blast on query fasta file"""
    from .columns import VALID
    if len(species) == 0:
        return

    myheader = None
    if header is not None:
        myheader = mkheader(header)

    df = blastall(
        query,
        species,
        release=cfg.release,
        best=best,
        with_seq=with_seq,
        header=myheader,
    )
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
@click.argument("species", nargs=-1)
@pass_config
def fasta_filenames(cfg: Config, species: Sequence[str]) -> None:
    """Find fasta filenames for plant species"""
    for info in find_fasta_names(species, release=cfg.release):
        if info.fasta is None:
            click.secho(f"{info.species}: no fasta!", fg="red", bold=True)
        else:
            click.secho(f"{info.species}: {info.fasta}")


@plants.command(name="species")
@pass_config
def species_cmd(cfg: Config) -> None:
    """Available species at ensembl"""
    sl = find_species(cfg.release)
    for s in sorted(sl):
        print(s)
