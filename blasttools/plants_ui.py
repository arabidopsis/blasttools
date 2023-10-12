from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence
import click
from .config import RELEASE

from .blast import (
    build,
    blastall,
    fetch_fastas,
    find_fasta_names,
    find_species,
)

from .cli import blast
from .blastapi import save_df


@dataclass
class Config:
    release: int = RELEASE


pass_config = click.make_pass_decorator(Config, ensure=True)


@blast.group(help=click.style("blast commands that understand Ensembl", fg="magenta"))
@click.option(
    "-r", "--release", default=RELEASE, show_default=True, help="release number"
)
@click.pass_context
def plants(ctx: click.Context, release: int):
    """Run blast on ensembl plant genomes"""
    ctx.obj = Config(release=release)


@plants.command(name="build")
@click.argument("species", nargs=-1)
@pass_config
def build_cmd(cfg: Config, species: Sequence[str]):
    """Download and build blast databases"""
    build(species, release=cfg.release)


@plants.command(name="blast")
@click.option("--out", help="output file", type=click.Path(dir_okay=False))
@click.option("--best", default=2, help="best (lowest) evalues [=0 take all]")
@click.option("--with_seq", is_flag=True, help="added sequence data to output")
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
):
    """Run blast"""
    df = blastall(query, species, release=cfg.release, best=best, with_seq=with_seq)
    if out is None:
        out = query + ".csv"

    click.secho(f"writing {out}", fg="green")
    save_df(df, out, index=False)


@plants.command(name="fetch-fastas")
@click.argument("species", nargs=-1)
@pass_config
def fetch_fasta_cmd(cfg: Config, species: Sequence[str]):
    """Fetch fasta files"""
    fetch_fastas(species, release=cfg.release)


@plants.command()
@click.argument("species", nargs=-1)
@pass_config
def fasta_names(cfg: Config, species: Sequence[str]) -> None:
    """Find fasta names for plant species"""
    for plant, name in zip(species, find_fasta_names(species, release=cfg.release)):
        click.secho(f"{plant}: {name}")


@plants.command(name="species")
@pass_config
def species_cmd(cfg: Config) -> None:
    """Available species at ensembl"""
    sl = find_species(cfg.release)
    sl = sorted(sl)
    for s in sl:
        print(s)
