from __future__ import annotations
import click
from dataclasses import dataclass
from .config import RELEASE
from typing import Sequence

from .blast import (
    build,
    blastall,
    fetch_fastas,
    find_fasta_names,
    find_species,
    merge_fasta,
)


@dataclass
class Config:
    release: int = RELEASE


pass_config = click.make_pass_decorator(Config, ensure=True)


@click.group()
@click.option("-r", "--release", default=RELEASE, show_default=True, help="release number")
@click.version_option()
@click.pass_context
def blast(ctx: click.Context, release: int):
    ctx.obj = Config(release=release)


@blast.command(name="build")
@click.argument("species", nargs=-1)
@pass_config
def build_cmd(cfg: Config, species: Sequence[str]):
    """Download and build blast databases"""
    build(species, release=cfg.release)


@blast.command(name="blast")
@click.option("--best", default=2, help="best (lowest) evalues [=0 take all]")
@click.option("--with_seq", is_flag=True, help="added sequence data to output")
@click.argument("query", type=click.Path(exists=True, dir_okay=False))
@click.argument("species", nargs=-1)
@pass_config
def blast_cmd(
    cfg: Config, query: str, species: Sequence[str], best: int, with_seq: bool
):
    """Run blast"""
    df = blastall(query, species, release=cfg.release, best=best, with_seq=with_seq)
    out = query + ".csv"
    click.secho(f"writing {out}", fg="green")
    df.to_csv(out, index=False)


@blast.command(name="fetch-fastas")
@click.argument("species", nargs=-1)
@pass_config
def fetch_fasta_cmd(cfg: Config, species: Sequence[str]):
    """Fetch fasta files"""
    fetch_fastas(species, release=cfg.release)


@blast.command()
@click.argument("species", nargs=-1)
@pass_config
def fasta_names(cfg: Config, species: Sequence[str]) -> None:
    """Find fasta names for plant species"""
    for plant, name in zip(species, find_fasta_names(species, release=cfg.release)):
        click.secho(f"{plant}: {name}")


@blast.command(name="species")
@pass_config
def species_cmd(cfg: Config) -> None:
    """Available species at ensembl"""
    sl = find_species(cfg.release)
    sl = sorted(sl)
    for s in sl:
        print(s)


@blast.command(name="merge-fasta")
@click.argument("fasta1", type=click.Path(exists=True, dir_okay=False))
@click.argument("fasta2", type=click.Path(exists=True, dir_okay=False))
@click.argument("outfasta", type=click.Path(exists=False, dir_okay=False))
def merge_fasta_cmd(fasta1: str, fasta2: str, outfasta: str) -> None:
    """merge 2 fasta files"""
    merge_fasta(fasta1, fasta2, outfasta)
