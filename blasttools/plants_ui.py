from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import click

from .cli import blast
from .config import RELEASE
from .options import blast_options


@dataclass
class Config:
    release: int = RELEASE


pass_config = click.make_pass_decorator(Config, ensure=True)

EPILOG = """
Fasta files and blast databases files are identified by their species name (see `blasttools plants species`).
All files are stored in subdirectories identified by the `--release` number.
e.g. `ensemblplants-{release}/`. Fasta files are downloaded from "ftp://ftp.ebi.ac.uk" as needed.
"""


@blast.group(
    help=click.style("blast commands that understand Ensembl", fg="magenta"),
    epilog=EPILOG,
)
@click.option(
    "-r",
    "--release",
    default=RELEASE,
    show_default=True,
    help="release number",
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
    from .plants import build

    build(species, release=cfg.release, path=builddir)


@plants.command(name="blast")
@click.option(
    "--out",
    help="output filename (default is to write <query>.csv)",
    type=click.Path(dir_okay=False),
)
@click.option(
    "-n",
    "--names",
    is_flag=True,
    help="use descriptive column names in output",
)
@click.option(
    "-c",
    "--columns",
    help="space or comma separated list of columns (see columns cmd for a list of valid columns)."
    "Prefix with a '+' to add the columns or '-' to remove them",
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
@blast_options
@click.argument("query", type=click.Path(exists=True, dir_okay=False))
@click.argument("species", nargs=-1)
@pass_config
def blast_cmd(
    cfg: Config,
    query: str,
    species: Sequence[str],
    out: str | None,
    names: bool,
    columns: str | None,
    builddir: str | None,
    all_db: bool,
    # blast options
    best: int,
    with_subject_seq: bool,
    num_threads: int,
    with_description: bool,
    expr: str,
    without_query_seq: bool,
    xml: bool,
    needs_translation: bool,
) -> None:
    """Run blast on query fasta file"""
    from .plants import available_species
    from .columns import VALID
    from .blastapi import BlastConfig
    from .blastapi import check_ext
    from .blastapi import mkheader
    from .blastapi import save_df
    from .blastapi import test_save
    from .plants import blastall

    if len(species) == 0 and not all_db:
        return
    if out is not None:
        check_ext(out)
        test_save(out)

    myheader = None
    if columns is not None:
        if xml:
            raise click.BadParameter(
                'Can\'t have "--columns" with "--xml"',
                param_hint="xml",
            )
        myheader = mkheader(columns)

    if all_db:
        species = available_species(cfg.release)
        if len(species) == 0:
            return
    config = BlastConfig(
        best=best,
        with_subject_seq=with_subject_seq,
        header=myheader,
        num_threads=num_threads,
        with_description=with_description,
        expr=expr,
        without_query_seq=without_query_seq,
        xml=xml,
        needs_translation=needs_translation,
    )
    df = blastall(query, species, release=cfg.release, path=builddir, config=config)
    if out is None:
        out = query + ".csv"
    if names:
        df.rename(columns=VALID, inplace=True)
    click.secho(f"writing {out}", fg="green")
    save_df(df, out, index=False)


@plants.command(name="fasta-fetch")
@click.option(
    "--cdna",
    is_flag=True,
    help="download the cdna (instead of peptide) version",
)
@click.argument("species", nargs=-1)
@pass_config
def fasta_fetch_cmd(cfg: Config, species: Sequence[str], cdna: bool) -> None:
    """Download fasta files from FTP site"""
    from .plants import fetch_fastas

    if not species:
        return
    download_dir = fetch_fastas(
        species,
        release=cfg.release,
        seqtype="cdna" if cdna else "pep",
    )
    dd = click.style(str(download_dir), fg="blue")
    s = "s" if len(species) > 1 else ""
    click.echo(f"downloaded file{s} into: {dd}/")


@plants.command("fasta-filenames")
@click.option("-f", "--full", is_flag=True, help="show full URL to file")
@click.argument("species", nargs=-1)
@pass_config
def fasta_filenames_cmd(
    cfg: Config,
    species: Sequence[str],
    full: bool,
) -> None:
    """Find fasta filenames for plant species"""
    from .plants import ENSEMBL
    from .plants import find_fasta_names

    for info in find_fasta_names(species, release=cfg.release):
        if info.fasta is None:
            click.secho(f"{info.species}: no fasta!", fg="red", bold=True)
        else:
            if full:
                fasta = ENSEMBL.format(
                    release=cfg.release,
                    plant=info.species,
                    file=info.fasta,
                )
                click.secho(f"{info.species}: {fasta}")
            else:
                click.secho(f"{info.species}: {info.fasta}")


@plants.command(name="species")
@pass_config
def species_cmd(cfg: Config) -> None:
    """Available species at Ensembl"""
    from .plants import find_species

    sl = find_species(cfg.release)
    for s in sorted(sl):
        click.echo(s)


@plants.command(name="releases")
@click.option("--max", "max_only", help="max version", is_flag=True)
def releases_cmd(max_only: bool) -> None:
    """Available species at Ensembl"""
    from .plants import find_releases

    rl = find_releases()
    if max_only:
        rl = rl[-1:]
    for s in rl:
        click.echo(s)


@plants.command(name="ortholog")
@click.option(
    "--out",
    help="output filename (default is to write <query>.csv)",
    type=click.Path(dir_okay=False),
)
@blast_options
@click.argument("query_species")
@click.argument("subject_species")
@pass_config
def ortholog_cmd(
    cfg: Config,
    query_species: str,
    subject_species: str,
    out: str | None,
    best: int,
    with_subject_seq: bool,
    num_threads: int,
    with_description: bool,
    expr: str,
    without_query_seq: bool,
    xml: bool,
    needs_translation: bool,
) -> None:
    """Create an ortholog DataFrame between two species"""
    from .blastapi import BlastConfig, check_ext, test_save, save_df
    from .plants import orthologs

    config = BlastConfig(
        best=best,
        with_subject_seq=with_subject_seq,
        header=None,
        num_threads=num_threads,
        with_description=with_description,
        expr=expr,
        without_query_seq=without_query_seq,
        xml=xml,
        needs_translation=needs_translation,
    )
    if out is not None:
        check_ext(out)
        test_save(out)
    df = orthologs(
        query_species,
        subject_species,
        release=cfg.release,
        config=config,
    )
    if df is None:
        raise click.ClickException("Can't build databases!")
    if out is None:
        out = f"{query_species}-{subject_species}.csv"
    click.secho(f'writing "{out}"', fg="green")
    save_df(df, out, index=False)
