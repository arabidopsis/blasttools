# see https://plants.ensembl.org/info/data/ftp/index.html
# or asia https://asia.ensembl.org/info/data/ftp/index.html
from __future__ import annotations

import ftplib
import glob
import subprocess
from dataclasses import dataclass
from pathlib import Path
from shutil import which
from typing import Iterator
from typing import Sequence

import click
import pandas as pd

from .blastapi import Blast6
from .blastapi import BlastConfig
from .blastapi import BlastDb
from .blastapi import check_expr
from .blastapi import fasta_to_df
from .blastapi import fetch_seq as fetch_seq_raw
from .blastapi import find_best
from .blastapi import remove_files
from .blastapi import safe_which

FASTAS_DIR = "ensemblgenomes/pub/release-{release}/plants/fasta/"
TOP = "ensemblgenomes/pub/release-{release}/plants/"
PEP_DIR = TOP + "fasta/{plant}/pep"
FTPURL = "ftp.ebi.ac.uk"
ENSEMBL = f"ftp://{FTPURL}/" + PEP_DIR + "/{file}"
SPECIES_TSV = "species_EnsemblPlants.txt"


def blast_dir(release: int) -> Path:
    return Path(f"ensemblblast-{release}")


@dataclass
class FileInfo:
    species: str
    fasta: str | None


def find_fasta_names(plants: Sequence[str], release: int) -> Iterator[FileInfo]:
    from .config import FTP_TIMEOUT

    with ftplib.FTP(FTPURL, timeout=FTP_TIMEOUT) as ftp:
        ftp.login()
        for plant in plants:
            dname = PEP_DIR.format(release=release, plant=plant)
            for n in ftp.nlst(dname):
                if n.endswith(".pep.all.fa.gz"):
                    yield FileInfo(plant, n[len(dname) + 1 :])
                    break
            else:
                yield FileInfo(plant, None)


def _wgetcmd(filename: str, quiet: bool) -> list[str]:
    wget = which("wget")
    if wget is not None:
        q = ["-q"] if quiet else []
        cmds = [
            wget,
            *q,
            "-O",
            filename,
        ]
    else:
        curl = safe_which("curl")
        q = ["--silent"] if quiet else []
        cmds = [
            curl,
            *q,
            "-o",
            filename,
        ]
    return cmds


def fetch_fasta(
    plant: str,
    filename: str,
    release: int,
    *,
    quiet: bool = False,
) -> subprocess.CompletedProcess[bytes]:
    url = ENSEMBL.format(release=release, plant=plant, file=filename)
    return fetch_file(filename, url, release, quiet=quiet)


def fetch_file(
    filename: str,
    url: str,
    release: int,
    *,
    quiet: bool = False,
) -> subprocess.CompletedProcess[bytes]:
    cwd = blast_dir(release)
    if not cwd.exists():
        cwd.mkdir(parents=True, exist_ok=True)

    cmds = _wgetcmd(filename, quiet)
    cmds.append(url)
    resp = subprocess.run(
        cmds,
        cwd=str(cwd),
        check=False,
    )

    if resp.returncode != 0:
        remove_files([cwd / filename])
    return resp


def mkblast(plant: str, fastafile: str | Path, release: int) -> bool:
    directory = blast_dir(release)
    db = directory / plant
    bdb = BlastDb(db)
    return bdb.run(directory / fastafile)


# def doblast(
#     queryfasta: str,
#     plant: str,
#     release: int,
#     header: Sequence[str] | None = None,
#     *,
#     path: str | None,
#     num_threads: int = 1,
# ) -> pd.DataFrame:
#     blastdir = blast_dir(release)
#     if path is not None:
#         blastdir = Path(path) / blastdir

#     return doblast6(
#         queryfasta,
#         str(blastdir / plant),
#         header=header,
#         num_threads=num_threads,
#     )


def has_blast_db(blastdir: Path, plant: str) -> bool:
    return len(glob.glob(str(blastdir / (plant + ".*")))) > 0


def has_fasta(blastdir: Path, filename: str) -> bool:
    return (blastdir / filename).exists()


def fetch_seq_df(
    seqids: Sequence[str],
    plant: str,
    release: int,
    *,
    path: str | None,
) -> pd.DataFrame:
    blastdir = blast_dir(release)
    if path is not None:
        blastdir = Path(path) / blastdir
    blastdb = str(blastdir / plant)

    return pd.DataFrame(
        [
            {"saccver": rec.id, "subject_seq": str(rec.seq)}
            for rec in fetch_seq_raw(seqids, blastdb)
        ],
    )


def find_species(release: int) -> list[str]:
    from .config import FTP_TIMEOUT

    with ftplib.FTP(FTPURL, timeout=FTP_TIMEOUT) as ftp:
        ftp.login()
        ftp.cwd(FASTAS_DIR.format(release=release))
        return list(ftp.nlst())


def find_species_file(release: int, *, quiet: bool = False) -> pd.DataFrame | None:
    url = f"ftp://{FTPURL}/{TOP.format(release=release)}/{SPECIES_TSV}"
    resp = fetch_file(SPECIES_TSV, url, release, quiet=quiet)
    if resp.returncode:
        return None
    bd = blast_dir(release)
    return pd.read_csv(bd / SPECIES_TSV, sep="\t")


def fetch_fastas(plants: Sequence[str], release: int) -> None:
    bd = blast_dir(release)

    bd.mkdir(parents=True, exist_ok=True)

    for info in find_fasta_names(plants, release):
        if info.fasta is None:
            click.echo(f"Can't find fasta for {info.species}!", err=True)
            continue

        if (bd / info.fasta).exists():
            continue
        click.echo(f"fetching {info.fasta}")
        r = fetch_fasta(info.species, info.fasta, release=release)
        if r.returncode:
            click.secho(f"failed to fetch {info.fasta}: {r}", fg="red", err=True)


def build(species: Sequence[str], release: int, *, path: str | None = None) -> bool:
    blastdir = blast_dir(release)
    if path is not None:
        blastdir = Path(path) / blastdir
    blastdir.mkdir(parents=True, exist_ok=True)

    plants = list(find_fasta_names(species, release=release))
    ret = True
    for info in plants:
        # we are OK if we have a blast database
        # but *NO* fasta file
        if has_blast_db(blastdir, info.species):
            continue
        if info.fasta is None:
            continue
        if has_fasta(blastdir, info.fasta):
            continue
        click.echo(f"fetching {info.fasta} for release: {release}")
        r = fetch_fasta(info.species, info.fasta, release)
        if r.returncode:
            ret = False
            click.secho(
                f"failed to fetch {info.fasta}",
                fg="red",
                bold=True,
                err=True,
            )

    for info in plants:
        if has_blast_db(blastdir, info.species):
            continue
        if info.fasta is None:
            ret = False
            click.secho(
                f'can\'t find fasta file for "{info.species}"',
                fg="red",
                bold=True,
                err=True,
            )
            continue
        click.echo(f"creating blast db for {info.species}")
        ok = mkblast(info.species, info.fasta, release)
        if not ok:
            ret = False
            click.secho(
                f"failed to create blast db for {info.species}",
                fg="red",
                bold=True,
                err=True,
            )
    return ret


def blastall(
    queryfasta: str | Path,
    species: Sequence[str],
    release: int,
    *,
    path: str | None = None,
    config: BlastConfig = BlastConfig(),
) -> pd.DataFrame:
    b6 = Blast6(config.header, num_threads=config.num_threads, blastp=config.blastp)

    check_expr(b6.header, config.expr)  # fail early

    qdf = fasta_to_df(queryfasta, with_description=config.with_description)
    if not qdf["id"].is_unique:
        raise click.ClickException(
            f'sequences IDs are not unique for query file "{queryfasta}"',
        )
    if config.without_query_seq:
        qdf.drop(columns=["seq"], inplace=True)
    res = []

    ok = build(species, release, path=path)
    if not ok:
        raise click.ClickException("Can't build blast databases(s)")

    blastdir = blast_dir(release)
    if path is not None:
        blastdir = Path(path) / blastdir
    for plant in species:
        rdf = b6.run(queryfasta, blastdir / plant)

        if config.with_seq and "saccver" in rdf.columns:
            saccver = list(rdf["saccver"])
            sdf = fetch_seq_df(saccver, plant, release, path=path)
            rdf = pd.merge(rdf, sdf, left_on="saccver", right_on="saccver")

        myrdf = find_best(rdf, qdf, nevalues=config.best, evalue_col=config.expr)

        myrdf["species"] = plant
        res.append(myrdf)
    ddf = pd.concat(res, axis=0)

    return ddf


def available_species(release: int) -> list[str]:
    db = blast_dir(release)
    ret = set()
    for path in db.glob("*.pot"):
        n, _ = path.name.split(".", maxsplit=1)
        ret.add(n)
    return list(ret)


def orthologs(
    query_species: str,
    subject_species: str,
    release: int,
    config: BlastConfig = BlastConfig(),
) -> pd.DataFrame | None:
    if not build([query_species, subject_species], release):
        return None
    names = list(find_fasta_names([query_species], release=release))
    if not names:
        return None
    fasta = names[0].fasta
    if fasta is None:
        return None
    blastdir = blast_dir(release)
    query_fasta = blastdir / fasta
    if not query_fasta.exists():
        return None

    return blastall(query_fasta, [subject_species], release=release, config=config)
