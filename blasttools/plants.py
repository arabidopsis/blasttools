# see https://plants.ensembl.org/info/data/ftp/index.html
# or asia https://asia.ensembl.org/info/data/ftp/index.html
from __future__ import annotations

from typing import Iterator, Sequence
from dataclasses import dataclass
import ftplib
import glob
import subprocess

from pathlib import Path

import click
import pandas as pd

from .blastapi import (
    safe_which,
    fasta_to_df,
    doblast6,
    fetch_seq as fetch_seq_raw,
    find_best,
)

from .blastapi import BlastDb

FASTAS_DIR = "ensemblgenomes/pub/release-{release}/plants/fasta/"
PEP_DIR = "ensemblgenomes/pub/release-{release}/plants/fasta/{plant}/pep"
FTPURL = "ftp.ebi.ac.uk"
ENSEMBL = f"ftp://{FTPURL}/" + PEP_DIR + "/{file}"


def blast_dir(release: int) -> Path:
    return Path(f"ensemblblast-{release}")


@dataclass
class FileInfo:
    species: str
    fasta: str | None


def find_fasta_names(plants: Sequence[str], release: int) -> Iterator[FileInfo]:
    with ftplib.FTP(FTPURL) as ftp:
        ftp.login()
        for plant in plants:
            dname = PEP_DIR.format(release=release, plant=plant)
            for n in ftp.nlst(dname):
                if n.endswith(".pep.all.fa.gz"):
                    yield FileInfo(plant, n[len(dname) + 1 :])
                    break
            else:
                yield FileInfo(plant, None)


def fetch_fasta(
    plant: str, filename: str, release: int
) -> subprocess.CompletedProcess[bytes]:
    curl = safe_which("curl")
    r = subprocess.run(
        [
            curl,
            "-o",
            filename,
            ENSEMBL.format(release=release, plant=plant, file=filename),
        ],
        cwd=str(blast_dir(release)),
        check=False,
    )
    return r


def mkblast(plant: str, fastafile: str, release: int) -> bool:
    db = blast_dir(release) / plant
    bdb = BlastDb(str(db))
    return bdb.run(fastafile)


def doblast(
    queryfasta: str, plant: str, release: int, header: Sequence[str] | None = None
) -> pd.DataFrame:
    blastdir = blast_dir(release)

    return doblast6(queryfasta, str(blastdir / plant), header=header)


def has_blast_db(blastdir: Path, plant: str) -> bool:
    return len(glob.glob(str(blastdir / (plant + ".*")))) > 0


def has_fasta(blastdir: Path, filename: str) -> bool:
    return (blastdir / filename).exists()


def fetch_seq_df(seqids: Sequence[str], plant: str, release: int) -> pd.DataFrame:
    blastdir = blast_dir(release)
    blastdb = str(blastdir / plant)

    return pd.DataFrame(
        [
            {"saccver": rec.id, "subject_seq": str(rec.seq)}
            for rec in fetch_seq_raw(seqids, blastdb)
        ]
    )


def find_species(release: int) -> list[str]:
    with ftplib.FTP(FTPURL) as ftp:
        ftp.login()
        ftp.cwd(FASTAS_DIR.format(release=release))
        return list(ftp.nlst())


def fetch_fastas(plants: Sequence[str], release: int) -> None:
    bd = blast_dir(release)

    bd.mkdir(parents=True, exist_ok=True)

    for info in find_fasta_names(plants, release):
        if info.fasta is None:
            print(f"can't find fasta for {info.species}!")
            continue

        if (bd / info.fasta).exists():
            continue
        print("fetching", info.fasta)
        r = fetch_fasta(info.species, info.fasta, release=release)
        if r.returncode:
            print(r)


def build(species: Sequence[str], release: int) -> bool:
    blastdir = blast_dir(release)
    blastdir.mkdir(parents=True, exist_ok=True)

    plants = list(find_fasta_names(species, release=release))
    ret = True
    for info in plants:
        if info.fasta is None:
            continue
        if has_fasta(blastdir, info.fasta):
            continue
        click.echo(f"fetching {info.fasta} for release: {release}")
        r = fetch_fasta(info.species, info.fasta, release)
        if r.returncode:
            ret = False
            click.secho(f"failed to fetch {info.fasta}", fg="red", err=True, bold=True)

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
    queryfasta: str,
    species: Sequence[str],
    release: int,
    best: int,
    with_seq: bool,
    header: Sequence[str] | None = None,
) -> pd.DataFrame:
    df = fasta_to_df(queryfasta)
    if not df["id"].is_unique:
        raise click.ClickException(
            f'sequences IDs are not unique for query file "{queryfasta}"'
        )
    res = []

    ok = build(species, release)
    if not ok:
        raise click.ClickException("can't build blast databases(s)")

    for plant in species:
        rdf = doblast(queryfasta, plant, release=release, header=header)

        if with_seq and "saccver" in rdf.columns:
            saccver = list(rdf["saccver"])
            sdf = fetch_seq_df(saccver, plant, release)
            rdf = pd.merge(rdf, sdf, left_on="saccver", right_on="saccver")

        myrdf = find_best(rdf, df, nevalues=best)
        myrdf["species"] = plant
        res.append(myrdf)
    ddf = pd.concat(res, axis=0)

    return ddf
