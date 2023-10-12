# see https://plants.ensembl.org/info/data/ftp/index.html
# or asia https://asia.ensembl.org/info/data/ftp/index.html
from __future__ import annotations

from typing import Iterator, Sequence

import ftplib
import glob
import subprocess

from pathlib import Path

import click
import pandas as pd  # type: ignore

from .blastapi import (
    safe_which,
    fasta_to_df,
    doblast6,
    fetch_seq as fetch_seq_raw,
    find_best,
)

from .blastapi import BlastDb

PEP_DIR = "ensemblgenomes/pub/release-{release}/plants/fasta/{plant}/pep"
FTPURL = "ftp.ebi.ac.uk"
ENSEMBL = f"ftp://{FTPURL}/" + PEP_DIR + "/{file}"


def blast_dir(release: int) -> Path:
    return Path(f"blast-{release}")


def find_fasta_names(plants: Sequence[str], release: int) -> Iterator[str | None]:
    with ftplib.FTP(FTPURL) as ftp:
        ftp.login()
        for plant in plants:
            dname = PEP_DIR.format(release=release, plant=plant)
            for n in ftp.nlst(dname):
                if n.endswith(".pep.all.fa.gz"):
                    yield n[len(dname) + 1 :]
                    break
            else:
                yield None


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
    queryfasta: str,
    plant: str,
    release: int,
) -> pd.DataFrame:
    blastdir = blast_dir(release)

    return doblast6(queryfasta, str(blastdir / plant))


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
        ftp.cwd(f"ensemblgenomes/pub/release-{release}/plants/fasta/")
        return list(ftp.nlst())


def fetch_fastas(plants: Sequence[str], release: int) -> None:
    bd = blast_dir(release)

    bd.mkdir(parents=True, exist_ok=True)

    for plant, filename in zip(plants, find_fasta_names(plants, release)):
        if filename is None:
            print(f"can't find fasta for {plant}!")
            continue

        if (bd / filename).exists():
            continue
        print("fetching", filename)
        r = fetch_fasta(plant, filename, release=release)
        if r.returncode:
            print(r)


def build(species: Sequence[str], release: int):
    blastdir = blast_dir(release)
    blastdir.mkdir(parents=True, exist_ok=True)

    plants = [
        (p, f)
        for p, f in zip(species, find_fasta_names(species, release))
        if f is not None
    ]

    for plant, filename in plants:
        if has_fasta(blastdir, filename):
            continue
        print(f"fetching {filename} for release: {release}")
        r = fetch_fasta(plant, filename, release)
        if r.returncode:
            click.secho(f"failed to fetch {filename}", fg="red", err=True)

    for plant, filename in plants:
        if has_blast_db(blastdir, plant):
            continue
        print("creating blast db for", plant)
        ok = mkblast(plant, filename, release)
        if not ok:
            print(f"failed to create blast db for {plant}")


def blastall(
    query: str, species: Sequence[str], release: int, best: int, with_seq: bool
) -> pd.DataFrame:
    df = fasta_to_df(query)
    res = []

    build(species, release)
    for plant in species:
        rdf = doblast(query, plant, release=release)

        if with_seq and "saccver" in rdf.columns:
            saccver = list(rdf["saccver"])
            sdf = fetch_seq_df(saccver, plant, release)
            rdf = pd.merge(rdf, sdf, left_on="saccver", right_on="saccver")

        myrdf = find_best(rdf, df, nevalues=best)
        myrdf["species"] = plant
        res.append(myrdf)
    ddf = pd.concat(res, axis=0)

    return ddf
