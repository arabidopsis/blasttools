#!/usr/bin/env python
# coding: utf-8

# ## Blasting for Nathan
#
# see https://plants.ensembl.org/info/data/ftp/index.html
# or asia https://asia.ensembl.org/info/data/ftp/index.html
from __future__ import annotations

from typing import Iterator, Sequence


import gzip
import os
import ftplib
import glob
import subprocess
from shutil import which
from pathlib import Path
from uuid import uuid4

import pandas as pd  # type: ignore

from Bio import SeqIO  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio.Seq import Seq  # type: ignore

PEP_DIR = "ensemblgenomes/pub/release-{release}/plants/fasta/{plant}/pep"
FTPURL = "ftp.ebi.ac.uk"
ENSEMBL = f"ftp://{FTPURL}/" + PEP_DIR + "/{file}"


curl = which("curl")
gunzip = which("gunzip")

makeblastdb = which("makeblastdb")
blastp = which("blastp")
blastdbcmd = which("blastdbcmd")

assert curl and makeblastdb and gunzip and blastp and blastdbcmd, "missing executables"


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
    assert curl is not None
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


def mkblast(plant: str, fastafile: str, release: int) -> subprocess.Popen[bytes]:
    assert gunzip is not None and makeblastdb is not None
    # -parse_seqids so that we can get sequences out from ids with blastdbcmd
    bd = blast_dir(release)
    with subprocess.Popen(
        [gunzip, "--stdout", fastafile], stdout=subprocess.PIPE, cwd=str(bd)
    ) as p1:
        with subprocess.Popen(
            [
                makeblastdb,
                "-in",
                "-",
                "-out",
                plant,
                "-input_type=fasta",
                "-dbtype",
                "prot",
                "-title",
                plant,
                "-parse_seqids",
            ],
            stdin=p1.stdout,
            cwd=str(bd),
        ) as p2:
            if p1.stdout:
                p1.stdout.close()
            p2.wait()
            p1.wait()
            return p2


def read_fasta(blastdir: Path, gzfasta: str) -> Iterator[SeqRecord]:
    yield from read_fastax(str(blastdir / gzfasta))


def read_fastax(path: str) -> Iterator[SeqRecord]:
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fp:
            yield from SeqIO.parse(fp, "fasta")
    else:
        with open(path, "rt", encoding="utf-8") as fp:
            yield from SeqIO.parse(fp, "fasta")


def fasta_to_df(path: str) -> pd.DataFrame:
    return pd.DataFrame(
        [dict(id=rec.id, seq=str(rec.seq).upper()) for rec in read_fastax(path)]
    )


# See `blastp -help | less`
# or http://scikit-bio.org/docs/0.5.4/generated/skbio.io.format.blast6.html
# https://www.ncbi.nlm.nih.gov/books/NBK279684/ seems out of date
# default header for -outfmt 6
HEADER = (
    "qaccver",
    "saccver",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)

QUERY = HEADER[0]
EVALUE = HEADER[-2]

RENAMES = {"length": "alignment length"}


def doblast(
    queryfasta: str,
    plant: str,
    release: int,
    header: Sequence[str] = HEADER,
) -> pd.DataFrame:
    blastdir = blast_dir(release)

    return doblastx(queryfasta, str(blastdir / plant), header=header)


def doblastx(
    queryfasta: str,
    blastdb: str,
    header: Sequence[str] = HEADER,
) -> pd.DataFrame:
    assert blastp is not None
    outfmt = f'6 {" ".join(header)}'
    out = f"{queryfasta}.tsv"
    r = subprocess.run(
        [
            blastp,
            "-outfmt",
            outfmt,
            "-query",
            queryfasta,
            "-db",
            blastdb,
            "-out",
            out,
        ],
        check=False,
    )
    if r.returncode:
        raise RuntimeError(f"can't run blastp using {queryfasta}")
    return pd.read_csv(out, header=0, sep="\t", names=header)


def find_best(
    blast_df: pd.DataFrame, df: pd.DataFrame, nevalues: int = 2, evalue: str = EVALUE
) -> pd.DataFrame:
    r = blast_df.groupby(QUERY)[evalue].nsmallest(nevalues).reset_index(level=0)
    myrdf = blast_df.loc[r.index].sort_values([QUERY, evalue], ascending=[True, True])
    myrdf = pd.merge(myrdf, df, left_on=QUERY, right_on="id")
    return myrdf


def has_blast_db(blastdir: Path, plant: str) -> bool:
    return len(glob.glob(str(blastdir / (plant + ".*")))) > 0


def has_fasta(blastdir: Path, filename: str) -> bool:
    return (blastdir / filename).exists()


def remove_files(files: list[str]) -> None:
    for f in files:
        try:
            os.remove(f)
        except OSError:
            pass


def fetch_seq(seqids: Sequence[str], plant: str, release: int) -> Iterator[SeqRecord]:
    assert blastdbcmd is not None
    u = uuid4()
    seqfile = f"{u}.seq"
    out = f"{u}.fasta"
    blastdir = blast_dir(release)
    with open(seqfile, "wt", encoding="utf-8") as fp:
        for seq in seqids:
            print(seq, file=fp)
    try:
        r = subprocess.run(
            [blastdbcmd, "-db", blastdir / plant, "-out", out, "-entry_batch", seqfile],
            check=False,
        )
        if r.returncode:
            raise RuntimeError("can't fetch sequences")
        with open(out, "rt", encoding="utf-8") as fp:
            yield from SeqIO.parse(fp, "fasta")
    finally:
        remove_files([seqfile, out])


def fetch_seq_df(seqids: Sequence[str], plant: str, release: int) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"saccver": rec.id, "subject_seq": str(rec.seq)}
            for rec in fetch_seq(seqids, plant, release)
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
        if filename is None:
            continue
        if has_fasta(blastdir, filename):
            continue
        print("fetching", filename)
        r = fetch_fasta(plant, filename, release)
        if r.returncode:
            print(r)

    for plant, filename in plants:
        if has_blast_db(blastdir, plant):
            continue
        print("creating blast db for", plant)
        r2 = mkblast(plant, filename, release)
        if r2.returncode:
            print("failed to create db")


def merge_fasta(fasta1: str, fasta2: str, out: str) -> None:
    def rf(pth: str) -> Iterator[SeqRecord]:
        with gzip.open(pth, "rt") as fp:
            yield from SeqIO.parse(fp, "fasta")

    df1 = fasta_to_df(fasta1)
    prefix1 = os.path.commonprefix(df1.id.to_list())
    df2 = fasta_to_df(fasta2)
    prefix2 = os.path.commonprefix(df2.id.to_list())
    prefix = os.path.commonprefix([prefix1, prefix2])
    df = pd.merge(
        df1, df2, how="outer", left_on="seq", right_on="seq", suffixes=("_1", "_2")
    )
    df = df.fillna("x")

    n = len(prefix)

    def nameit(s):
        i, j = s.id_1, s.id_2
        if i == "x":
            return f"{prefix}-{i}-{j[n:]}"

        if j == "x":
            return f"{prefix}-{i[n:]}-{j}"

        return f"{prefix}-{i[n:]}-{j[n:]}"

    df["id"] = df[["id_1", "id_2"]].apply(nameit, axis=1)

    write_fasta(df, out)


def write_fasta(df: pd.DataFrame, filename: str) -> None:
    r = df[["id", "seq"]]
    if filename.endswith(".gz"):
        with gzip.open(filename, "wt", encoding="utf-8") as fp:
            for _, d in r.iterrows():
                rec = SeqRecord(id=d.id, seq=Seq(d.seq), description="")
                # print(rec.seq)
                SeqIO.write(rec, fp, format="fasta")
    else:
        with open(filename, "wt", encoding="utf-8") as fp:
            for _, d in r.iterrows():
                rec = SeqRecord(id=d.id, seq=Seq(d.seq), description="")
                # print(rec.seq)
                SeqIO.write(rec, fp, format="fasta")


def blastall(
    query: str, species: Sequence[str], release: int, best: int, with_seq: bool
) -> pd.DataFrame:
    df = read_fastax(query)
    res = []

    build(species, release)
    for plant in species:
        rdf = doblast(query, plant, release=release)

        if with_seq and "saccver" in rdf.columns:
            saccver = list(rdf["saccver"])
            sdf = fetch_seq_df(saccver, plant, release)
            rdf = pd.merge(rdf, sdf, left_on="saccver", right_on="saccver")
        renames = {k: v for k, v in RENAMES.items() if k in rdf.columns}
        if renames:
            rdf = rdf.rename(columns=renames)
        myrdf = find_best(rdf, df, nevalues=best)
        myrdf["species"] = plant
        res.append(myrdf)
    ddf = pd.concat(res, axis=0)

    return ddf
