from __future__ import annotations

from typing import Iterator, Sequence

from functools import cache
import gzip
import os

import subprocess
from shutil import which
from pathlib import Path
from uuid import uuid4

import click
import pandas as pd  # type: ignore

from Bio import SeqIO  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio.Seq import Seq  # type: ignore


@cache
def safe_which(cmd: str) -> str:
    r = which(cmd)
    if r is None:
        raise click.Abort(f"can't find application: {cmd}")
    return r


def read_fasta(path: str) -> Iterator[SeqRecord]:
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fp:
            yield from SeqIO.parse(fp, "fasta")
    else:
        with open(path, "rt", encoding="utf-8") as fp:
            yield from SeqIO.parse(fp, "fasta")


def fasta_to_df(path: str, with_description: bool = False) -> pd.DataFrame:
    def todict1(rec):
        return dict(id=rec.id, seq=str(rec.seq).upper())

    def todict2(rec):
        return dict(id=rec.id, seq=str(rec.seq).upper(), description=rec.description)

    todict = todict2 if with_description else todict1
    return pd.DataFrame([todict(rec) for rec in read_fasta(path)])


class BlastDb:
    def __init__(self, database: str | Path):
        self.database = Path(database)

    def run(self, fastafile: str | Path) -> bool:
        makeblastdb = safe_which("makeblastdb")
        cwd = self.database.parent
        out = self.database.name
        fastafile = Path(fastafile)
        title = fastafile.name

        if fastafile.name.endswith(".gz"):
            cmd = [safe_which("gunzip"), "--stdout", str(fastafile)]
        else:
            cmd = [safe_which("cat"), str(fastafile)]
        # -parse_seqids so that we can get sequences out from ids with blastdbcmd
        with subprocess.Popen(
            cmd, stdout=subprocess.PIPE, cwd=str(fastafile.parent)
        ) as p1:
            with subprocess.Popen(
                [
                    makeblastdb,
                    "-in",
                    "-",
                    "-out",
                    out,
                    "-input_type=fasta",
                    "-dbtype",
                    "prot",
                    "-title",
                    title,
                    "-parse_seqids",
                    "-blastdb_version",  # https://www.biostars.org/p/390220/
                    "5",
                ],
                stdin=p1.stdout,
                cwd=str(cwd),
            ) as p2:
                if p1.stdout:
                    p1.stdout.close()
                p2.wait()
                p1.wait()

                return p2.returncode == 0


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


class Blast6:
    def __init__(self, header: Sequence[str] = HEADER):
        self.header = header

    def run(
        self,
        queryfasta: str,
        blastdb: str,
    ) -> pd.DataFrame:
        blastp = safe_which("blastp")
        outfmt = f'6 {" ".join(self.header)}'
        out = f"{uuid4()}.tsv"
        try:
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
            return out6_to_df(out, self.header)
        finally:
            remove_files([out])


def doblast6(
    queryfasta: str, blastdb: str, header: Sequence[str] = HEADER
) -> pd.DataFrame:
    b6 = Blast6(header)
    return b6.run(queryfasta, blastdb)


def out6_to_df(tsvfile: str, header=HEADER) -> pd.DataFrame:
    return pd.read_csv(tsvfile, header=0, sep="\t", names=header)


def write_fasta(df: pd.DataFrame, filename: str) -> None:
    r = df[["id", "seq"]]
    if filename.endswith(".gz"):
        with gzip.open(filename, "wt", encoding="utf-8") as fp:
            for row in r.itertuples():
                rec = SeqRecord(id=row.id, seq=Seq(row.seq), description="")
                # print(rec.seq)
                SeqIO.write(rec, fp, format="fasta")
    else:
        with open(filename, "wt", encoding="utf-8") as fp:
            for row in r.itertuples():
                rec = SeqRecord(id=row.id, seq=Seq(row.seq), description="")
                # print(rec.seq)
                SeqIO.write(rec, fp, format="fasta")


def remove_files(files: list[str]) -> None:
    for f in files:
        try:
            os.remove(f)
        except OSError:
            pass


def fetch_seq(seqids: Sequence[str], blastdb: str) -> Iterator[SeqRecord]:
    blastdbcmd = safe_which("blastdbcmd")
    u = uuid4()
    seqfile = f"{u}.seq"
    out = f"{u}.fasta"
    with open(seqfile, "wt", encoding="utf-8") as fp:
        for seq in seqids:
            print(seq, file=fp)
    try:
        r = subprocess.run(
            [blastdbcmd, "-db", blastdb, "-out", out, "-entry_batch", seqfile],
            check=False,
        )
        if r.returncode:
            raise RuntimeError("can't fetch sequences")
        with open(out, "rt", encoding="utf-8") as fp:
            yield from SeqIO.parse(fp, "fasta")
    finally:
        remove_files([seqfile, out])


def merge_fasta(fasta1: str, fasta2: str, out: str) -> None:
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


def find_best(
    blast_df: pd.DataFrame,
    query_df: pd.DataFrame,
    nevalues: int = 2,
    evalue: str = EVALUE,
) -> pd.DataFrame:
    if nevalues > 0:
        r = blast_df.groupby(QUERY)[evalue].nsmallest(nevalues).reset_index(level=0)
        myrdf = blast_df.loc[r.index].sort_values(
            [QUERY, evalue], ascending=[True, True]
        )
    else:
        myrdf = blast_df.sort_values([QUERY, evalue], ascending=[True, True])
    myrdf = pd.merge(myrdf, query_df, left_on=QUERY, right_on="id")
    return myrdf


def fetch_seq_df(seqids: Sequence[str], database: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"saccver": rec.id, "subject_seq": str(rec.seq)}
            for rec in fetch_seq(seqids, database)
        ]
    )


OKEXT = {
    "csv",
    "xlsx",
    "feather",
    "parquet",
    "pkl",
    "pickle",
    "hdf",
    "h5",
}


def get_ext(filename: Path) -> str | None:
    name = filename.name
    if name.endswith(".gz"):
        name, _ = name.rsplit(".", 1)

    if "." in name:
        ext = name.rsplit(".", 1)[-1]
    else:
        ext = None
    return ext


def save_df(
    df: pd.DataFrame,
    filename: str | Path,
    index: bool = False,
    default: str = "csv",
    key: str = "blast",
) -> None:
    filename = Path(filename)
    ext = get_ext(filename)
    if ext is None:
        ext = default

    if ext == "csv":
        df.to_csv(filename, index=index)
    elif ext == "xlsx":
        df.to_excel(filename, index=index)
    elif ext == "feather":
        df.to_feather(filename, index=index)
    elif ext == "parquet":
        df.to_parquet(filename, index=index)
    elif ext in {"pkl", "pickle"}:
        df.to_pickle(filename)
    elif ext in {"hdf", "h5"}:
        df.to_hdf(filename, key)
    else:
        raise ValueError(f"unknown file extension: {filename}")


def toblastdb(blastdbs: Sequence[str]) -> list[str]:
    s = {b.rsplit(".", 1)[0] for b in blastdbs}
    return sorted(s)


def blastall(
    query: str, blastdbs: Sequence[str], best: int, with_seq: bool
) -> pd.DataFrame:
    df = fasta_to_df(query)
    res = []

    b6 = Blast6()
    for blastdb in blastdbs:
        rdf = b6.run(query, blastdb)

        if with_seq and "saccver" in rdf.columns:
            saccver = list(rdf["saccver"])
            sdf = fetch_seq_df(saccver, blastdb)
            rdf = pd.merge(rdf, sdf, left_on="saccver", right_on="saccver")

        myrdf = find_best(rdf, df, nevalues=best)
        myrdf["blastdb"] = Path(blastdb).name
        res.append(myrdf)
    ddf = pd.concat(res, axis=0)

    return ddf


def buildall(fastafiles: Sequence[str], builddir: str | Path | None = None) -> None:
    if builddir is not None:
        builddir = Path(builddir)
        if not builddir.exists():
            builddir.mkdir(exist_ok=True, parents=True)
    for fastafile in fastafiles:
        fa = Path(fastafile)
        name = fa.name
        name, _ = name.split(".", 1)
        name = name.lower()
        db = builddir / name if builddir else fa.parent / name
        b = BlastDb(db)
        ok = b.run(fa)
        if not ok:
            raise RuntimeError(f"can't build database with {fastafile}")
