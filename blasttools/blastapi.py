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
import pandas as pd
from pandas.errors import UndefinedVariableError

from Bio import SeqIO  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from Bio.Seq import Seq  # type: ignore


@cache
def safe_which(cmd: str) -> str:
    r = which(cmd)
    if r is None:
        raise click.ClickException(f'can\'t find application: "{cmd}"')
    return r


def read_fasta(path: str) -> Iterator[SeqRecord]:
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fp:
            yield from SeqIO.parse(fp, "fasta")
    else:
        with open(path, "rt", encoding="utf-8") as fp:
            yield from SeqIO.parse(fp, "fasta")


def has_pdatabase(path: str) -> bool:
    if "." not in path:
        path += ".pdb"
    return Path(path).exists()


def fasta_to_df(path: str, with_description: bool = False) -> pd.DataFrame:
    def todict1(rec: SeqRecord) -> dict[str, str]:
        return dict(id=rec.id, seq=str(rec.seq).upper())

    def todict2(rec: SeqRecord) -> dict[str, str]:
        return dict(id=rec.id, seq=str(rec.seq).upper(), description=rec.description)

    todict = todict2 if with_description else todict1
    return pd.DataFrame([todict(rec) for rec in read_fasta(path)])


class BlastDb:
    def __init__(self, database: str | Path, *, blastp: bool = True):
        self.database = Path(database)
        self.blastp = blastp

    def run(self, fastafile: str | Path) -> bool:
        makeblastdb = safe_which("makeblastdb")
        cwd = self.database.parent
        out = self.database.name
        fastafile = Path(fastafile)
        title = fastafile.name

        if fastafile.name.endswith(".gz"):
            cmd = [safe_which("gunzip"), "--stdout", fastafile.name]
        else:
            cmd = [safe_which("cat"), fastafile.name]
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
                    "prot" if self.blastp else "nucl",
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
# with added 'gaps'
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


def mkheader(header: str) -> Sequence[str]:
    from .columns import VALID

    add = header.startswith("+")
    sub = header.startswith("-")
    h = header[1:] if add or sub else header
    hl = h.strip().split()
    if sub:
        hl = [h for h in HEADER if h not in hl]
    elif add:
        hl = list(HEADER) + [h for h in hl if h not in HEADER]
    unknown = set(hl) - set(VALID)
    if unknown:
        raise click.ClickException(f"unknown headers \"{' '.join(unknown)}\"")
    return tuple(hl)


class Blast6:
    def __init__(
        self,
        header: Sequence[str] | None = None,
        num_threads: int = 1,
        blastp: bool = True,
    ):
        if header is None:
            header = HEADER
        self.header = header
        self.num_threads = num_threads
        self.blastp = blastp

    def get_blast(self) -> str:
        return safe_which("blastp") if self.blastp else safe_which("blastn")

    def run(
        self,
        queryfasta: str,
        blastdb: str,
    ) -> pd.DataFrame:
        blast = self.get_blast()
        outfmt = f'6 {" ".join(self.header)}'
        out = f"{uuid4()}.tsv"
        try:
            r = subprocess.run(
                [
                    blast,
                    "-outfmt",
                    outfmt,
                    "-query",
                    queryfasta,
                    "-db",
                    blastdb,
                    "-out",
                    out,
                    "-num_threads",
                    str(self.num_threads),
                ],
                check=False,
            )
            if r.returncode:
                b = "blastp" if self.blastp else "blastn"
                raise click.ClickException(f"can't run {b} using {queryfasta}")
            return out6_to_df(out, self.header)
        finally:
            remove_files([out])


def doblast6(
    queryfasta: str,
    blastdb: str,
    header: Sequence[str] | None = None,
    *,
    num_threads: int = 1,
) -> pd.DataFrame:
    b6 = Blast6(header, num_threads=num_threads)
    return b6.run(queryfasta, blastdb)


def out6_to_df(tsvfile: str, header: Sequence[str] = HEADER) -> pd.DataFrame:
    return pd.read_csv(tsvfile, header=0, sep="\t", names=list(header))


def write_fasta(df: pd.DataFrame, filename: str) -> None:
    if "description" in df.columns:
        r = df[["id", "seq", "description"]]
        wd = True
    else:
        r = df[["id", "seq"]]
        wd = False

    def toseq(row):
        d = row.description if wd else ""
        return SeqRecord(id=row.id, seq=Seq(row.seq), description=d)

    if filename.endswith(".gz"):
        with gzip.open(filename, "wt", encoding="utf-8") as fp:
            for row in r.itertuples():
                SeqIO.write(toseq(row), fp, format="fasta")
    else:
        with open(filename, "wt", encoding="utf-8") as fp:
            for row in r.itertuples():
                SeqIO.write(toseq(row), fp, format="fasta")


def remove_files(files: list[str | Path]) -> None:
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
    try:
        with open(seqfile, "wt", encoding="utf-8") as fp:
            for seq in seqids:
                print(seq, file=fp)
        r = subprocess.run(
            [blastdbcmd, "-db", blastdb, "-out", out, "-entry_batch", seqfile],
            check=False,
        )
        if r.returncode:
            raise click.ClickException("can't fetch sequences")
        with open(out, "rt", encoding="utf-8") as fp:
            yield from SeqIO.parse(fp, "fasta")
    finally:
        remove_files([seqfile, out])


def merge_fasta(
    fasta1: str, fasta2: str, out: str, with_description: bool = True
) -> None:
    df1 = fasta_to_df(fasta1, with_description=with_description)
    df2 = fasta_to_df(fasta2, with_description=with_description)
    prefix1 = os.path.commonprefix(df1.id.to_list())
    prefix2 = os.path.commonprefix(df2.id.to_list())
    prefix = os.path.commonprefix([prefix1, prefix2])

    FILLNA = "x"
    df = pd.merge(
        df1, df2, how="outer", left_on="seq", right_on="seq", suffixes=("_1", "_2")
    )
    df = df.fillna(FILLNA)

    n = len(prefix)

    def nameit(s):
        i, j = s.id_1, s.id_2
        if i == FILLNA:
            return f"{prefix}-{i}-{j[n:]}"

        if j == FILLNA:
            return f"{prefix}-{i[n:]}-{j}"

        return f"{prefix}-{i[n:]}-{j[n:]}"

    df["id"] = df[["id_1", "id_2"]].apply(nameit, axis=1)

    def desc(s):
        i, j = s.description_1, s.description_2
        if i == j:
            return i
        if i == FILLNA:
            return j
        if j == FILLNA:
            return i
        return f"{i} | {j}"

    df["description"] = df[["description_1", "description_2"]].apply(desc, axis=1)
    write_fasta(df, out)


EVAL_COL = "__xxxx__"


def find_best(
    blast_df: pd.DataFrame,
    query_df: pd.DataFrame,
    nevalues: int = 2,
    evalue_col: str = EVALUE,
    query_col: str = QUERY,
) -> pd.DataFrame:
    if evalue_col not in blast_df:
        blast_df[EVAL_COL] = blast_df.eval(query_col)

    if nevalues > 0:
        r = (
            blast_df.groupby(query_col)[evalue_col]
            .nsmallest(nevalues)
            .reset_index(level=0)
        )
        myrdf = blast_df.loc[r.index].sort_values(
            [query_col, evalue_col], ascending=[True, True]
        )
    else:
        myrdf = blast_df.sort_values([query_col, evalue_col], ascending=[True, True])
    myrdf = pd.merge(myrdf, query_df, left_on=query_col, right_on="id")
    todrop = ["id"]
    if EVAL_COL in myrdf:
        todrop.append(EVAL_COL)
    myrdf.drop(columns=todrop, inplace=True)
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


def check_ext(filename: str) -> None:
    ext = get_ext(Path(filename))
    if ext is not None and ext in OKEXT:
        return
    ex = ",".join(OKEXT)
    raise click.ClickException(
        f'unknown output type for file: "{filename}": require extension to be one of .{{{ex}}}'
    )


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
    try:
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
            click.secho(
                f'unknown file extension for "{filename}", saving as csv',
                fg="red",
                bold=True,
                err=True,
            )
            df.to_csv(filename, index=index)
    except ModuleNotFoundError as exc:
        csvf = str(filename) + ".csv"
        click.secho(
            f"Can't save as {ext} ({exc}). Will save as *CSV* to this {csvf}",
            err=True,
            fg="red",
            bold=True,
        )
        df.to_csv(csvf, index=False)


def read_df(
    filename: str | Path,
    key: str = "blast",
) -> pd.DataFrame | None:
    filename = Path(filename)
    ext = get_ext(filename)
    if ext is None:
        return None

    elif ext == "xlsx":
        return pd.read_excel(filename)
    elif ext == "feather":
        return pd.read_feather(filename)
    elif ext == "parquet":
        return pd.read_parquet(filename)
    elif ext in {"pkl", "pickle"}:
        return pd.read_pickle(filename)
    elif ext in {"hdf", "h5"}:
        return pd.read_hdf(filename, key)  # type: ignore

    # if ext == "csv":
    return pd.read_csv(filename)


def toblastdb(blastdbs: Sequence[str]) -> list[str]:
    s = {b.rsplit(".", maxsplit=1)[0] for b in blastdbs}
    return sorted(s)


def check_expr(headers: Sequence[str], expr: str):
    df = pd.DataFrame({col: [] for col in headers})
    try:
        df.eval(expr)
    except UndefinedVariableError as exc:
        raise click.BadParameter(
            f"expression supplied references unknown column(s): {exc}",
            param_hint="expr",
        ) from exc
    except SyntaxError as exc:
        raise click.BadParameter(str(exc), param_hint="expr") from exc


def blastall(
    queryfasta: str,
    blastdbs: Sequence[str],
    best: int,
    with_seq: bool,
    *,
    num_threads: int = 1,
    blastp: bool = True,
    header: Sequence[str] | None = None,
    with_description: bool = False,
    expr: str = EVALUE,
) -> pd.DataFrame:
    df = fasta_to_df(queryfasta, with_description=with_description)
    if not df["id"].is_unique:
        raise click.ClickException(
            f'sequences IDs are not unique for query file "{queryfasta}"'
        )
    res = []

    b6 = Blast6(header, num_threads=num_threads, blastp=blastp)

    h = b6.header
    check_expr(h, expr)
    for blastdb in blastdbs:
        rdf = b6.run(queryfasta, blastdb)

        if with_seq and "saccver" in rdf.columns:
            saccver = list(rdf["saccver"])
            sdf = fetch_seq_df(saccver, blastdb)
            rdf = pd.merge(rdf, sdf, left_on="saccver", right_on="saccver")

        myrdf = find_best(rdf, df, nevalues=best, evalue_col=expr)
        myrdf["blastdb"] = Path(blastdb).name
        res.append(myrdf)
    ddf = pd.concat(res, axis=0)

    return ddf


def concat_fasta(fastafiles: Sequence[str], out: str) -> None:
    if out.endswith(".gz"):
        of = gzip.open(out, "wt")
    else:
        of = open(out, "wt", encoding="utf-8")
    with of:
        for fa in fastafiles:
            for rec in read_fasta(fa):
                SeqIO.write(rec, of, "fasta")


def buildall(
    fastafiles: Sequence[str],
    builddir: str | Path | None = None,
    merge: str | None = None,
    *,
    blastp: bool = True,
) -> None:
    if builddir is not None:
        builddir = Path(builddir)
        if not builddir.exists():
            builddir.mkdir(exist_ok=True, parents=True)
    out = None
    if merge:
        out = f"{merge}.{uuid4()}.fa.gz"
        concat_fasta(fastafiles, out)
        fastafiles = [out]
    try:
        for fastafile in fastafiles:
            fa = Path(fastafile)
            name = fa.name
            name, _ = name.split(".", maxsplit=1)
            name = name.lower()
            db = builddir / name if builddir else fa.parent / name
            b = BlastDb(db, blastp=blastp)
            ok = b.run(fa)
            if not ok:
                raise click.ClickException(f"can't build database with {fastafile}")
    finally:
        if out:
            remove_files([out])
