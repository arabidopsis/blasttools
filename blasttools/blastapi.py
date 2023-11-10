from __future__ import annotations

import gzip
import os
import subprocess
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from shutil import which
from typing import Iterator
from typing import Sequence
from uuid import uuid4

import click
import pandas as pd
from Bio import SeqIO  # type: ignore
from Bio.Seq import Seq  # type: ignore
from Bio.SeqRecord import SeqRecord  # type: ignore
from pandas.errors import UndefinedVariableError


@cache
def safe_which(cmd: str) -> str:
    r = which(cmd)
    if r is None:
        raise click.ClickException(f'can\'t find application: "{cmd}"')
    return r


def read_fasta(path: str | Path) -> Iterator[SeqRecord]:
    path = Path(path)
    if path.name.endswith(".gz"):
        with gzip.open(path, "rt") as fp:
            yield from SeqIO.parse(fp, "fasta")
    else:
        with open(path, encoding="utf-8") as fp:
            yield from SeqIO.parse(fp, "fasta")


def has_pdatabase(path: str) -> bool:
    if "." not in path:
        path += ".pdb"
    return Path(path).exists()


def fasta_to_df(path: str | Path, with_description: bool = False) -> pd.DataFrame:
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
            cmd,
            stdout=subprocess.PIPE,
            cwd=str(fastafile.parent),
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
        queryfasta: str | Path,
        blastdb: str | Path,
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
                    str(queryfasta),
                    "-db",
                    str(blastdb),
                    "-out",
                    out,
                    "-num_threads",
                    str(self.num_threads),
                ],
                check=False,
            )
            if r.returncode:
                b = "blastp" if self.blastp else "blastn"
                raise click.ClickException(f"Can't run {b} using {queryfasta}")
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

    def toseq(row) -> SeqRecord:
        d = row.description if wd else ""
        return SeqRecord(id=row.id, seq=Seq(row.seq), description=d)

    write_fasta_iter(Path(filename), (toseq(row) for row in r.itertuples()))


def write_fasta_iter(filename: Path, records: Iterator[SeqRecord]) -> None:
    if filename.name.endswith(".gz"):
        with gzip.open(filename, "wt", encoding="utf-8") as fp:
            for rec in records:
                SeqIO.write(rec, fp, format="fasta")
    else:
        with open(filename, "w", encoding="utf-8") as fp:
            for rec in records:
                SeqIO.write(rec, fp, format="fasta")


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
        with open(seqfile, "w", encoding="utf-8") as fp:
            for seq in seqids:
                print(seq, file=fp)
        r = subprocess.run(
            [blastdbcmd, "-db", blastdb, "-out", out, "-entry_batch", seqfile],
            check=False,
        )
        if r.returncode:
            raise click.ClickException("Can't fetch sequences")
        with open(out, encoding="utf-8") as fp:
            yield from SeqIO.parse(fp, "fasta")
    finally:
        remove_files([seqfile, out])


def merge_fasta(
    fasta1: str,
    fasta2: str,
    out: str,
    with_description: bool = True,
) -> None:
    df1 = fasta_to_df(fasta1, with_description=with_description)
    df2 = fasta_to_df(fasta2, with_description=with_description)
    prefix1 = os.path.commonprefix(df1.id.to_list())
    prefix2 = os.path.commonprefix(df2.id.to_list())
    prefix = os.path.commonprefix([prefix1, prefix2])

    FILLNA = "x"
    df = pd.merge(
        df1,
        df2,
        how="outer",
        left_on="seq",
        right_on="seq",
        suffixes=("_1", "_2"),
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


TMP_EVAL_COL = "__xxxx__"


def find_best(
    blast_df: pd.DataFrame,
    query_df: pd.DataFrame,
    nevalues: int = 2,
    evalue_col: str = EVALUE,  # or qstart - qend + mismatch
    id_col: str = QUERY,
    query_id_col: str = "id",
) -> pd.DataFrame:
    if evalue_col not in blast_df:
        blast_df[TMP_EVAL_COL] = blast_df.eval(
            evalue_col,
        )  # expression like 'qstart - qend'
        evalue_col = TMP_EVAL_COL
    if nevalues > 0:
        r = (
            blast_df.groupby(id_col)[evalue_col]
            .nsmallest(nevalues)
            .reset_index(level=0)
        )
        myrdf = blast_df.loc[r.index].sort_values(
            [id_col, evalue_col],
            ascending=[True, True],
        )
    else:
        myrdf = blast_df.sort_values([id_col, evalue_col], ascending=[True, True])
    myrdf = pd.merge(myrdf, query_df, left_on=id_col, right_on=query_id_col)
    todrop = [query_id_col]
    if TMP_EVAL_COL in myrdf:
        todrop.append(TMP_EVAL_COL)
    myrdf.drop(columns=todrop, inplace=True)
    return myrdf


def fetch_seq_df(seqids: Sequence[str], database: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"saccver": rec.id, "subject_seq": str(rec.seq)}
            for rec in fetch_seq(seqids, database)
        ],
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
        f'unknown output type for file: "{filename}": require extension to be one of .{{{ex}}}',
    )


def try_save(
    ext: str,
    df: pd.DataFrame,
    filename: str | Path,
    index: bool = False,
    key: str = "blast",
) -> bool:
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
        return False
    return True


def test_save(filename: str | Path):
    import tempfile

    ext = get_ext(Path(filename))
    if ext is None:
        return  # CSV
    df = pd.DataFrame(dict(x=[1]))
    with tempfile.NamedTemporaryFile(suffix="." + ext) as fp:
        try:
            try_save(ext, df, fp.name)
        except ModuleNotFoundError as exc:
            raise click.ClickException(
                f"Can't save DataFrame as {filename} ({exc})",
            ) from exc


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
        ok = try_save(ext, df, filename, index, key)
        if ok:
            return

        click.secho(
            f'unknown file extension for "{filename}", saving as csv',
            fg="red",
            bold=True,
            err=True,
        )
        df.to_csv(filename, index=index)
    except ModuleNotFoundError as exc:
        csvf = str(filename) + ".csv"
        msg1 = click.style(
            f"Can't save as file type: {ext} ({exc}).",
            fg="red",
            bold=True,
        )
        msg2 = click.style(f'Will save as *CSV* to "{csvf}".', fg="green", bold=True)
        click.echo(
            f"{msg1} {msg2}",
            err=True,
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


def check_expr(headers: Sequence[str], expr: str) -> None:
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


@dataclass
class BlastConfig:
    best: int = 0
    with_seq: bool = False
    header: Sequence[str] | None = None
    # path: str | None = None
    num_threads: int = 1
    with_description: bool = True
    expr: str = EVALUE
    blastp: bool = True
    without_query_seq: bool = False


def blastall(
    queryfasta: str,
    blastdbs: Sequence[str],
    *,
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
    for blastdb in blastdbs:
        rdf = b6.run(queryfasta, blastdb)

        if config.with_seq and "saccver" in rdf.columns:
            saccver = list(rdf["saccver"])
            sdf = fetch_seq_df(saccver, blastdb)
            rdf = pd.merge(rdf, sdf, left_on="saccver", right_on="saccver")

        myrdf = find_best(rdf, qdf, nevalues=config.best, evalue_col=config.expr)
        myrdf["blastdb"] = Path(blastdb).name
        res.append(myrdf)
    ddf = pd.concat(res, axis=0, ignore_index=True)

    return ddf


def concat_fasta(fastafiles: Sequence[str], out: str) -> None:
    if out.endswith(".gz"):
        of = gzip.open(out, "wt")
    else:
        of = open(out, "w", encoding="utf-8")
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
                raise click.ClickException(f"Can't build database with {fastafile}")
    finally:
        if out:
            remove_files([out])
