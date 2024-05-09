from __future__ import annotations

import gzip
import os
import subprocess
import sys
from dataclasses import dataclass
from functools import cache
from pathlib import Path
from shutil import which
from typing import Any
from typing import Iterator
from typing import Sequence
from uuid import uuid4

import click
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pandas.errors import UndefinedVariableError


@cache
def safe_which(cmd: str) -> str:
    r = which(cmd)
    if r is None:
        raise click.ClickException(f'can\'t find application: "{cmd}" in PATH')
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


def get_cat(filename: Path) -> list[str]:
    if filename.name.endswith(".gz"):
        return [safe_which("gunzip"), "--stdout"]
    if filename.name.endswith(".bz2"):
        return [safe_which("bzcat")]
    return [safe_which("cat")]


def fasta_to_df(
    path: str | Path,
    with_description: bool = False,
    without_query_seq: bool = False,
) -> pd.DataFrame:
    """Read a Fasta file into a pandas dataframe"""

    def todict1(rec: SeqRecord) -> dict[str, str]:
        return dict(id=rec.id or "", seq=str(rec.seq).upper())

    def todict2(rec: SeqRecord) -> dict[str, str]:
        return dict(
            id=rec.id or "",
            seq=str(rec.seq).upper(),
            description=rec.description,
        )

    todict = todict2 if with_description else todict1
    qdf = pd.DataFrame([todict(rec) for rec in read_fasta(path)])

    if not qdf["id"].is_unique:
        raise click.ClickException(
            f'sequences IDs are not unique for query file "{path}"',
        )
    if without_query_seq:
        qdf.drop(columns=["seq"], inplace=True)
    return qdf


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
        cat_cmd = get_cat(fastafile)
        cat_cmd.append(fastafile.name)

        # -parse_seqids so that we can get sequences out from ids with blastdbcmd
        with subprocess.Popen(
            cat_cmd,
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

    hl = []
    for v in h.strip().split():
        if "," in v:
            hl.extend([s.strip() for s in v.split(",")])
        else:
            hl.append(v)
    if sub:
        hl = [h for h in HEADER if h not in hl]
    elif add:
        _hl = list(HEADER)
        for h in hl:
            if h not in _hl:  # remove duplicates
                _hl.append(h)
        hl = _hl
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

    def blastall(
        self,
        qdf: pd.DataFrame,
        queryfasta: str | Path,
        blastdbs: Sequence[tuple[str, str | Path]],
        *,
        config,
    ) -> pd.DataFrame:
        res = []
        for species, blastdb in blastdbs:
            rdf = self.run(
                queryfasta,
                blastdb,
                needs_translation=config.needs_translation,
            )

            if config.with_subject_seq and "saccver" in rdf.columns:
                saccver = list(set(rdf["saccver"]))
                sdf = fetch_seq_df(saccver, blastdb)
                rdf = pd.merge(rdf, sdf, left_on="saccver", right_on="saccver")

            myrdf = find_bestx(rdf, qdf, nbest=config.best, expr=config.expr)
            myrdf["species"] = species
            res.append(myrdf)
        ddf = pd.concat(res, axis=0, ignore_index=True)

        return ddf

    def get_blast_args(self) -> list[str]:
        if "BLAST_ARGS" not in os.environ:
            return []
        return os.environ["BLAST_ARGS"].split()

    def get_output(self) -> str:
        return f"{uuid4()}.tsv"

    def run(
        self,
        queryfasta: str | Path,
        blastdb: str | Path,
        *,
        needs_translation: bool = False,
    ) -> pd.DataFrame:
        blast = self.get_blast()
        outfmt = f'6 {" ".join(self.header)}'
        out = self.get_output()
        queryfasta = Path(queryfasta)
        if not needs_translation:
            cat_cmd = get_cat(queryfasta)
            cat_cmd.append(queryfasta.name)
        else:
            cat_cmd = [sys.executable, "-m", "blasttools.translate", queryfasta.name]
        try:
            with subprocess.Popen(
                cat_cmd,
                stdout=subprocess.PIPE,
                cwd=str(queryfasta.parent),
            ) as p1:
                with subprocess.Popen(
                    [
                        blast,
                        "-outfmt",
                        outfmt,
                        "-query",
                        "-",
                        "-db",
                        str(blastdb),
                        "-out",
                        out,
                        "-num_threads",
                        str(self.num_threads),
                        *self.get_blast_args(),
                    ],
                    stdin=p1.stdout,
                ) as p2:
                    if p1.stdout:
                        p1.stdout.close()
                    p2.wait()
                    p1.wait()

                    if p2.returncode:
                        b = "blastp" if self.blastp else "blastn"
                        raise click.ClickException(f"Can't run {b} using {queryfasta}")
                    return out6_to_df(out, self.header)
        finally:
            remove_files([out])


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


def fetch_seq(seqids: Sequence[str], blastdb: str | Path) -> Iterator[SeqRecord]:
    blastdbcmd = safe_which("blastdbcmd")
    u = uuid4()
    seqfile = f"{u}.seq"
    out = f"{u}.fasta"
    try:
        with open(seqfile, "w", encoding="utf-8") as fp:
            for seq in seqids:
                print(seq, file=fp)
        r = subprocess.run(
            [blastdbcmd, "-db", str(blastdb), "-out", out, "-entry_batch", seqfile],
            check=False,
        )
        if r.returncode:
            raise click.ClickException("Can't fetch sequences")
        with open(out, encoding="utf-8") as fp:
            yield from SeqIO.parse(fp, "fasta")
    finally:
        remove_files([seqfile, out])


def fasta_xref(fasta1: str, fasta2: str, fillna: str = "_") -> pd.DataFrame:
    df1 = fasta_to_df(fasta1, with_description=False)
    df2 = fasta_to_df(fasta2, with_description=False)

    df = pd.merge(
        df1,
        df2,
        how="outer",
        left_on="seq",
        right_on="seq",
        suffixes=("", "_1"),
    )
    df = df.fillna(fillna)
    df = df[["id", "id_1"]]
    df.rename(columns={"id_1": "other"}, inplace=True)
    return df


def fasta_merge(
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
    nbest: int = 2,
    expr: str = EVALUE,  # or qstart - qend + mismatch
    group_by_col: str = QUERY,
) -> pd.DataFrame:
    if expr not in blast_df:
        blast_df[TMP_EVAL_COL] = blast_df.eval(
            expr,
        )  # expression like 'qstart - qend'
        expr = TMP_EVAL_COL
    if blast_df[expr].dtype.kind == "b":
        # if its a boolean expression the just filter
        myrdf = blast_df[blast_df[expr]]
    else:
        if nbest > 0:
            r = (
                blast_df[[group_by_col, expr]]
                .groupby(group_by_col)[expr]
                .nsmallest(nbest)
                .reset_index(level=0)
            )
            myrdf = blast_df.loc[r.index]
        else:
            myrdf = blast_df

    myrdf = myrdf.sort_values(
        [group_by_col, expr],
        ascending=[True, True],
    )

    if TMP_EVAL_COL in myrdf.columns:
        myrdf.drop(columns=[TMP_EVAL_COL], inplace=True)
    return myrdf


def find_bestx(
    blast_df: pd.DataFrame,
    query_df: pd.DataFrame,
    nbest: int = 2,
    expr: str = EVALUE,  # or qstart - qend + mismatch
    group_by_col: str = QUERY,
    query_id_col: str = "id",
) -> pd.DataFrame:
    myrdf = find_best(
        blast_df,
        nbest,
        expr=expr,
        group_by_col=group_by_col,
    )
    myrdf = pd.merge(myrdf, query_df, left_on=group_by_col, right_on=query_id_col)
    if group_by_col != query_id_col:
        myrdf.drop(columns=[query_id_col], inplace=True)
    return myrdf


def fetch_seq_df(seqids: Sequence[str], database: str | Path) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"saccver": rec.id, "subject_seq": str(rec.seq)}
            for rec in fetch_seq(seqids, database)
        ],
    )


OKEXT = {"csv", "xlsx", "feather", "parquet", "pkl", "pickle", "hdf", "h5", "hd5"}


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
    elif ext in {"hdf", "h5", "hd5"}:
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
        except (ModuleNotFoundError, ImportError) as exc:
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
    if not filename.parent.exists():
        filename.parent.mkdir(parents=True, exist_ok=True)
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

    if ext == "xlsx":
        return pd.read_excel(filename)
    elif ext == "feather":
        return pd.read_feather(filename)
    elif ext == "parquet":
        return pd.read_parquet(filename)
    elif ext in {"pkl", "pickle"}:
        return pd.read_pickle(filename)
    elif ext in {"hdf", "h5", "hd5"}:
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
    """retain only n `best` sequences according to `expr`"""
    with_subject_seq: bool = False
    """put subject sequence into resulting dataframe (as subject_seq)"""
    header: Sequence[str] | None = None
    """Blast columns to generate"""
    # path: str | None = None
    num_threads: int = 1
    with_description: bool = True
    """add query description field (from fasta)"""
    expr: str = EVALUE
    """expression or column that defines 'best' (i.e. lowest) blast hit"""
    blastp: bool = True
    without_query_seq: bool = False
    """don't retain query sequence column (as seq)"""
    xml: bool = False
    """Use blast xml output to obtain match, query, sbjct sequences"""
    needs_translation: bool = False
    """query fasta contains mixed rna/dna/rna sequences too"""


def blastall(
    queryfasta: str,
    blastdbs: Sequence[str],
    *,
    config: BlastConfig = BlastConfig(),
) -> pd.DataFrame:
    from .blastxml import BlastXML

    b6: Blast6
    if config.xml:
        b6 = BlastXML(
            config.header,
            num_threads=config.num_threads,
            blastp=config.blastp,
        )
    else:
        b6 = Blast6(config.header, num_threads=config.num_threads, blastp=config.blastp)
    if config.expr not in b6.header:
        check_expr(b6.header, config.expr)  # fail early

    qdf = fasta_to_df(
        queryfasta,
        with_description=config.with_description,
        without_query_seq=config.without_query_seq,
    )

    return b6.blastall(
        qdf,
        queryfasta,
        [(Path(db).name, db) for db in blastdbs],
        config=config,
    )


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


def list_out(it: Iterator[Any], *, use_null: bool = False) -> None:
    if not use_null:
        for path in it:
            click.echo(str(path))
    else:
        for i, path in enumerate(it):
            if i != 0:
                click.echo("\0", nl=False)
            click.echo(str(path), nl=False)
