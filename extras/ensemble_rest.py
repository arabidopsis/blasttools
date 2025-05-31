from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import cast
from typing import Iterable
from typing import Iterator
from typing import Literal
from typing import NamedTuple
from typing import overload

import click
import pandas as pd
import requests

ENSEMBLE_REST = "https://rest.ensembl.org"


@dataclass(kw_only=True)
class Replacement:
    stable_id: str
    score: float


@dataclass(kw_only=True)
class EnsembleArchive:
    release: str
    id: str
    is_current: str
    possible_replacement: list[Replacement]
    assembly: str
    peptide: str | None
    latest: str
    type: str
    version: int


@overload
def archive_ids(
    ids: Iterable[str],
    verbose: Literal[True],
) -> tuple[list[EnsembleArchive], str | list[Any]]: ...


@overload
def archive_ids(
    ids: Iterable[str],
    verbose: Literal[False],
) -> list[EnsembleArchive]: ...
@overload
def archive_ids(ids: Iterable[str]) -> list[EnsembleArchive]: ...


def archive_ids(ids: Iterable[str], verbose: bool = False):
    data = {"id": list(ids)}

    ext = "/archive/id"
    # headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    r = requests.post(ENSEMBLE_REST + ext, json=data)
    if r.status_code != 200:
        if verbose:
            return [], r.text
        return []
    ret: list[EnsembleArchive] = []
    jsn = r.json()
    for d in jsn:
        d["possible_replacement"] = [
            Replacement(**i) for i in d["possible_replacement"]
        ]
        ret.append(EnsembleArchive(**d))

    if verbose:
        return ret, jsn
    return ret


# fake type for row itertuples()
class CRow(NamedTuple):
    v11: str
    v11_confidence: str
    v21: str
    v21_confidence: str


def get_21_conversions() -> pd.DataFrame:
    iwgsc = pd.read_csv(
        "https://urgi.versailles.inrae.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v2.1/iwgsc_refseq_all_correspondances.zip",
        sep=" ",
    )

    # 02 and 6 numbers
    TRAES02 = r"^(TraesCS(?:[1234567][ABD]|U)02[G][0-9]{6})(LC)?$"
    v11a = iwgsc["v1.1"].str.extract(TRAES02, expand=True)
    v11a.columns = ["v11", "v11_confidence"]
    v11a.loc[v11a.v11_confidence.isna(), "v11_confidence"] = "HC"

    assert not v11a.isna().any().any()

    # 03 and 7 numbers
    TRAES03 = r"^(TraesCS(?:[1234567][ABD]|U)03[G][0-9]{7})(LC)?$"
    v21a = iwgsc["v2.1"].str.extract(TRAES03, expand=True)
    v21a.columns = ["v21", "v21_confidence"]
    v21a.loc[v21a.v21_confidence.isna(), "v21_confidence"] = "HC"

    assert not v21a.isna().any().any()

    convert = pd.concat([v11a, v21a], axis=1)
    convert = convert.drop_duplicates(ignore_index=True)
    return convert


class GeneFamily(NamedTuple):
    description: str
    transcripts: list[str]


def read_genes(filename: str) -> list[GeneFamily]:
    """FASTA like file

    With `>description`
    followed by space|newline separated transcript ids
    """
    current: list[str] = []
    idata: list[GeneFamily] = []
    desc = None
    with open(filename, encoding="utf8") as fp:
        for line in fp:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if desc and current:
                    idata.append(GeneFamily(desc, current))
                current = []
                desc = line[1:]
            else:
                for i in line.split():
                    if i:
                        current.append(i)
    if desc and current:
        idata.append(GeneFamily(desc, current))
    return idata


def findall(
    gene_names: list[GeneFamily],
    convert: pd.DataFrame,
    *,
    sleep: float = 0.0,
    low_confidence: bool = False,
) -> Iterator[tuple[str, ...]]:
    n = len(gene_names)
    V11 = convert["v11"].str.upper()
    HC = convert["v21_confidence"] == "HC"
    for idx, (desc, ids) in enumerate(gene_names):
        if sleep and idx != 0:
            time.sleep(sleep)
        idss = set(ids)
        click.secho(f"{desc}[{len(idss)}]: {idx + 1}/{n}", fg="yellow", nl=False)
        res = archive_ids(idss)
        click.secho(" done", fg="green")
        for r in res:
            if r.id in idss:
                idss.remove(r.id)
            elif r.latest in idss:
                idss.remove(r.latest)
            for p in r.possible_replacement:
                q = p.stable_id.upper() == V11
                if not low_confidence:
                    q = q & HC
                v2ids = convert[q][["v21", "v21_confidence"]]
                if len(v2ids) == 0:
                    yield (desc, r.id, p.stable_id, "missing", "")
                else:
                    for t in v2ids.itertuples():
                        y = cast(CRow, t)
                        yield (desc, r.id, p.stable_id, y.v21, y.v21_confidence)
        if idss:
            for iid in idss:
                yield (desc, iid, "missing", "missing", "")


def get_v21_ids(
    idfilename: str,
    *,
    sleep: float = 0.0,
    low_confidence: bool = False,
) -> pd.DataFrame:
    gene_names = read_genes(idfilename)
    convert = get_21_conversions()
    columns = [
        "description",
        "transcript_id",
        "stable_id_v1",
        "stable_id_v2",
        "confidence",
    ]
    df = pd.DataFrame(
        list(findall(gene_names, convert, sleep=sleep, low_confidence=low_confidence)),
        columns=columns,
    )

    return df


def convert(idfilename: str, low_confidence: bool = False) -> pd.DataFrame:
    gene_names = read_genes(idfilename)
    convert = get_21_conversions()
    columns = ["description", "stable_id_v1", "stable_id_v2", "confidence"]
    return pd.DataFrame(
        list(run_conversion(gene_names, convert, low_confidence=low_confidence)),
        columns=columns,
    )


def run_conversion(
    gene_names: list[GeneFamily],
    convert: pd.DataFrame,
    low_confidence: bool = False,
) -> Iterator[tuple[str, ...]]:
    hc = convert["v21_confidence"] == "HC"
    V11 = convert["v11"].str.upper()
    for desc, ids in gene_names:
        for iid in ids:
            if "." in iid:
                gene, ext = iid.split(".")
                ext = "." + ext
            else:
                gene, ext = iid, ""
            q = gene.upper() == V11
            if not low_confidence:
                q = q & hc
            v2ids = convert[q][["v21", "v21_confidence"]]
            if len(v2ids) == 0:
                yield (desc, iid, "missing", "")
            else:
                for t in v2ids.itertuples():
                    y = cast(CRow, t)
                    yield (desc, iid, y.v21 + ext, y.v21_confidence)


@click.group()
def cli():
    pass


@cli.command(name="transcripts")
@click.option(
    "--sleep",
    default=1.0,
    help="wait sleep seconds between requests",
    show_default=True,
)
@click.option(
    "-o",
    "--out",
    help="output CSV file",
    type=click.Path(dir_okay=False, file_okay=True),
)
@click.option(
    "--low-confidence",
    is_flag=True,
    help="include low confidence genes",
    show_default=True,
)
@click.argument(
    "idfilename",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
)
def run_cmd(
    idfilename: str,
    sleep: float,
    low_confidence: bool,
    out: str | Path | None,
) -> None:
    "find some ids"
    pd = get_v21_ids(idfilename, sleep=sleep, low_confidence=low_confidence)
    if not low_confidence:
        pd = pd.drop(columns=["confidence"])
    if not out:
        pth = Path(idfilename)
        out = pth.parent / (pth.stem + ".csv")
    click.secho(f'writing: "{out}"', fg="green")
    pd.to_csv(out, index=False)


@cli.command("conv1v2")
@click.option(
    "-o",
    "--out",
    help="output CSV file",
    type=click.Path(dir_okay=False, file_okay=True),
)
@click.option(
    "--low-confidence",
    is_flag=True,
    help="include low confidence genes",
    show_default=True,
)
@click.argument(
    "idfilename",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
)
def cvt_cmd(idfilename: str, low_confidence: bool, out: str | Path | None) -> None:
    "convert v1.1 IDs to v2.1"
    pd = convert(idfilename, low_confidence=low_confidence)
    if not low_confidence:
        pd = pd.drop(columns=["confidence"])
    if not out:
        pth = Path(idfilename)
        out = pth.parent / (pth.stem + ".csv")
    click.secho(f'writing: "{out}"', fg="green")
    pd.to_csv(out, index=False)


if __name__ == "__main__":
    cli()
