from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import Iterable
from typing import Literal
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


def get_21_conversions() -> pd.DataFrame:
    iwgsc = pd.read_csv(
        "https://urgi.versailles.inrae.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v2.1/iwgsc_refseq_all_correspondances.zip",
        sep=" ",
    )

    # iwgsc = pd.read_csv('iwgsc_refseq_all_correspondances.csv', sep=' ')

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


def read_genes(filename: str) -> list[tuple[str, list[str]]]:
    current: list[str] = []
    idata: list[tuple[str, list[str]]] = []
    desc = None
    with open(filename, encoding="utf8") as fp:
        for r in fp:
            r = r.strip()
            if not r:
                continue
            if r.startswith(">"):
                if desc:
                    idata.append((desc, current))
                current = []
                desc = r[1:]
            else:
                current.append(r)
    if desc and current:
        idata.append((desc, current))
    return idata


def findall(
    gene_names: list[tuple[str, list[str]]],
    convert: pd.DataFrame,
    *,
    sleep: float = 0.0,
    lc: bool = False,
):
    n = len(gene_names)
    for idx, (desc, ids) in enumerate(gene_names):
        idss = set(ids)
        click.secho(f"{desc}[{len(idss)}]: {idx + 1}/{n}", fg="yellow", nl=False)
        res = archive_ids(idss)
        click.secho(" done", fg="green")
        for r in res:
            idss.remove(r.id)
            for p in r.possible_replacement:
                q = p.stable_id == convert["v11"]
                if not lc:
                    q = q & (convert["v21_confidence"] == "HC")
                v2ids = convert[q][["v21", "v21_confidence"]]
                if len(v2ids) == 0:
                    yield (desc, r.id, p.stable_id, "missing", "")
                else:
                    for t in v2ids.itertuples():
                        yield (desc, r.id, p.stable_id, t.v21, t.v21_confidence)
        if idss:
            for iid in ids:
                yield (desc, iid, "missing", "missing", "")
        if sleep and idx != n - 1:
            time.sleep(sleep)


def run(idfilename: str, *, sleep: float = 0.0, lc: bool = False) -> pd.DataFrame:
    gene_names = read_genes(idfilename)
    convert = get_21_conversions()
    columns = [
        "description",
        "transcript_id",
        "stable_id_v1",
        "stabled_id_v2",
        "confidence",
    ]
    df = pd.DataFrame(
        list(findall(gene_names, convert, sleep=sleep, lc=lc)),
        columns=columns,
    )

    return df


@click.command()
@click.option(
    "--sleep",
    default=1.0,
    help="wait sleep seconds between requests",
    show_default=True,
)
@click.option(
    "--lc",
    is_flag=True,
    help="include low  genes",
    show_default=True,
)
@click.argument(
    "idfilename",
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
)
def run_cmd(idfilename: str, sleep: float, lc: bool) -> None:
    "find some ids"
    pd = run(idfilename, sleep=sleep, lc=lc)
    pth = Path(idfilename)
    pd.to_csv(pth.parent / (pth.stem + ".csv"), index=False)


if __name__ == "__main__":
    run_cmd()
