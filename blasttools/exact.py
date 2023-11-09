from __future__ import annotations

from collections.abc import Sequence
from collections import defaultdict
from typing import Any

import click
import pandas as pd
from .blastapi import fasta_to_df, read_fasta, save_df
from .cli import blast


def exact(query: str, subjects: Sequence[str]) -> pd.DataFrame:
    df = fasta_to_df(query, with_description=True)
    dd = defaultdict(list)
    for subject in subjects:
        for rec in read_fasta(subject):
            for r in df.itertuples():
                if str(r.seq).upper() in str(rec.seq).upper():
                    dd[r.id].append(rec.id)

    ret: dict[str, list[Any]] = dict(id=[],  description=[], exact=[], subject=[])
    for r in df.itertuples():
        found = r.id in dd
        ret["id"].append(r.id)
        ret["description"].append(r.description)
        ret["exact"].append(found)
        if found:
            ret["subject"].append(", ".join(dd[r.id]))
        else:
            ret["subject"].append("")

    return pd.DataFrame(ret)


@blast.command("exact")
@click.option(
    "--out",
    help="output filename (default is to write <query>.csv)",
    type=click.Path(dir_okay=False),
)
@click.argument("query", type=click.Path(exists=True, dir_okay=False))
@click.argument("fastas", nargs=-1)
def exact_cmd(query: str, fastas: Sequence[str], out: str | None) -> None:
    "Try and match sequences exactly"
    df = exact(query, fastas)
    if out is None:
        out = query + ".csv"
    save_df(df, out, index=False)
