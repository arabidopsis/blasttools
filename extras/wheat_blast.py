from __future__ import annotations

import re
from pathlib import Path

import click
import pandas as pd

HUI_XLSX = "Annotation_Wheat_proteins_Hui.xlsx"


def hui_xlsx(annotations: str = HUI_XLSX) -> pd.DataFrame:
    hui = pd.read_excel(annotations)
    hui["uprotein_id"] = hui["protein_id"].str.upper()
    prots = hui[hui["uprotein_id"].str.match("^TRAES.*$")][["uprotein_id"]]
    assert prots["uprotein_id"].str.match("^TRAESCS.*$").all()
    assert hui["uprotein_id"].is_unique
    return hui


def run_wheat(
    filename: str,
    pident: float = 95.0,
    annotations: str = HUI_XLSX,
) -> pd.DataFrame:
    df = pd.read_pickle(filename)
    df["uprotein_id"] = df["saccver"].str.upper()
    df["uprotein_id"] = df["uprotein_id"].str.extract(
        r"^(.*?)(\.cds1)?$",
        flags=re.IGNORECASE,
    )[0]
    hui = hui_xlsx(annotations)
    mdf = df.merge(hui, left_on="uprotein_id", right_on="uprotein_id", how="left")

    b99 = mdf[mdf["pident"] > pident]

    u = b99["saccver"].drop_duplicates()
    g = mdf[~mdf["saccver"].isin(u)]
    v = (
        g[["saccver", "pident"]]
        .groupby("saccver")["pident"]
        .nlargest(2)
        .reset_index(level=0)
    )
    g = g.loc[v.index]
    odf = pd.concat([b99, g], axis=0)
    odf = odf.drop(columns=["uprotein_id", "subject_seq", "seq", "species"])
    cols = [
        "qaccver",
        "saccver",
        "protein_id",
        "gene_id",
        "protein_name",
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
        "description",
        "protein_id",
        "gene_id",
        "match_to_refSeq",
        "ort_source",
        "Athaliana ortholog",
        "SUBA_live",
        "SUBA_short1",
        "SUBA_short2",
        "SUBA_short3",
        "percent identity",
        "ort_homology_type",
        "Blast e-val",
        "functional_group",
    ]
    odf = odf[cols]
    return odf


@click.command()
@click.option(
    "--annotations",
    default=HUI_XLSX,
    help="Hui's annotation file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    show_default=True,
)
@click.option("--pident", default=95.0, help="percent identity")
@click.argument(
    "picklefile",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
)
def to_csv(picklefile: str, pident: float = 95.0, annotations: str = HUI_XLSX):
    odf = run_wheat(picklefile, pident, annotations)
    out = Path(picklefile)
    out = out.parent / (out.stem + ".csv.gz")
    click.secho(f'writing "{out}"')
    odf.to_csv(out, index=False)


if __name__ == "__main__":
    to_csv()
