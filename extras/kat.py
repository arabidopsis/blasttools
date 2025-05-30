from __future__ import annotations

from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import Seq
from Bio.SeqRecord import SeqRecord


@click.group()
def kat():
    pass


@kat.command()
@click.argument("csvfile", type=click.Path(exists=True, dir_okay=False))
def convert(csvfile):
    """convert CSV file to fasta"""
    df = pd.read_csv(csvfile)
    assert df.columns == ["Gene.stable.ID", "Protein_stable_ID", "AA_combined_f"]
    df.columns = ["gene", "id", "seq"]
    assert not df["id"].str.match(r"\s+").any(), "ids contain spaces"
    assert df["id"].is_unique, "non unique ids"
    out = f"{csvfile}.fasta"
    click.secho(f"writing: {out}", fg="green")
    with open(out, "w", encoding="ascii") as fp:
        for t in df.itertuples():
            rec = SeqRecord(Seq(t.seq), t.id, description=t.gene)
            SeqIO.write(rec, fp, "fasta")


@kat.command()
@click.option("--toxlsx", is_flag=True, help="output Excel")
@click.argument("pklfile", type=click.Path(exists=True, dir_okay=False))
def tocsv(pklfile: str, toxlsx: bool):
    """Convert blasttools output pickle file to CSV or Excel"""
    path = Path(pklfile)
    ext = ".xlsx" if toxlsx else ".csv"
    out = path.parent / (path.stem + ext)
    df = pd.read_pickle(path)
    click.secho(f"writing: {out}", fg="green")
    df.to_csv(out)


kat()
