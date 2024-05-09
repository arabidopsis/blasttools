from __future__ import annotations

import sys

import click
from Bio import SeqIO

#
# only called by blast command


@click.command()
@click.argument("fasta", type=click.Path(dir_okay=False))
def translate(fasta: str) -> None:
    from .utils import translate as trans

    SeqIO.write(trans(fasta), sys.stdout, "fasta")


if __name__ == "__main__":
    translate()  # pylint disable=no-value-for-parameter
