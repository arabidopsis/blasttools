from __future__ import annotations

import click


@click.command()
@click.argument("fasta", type=click.Path(dir_okay=False))
def translate(fasta: str) -> None:
    import sys
    from .utils import is_rna
    from Bio import SeqIO

    from .blastapi import read_fasta

    for rec in read_fasta(fasta):
        if is_rna(str(rec.seq)):
            rec.seq = rec.seq.transcribe().translate()
            # print('translated', rec.id,rec.seq, file=sys.stderr)
        SeqIO.write([rec], sys.stdout, "fasta")


if __name__ == "__main__":
    translate()  # pylint disable=no-value-for-parameter
