from __future__ import annotations

import click


def blast_options(f):

    f = click.option(
        "--needs-translation",
        is_flag=True,
        help="query fasta contains rna/dna/cdna sequences too",
    )(f)
    f = click.option(
        "--without-query-seq",
        is_flag=True,
        help="don't output query sequence",
    )(f)
    f = click.option("--xml", is_flag=True, help="run with xml output (for matches)")(f)
    f = click.option(
        "-t",
        "--num-threads",
        help="number of threads to use for blast",
        default=1,
    )(f)
    f = click.option(
        "--expr",
        help="expression to minimize when looking for --best. e.g. ",
        default="evalue",
        show_default=True,
    )(f)
    f = click.option(
        "-d",
        "--with-description",
        is_flag=True,
        help="include query description in output",
    )(f)
    f = click.option(
        "--with-seq",
        "with_subject_seq",
        is_flag=True,
        help="add subject sequence data to output",
    )(f)
    f = click.option(
        "--best",
        default=0,
        help="best (lowest) evalues [=0 take all]  (see also --expr)",
    )(f)
    return f
