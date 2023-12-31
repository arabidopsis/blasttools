# valid columns for -outfmt 6
# taken from `blastp -help` output
from __future__ import annotations

# in normal blasting
# saccver == sseqid == sacc
# qaccver == qseqid == qacc
VALID = {
    "qseqid": "Query Seq-id",
    "qgi": "Query GI",
    "qacc": "Query accesion",
    "qaccver": "Query accesion.version",
    "qlen": "Query sequence length",
    "sseqid": "Subject Seq-id",
    "sallseqid": "All subject Seq-id(s), separated by a ';'",
    "sgi": "Subject GI",
    "sallgi": "All subject GIs",
    "sacc": "Subject accession",
    "saccver": "Subject accession.version",
    "sallacc": "All subject accessions",
    "slen": "Subject sequence length",
    "qstart": "Start of alignment in query",
    "qend": "End of alignment in query",
    "sstart": "Start of alignment in subject",
    "send": "End of alignment in subject",
    "qseq": "Aligned part of query sequence",
    "sseq": "Aligned part of subject sequence",
    "evalue": "Expect value",
    "bitscore": "Bit score",
    "score": "Raw score",
    "length": "Alignment length",
    "pident": "Percentage of identical matches",
    "nident": "Number of identical matches",
    "mismatch": "Number of mismatches",
    "positive": "Number of positive-scoring matches",
    "gapopen": "Number of gap openings",
    "gaps": "Total number of gaps",
    "ppos": "Percentage of positive-scoring matches",
    "frames": "Query and subject frames separated by a '/'",
    "qframe": "Query frame",
    "sframe": "Subject frame",
    "btop": "Blast traceback operations (BTOP)",
    "staxid": "Subject Taxonomy ID",
    "ssciname": "Subject Scientific Name",
    "scomname": "Subject Common Name",
    "sblastname": "Subject Blast Name",
    "sskingdom": "Subject Super Kingdom",
    "staxids": "unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)",
    "sscinames": "unique Subject Scientific Name(s), separated by a ';'",
    "scomnames": "unique Subject Common Name(s), separated by a ';'",
    "sblastnames": "unique Subject Blast Name(s), separated by a ';' (in alphabetical order)",
    "sskingdoms": "unique Subject Super Kingdom(s), separated by a ';' (in alphabetical order) ",
    "stitle": "Subject Title",
    "salltitles": "All Subject Title(s), separated by a '<>'",
    "sstrand": "Subject Strand",
    "qcovs": "Query Coverage Per Subject",
    "qcovhsp": "Query Coverage Per HSP",
    "qcovus": "Query Coverage Per Unique Subject (blastn only)",
}
