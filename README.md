# blasttools

tools for running blasts. i.e. turning blast queries into pandas dataframes.

Install with

```sh
python -m pip install -U 'git+https://github.com/arabidopsis/blasttools.git'
```

Once installed you can update with `blasttools update`

## Common usages:

Build some blast databases from Ensembl Plants

```sh
blasttools plants --release=40 build triticum_aestivum zea_mays
```

Blast again my.fasta and save dataframe as pickle file.

```sh
blasttools plants blast --out=dataframe.pkl my.fasta triticum_aestivum zea_mays
```
