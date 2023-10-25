# blasttools

tools for running blasts. i.e. turning blast queries into pandas dataframes.

Install with

```sh
python -m pip install -U 'git+https://github.com/arabidopsis/blasttools.git'
```

Once installed you can update with `blasttools update`

## Common Usages:

Build some blast databases from Ensembl Plants

```sh
blasttools plants --release=40 build triticum_aestivum zea_mays
```

Blast against my.fasta and save dataframe as pickle file.

```sh
blasttools plants blast --out=dataframe.pkl my.fasta triticum_aestivum zea_mays
```

Get your blast data!

```python
import pandas as pd
df = pd.read_pickle('dataframe.pkl')
```

## Parallelization

When blasting you can specify `--num-threads` which is passed directly to the
blast command. If you want to parallelize over species or fasta files currently
I suggest you use [GNU Parallel](https://www.gnu.org/software/parallel/) [[Tutorial](https://blog.ronin.cloud/gnu-parallel/)]
e.g. build blast databases concurrently:

```sh
parallel blasttools build ::: *.fa.gz
```

Or build *everything*!

```sh
parallel blasttools plants build ::: $(blasttools plants species)
```
