# blasttools

tools for running blasts. i.e. turning blast queries into pandas dataframes.

Install with

```sh
python -m pip install -U 'git+https://github.com/arabidopsis/blasttools.git'
```

Once installed you can update with `blasttools update`

## Common Usages:

Build some blast databases from Ensembl Plants.

```sh
blasttools plants --release=40 build triticum_aestivum zea_mays
```

Find out what species are available:

```sh
blasttools plants --release=40 species
```

Blast against my.fasta and save dataframe as pickle file.

```sh
blasttools plants blast --out=dataframe.pkl my.fasta triticum_aestivum zea_mays
# OR just all downloaded databases
blasttools plants blast --all --out=dataframe.pkl my.fasta
```

Get your blast data!

```python
import pandas as pd
df = pd.read_pickle('dataframe.pkl')
```

## Parallelization

When blasting, you can specify `--num-threads` which is passed directly to the
blast command. If you want to parallelize over species, databases or fasta files, currently
I suggest you use [GNU Parallel](https://www.gnu.org/software/parallel/) [[Tutorial](https://blog.ronin.cloud/gnu-parallel/)]
e.g. build blast databases concurrently:

```sh
parallel blasttools build ::: *.fa.gz
```

Or blast _everything_!

```sh
species=$(blasttools plants species)
parallel blasttools plants build ::: $species
# must have different output files here...
parallel blasttools plants blast --out=my{}.pkl my.fasta ::: $species
# or in batches of 4 species at a time
parallel -N4 blasttools plants blast --out=my{#}.pkl my.fasta ::: $species
```

Remember if you parallelize your blasts _and_ use `--num-threads > 1`
then you are probably going to be fighting for cpu time
amongst yourselves!
