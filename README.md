# blasttools

Commands for turning blast queries into pandas dataframes.

Blast against any built blast databases

```sh
blasttools blast --out=my.pkl query.fasta my_blastdbs_dir/*.pot
```

## Install

Install with

```sh
python -m pip install -U blasttools
# *OR*
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

Blast against `my.fasta` and save dataframe as a pickle file (the default is to
save as a csv file named `my.fasta.csv`).

```sh
blasttools plants blast --out=dataframe.pkl my.fasta triticum_aestivum zea_mays
```

Get your blast data!

```python
import pandas as pd
df = pd.read_pickle('dataframe.pkl')
```

## Parallelization

When blasting, you can specify `--num-threads` which is passed directly to the
underlying blast command. If you want to parallelize over species, databases or fasta files,
I suggest you use [GNU Parallel](https://www.gnu.org/software/parallel/) [[Tutorial](https://blog.ronin.cloud/gnu-parallel/)].

`parallel` has a much better set of options for controlling how the parallelization works
and is also quite simple for simple things.

e.g. build blast databases from a set of fasta files concurrently:

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
parallel -N4 blasttools plants blast --out='my{#}.pkl' my.fasta ::: $species
```

Then gather them all together...

```sh
blasttools concat --out=alldone.xlsx my*.pkl && rm my*.pkl
```

or programmatically:

```python
from glob import glob
import pandas as pd
df = pd.concat([pd.read_pickle(f) for f in glob('my*.pkl')], ignore_index=True)
```

Remember: if you parallelize your blasts _and_ use `--num-threads > 1`
then you are probably going to be fighting for cpu time
amongst yourselves!

## Best matches

Usually if you want the top/best `--best=3` will select the _lowest_ evalue's for
each query sequence. However if you want say the best to, say, be the longest query match
then you can add `--expr='qstart - qend'`. (Remember we are looking for the lowest values).

## XML

Blast offers an xml (`--xml`) output format that adds `query`, `match`, `sbjct` strings. The other
fields are equivalent to adding `--columns='+score gaps nident positive qlen slen'`.

It also offers a way to display the blast match as a pairwise alignment.

```python
from blasttools.blastxml import hsp_match
df = pd.read_csv('results.csv')
df['alignment'] = df.apply(hsp_match, axis=1)
print(df.iloc[0].alignment)
```
