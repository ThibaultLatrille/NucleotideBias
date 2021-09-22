# Data Empirical

## Replicate experiments

Each folder in `DataEmpirical` contains an alignment in fasta format named `alignment.fasta` and a tree topology (branch lengths are estimated) in newick format named `tree.newick`.
Snakemake loops through all folder containing these two files and then run the MG, MF and GTR (on third positions) models for each folder.
The results of the different models are pooled in the file `results.tsv` for each folder.
Simply run `snakemake` to replicate the experiments shown in the manuscript (table 1):
```
cd DataEmpirical
snakemake -j 8 -k --printshellcmds
```

## Run your own analysis
To run your own analysis, create a folder containing an alignment in fasta format named `alignment.fasta` and a tree topology (branch lengths are estimated) in newick format named `tree.newick`.
Then simply run `snakemake` in the folder `DataEmpirical`:
```
cd DataEmpirical
snakemake -j 8 -k --printshellcmds
```
Snakemake will run the MG, MF and GTR (on third positions) models on your alignement given the tree.
The results of the different models are pooled in the file `results.tsv` in your folder.

You can also remove folders for which you don't snakemake to run.
