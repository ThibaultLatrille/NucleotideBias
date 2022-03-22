**An improved codon modeling approach for accurate estimation of the mutation bias,**\
Thibault Latrille, Nicolas Lartillot,\
*Molecular Biology and Evolution*, Volume 39, Issue 2, February 2022, msac005,\
https://doi.org/10.1093/molbev/msac005

--- 
This repository is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript.
The experiments can either run on a local computer or in a cluster configuration (slurm).

The experiments are meant to run on Linux/Unix/MacOS operating systems.

If problems and/or questions are encountered, feel free to [open issues](https://github.com/ThibaultLatrille/NucleotideBias/issues).

## 0. Local copy
Clone the repository and `cd` to the dir.
```
git clone https://github.com/ThibaultLatrille/NucleotideBias
cd NucleotideBias
```

## 1. Installation

### Installation on debian

Install python3 packages
```
sudo apt install -qq -y python3-dev python3-pip
pip3 install jinja2 snakemake numpy matplotlib pandas ete3 statsmodels --user
```
Install Hyphy from http://hyphy.org/installation/ using `miniconda`
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
conda config --add channels bioconda
conda install -c bioconda hyphy
```
## 2. Replicate experiments

To replicate figure 1, 2 and 3 of the manuscript on simulated dataset, the sub-folder [*DataSimulated*](https://github.com/ThibaultLatrille/NucleotideBias/tree/master/DataSimulated) contains a [README.md](https://github.com/ThibaultLatrille/NucleotideBias/tree/master/DataSimulated/README.md) with the necessary instructions.

To replicate table 1 of the manuscript on empirical sequence data, the sub-folder [*DataEmpirical*](https://github.com/ThibaultLatrille/NucleotideBias/tree/master/DataEmpirical) contains a [README.md](https://github.com/ThibaultLatrille/NucleotideBias/tree/master/DataEmpirical/README.md) with the necessary instructions.

## 3. Run your own experiments

Instructions can be found the sub-folder [*DataEmpirical*](https://github.com/ThibaultLatrille/NucleotideBias/tree/master/DataEmpirical).

## 4. Add features or debug in the python scripts
You made modifications to one of the python script, a notebook, this README.md, or you added new features.
You wish this work benefits to all (futur) users of this repository?
Please, feel free to open a [pull-request](https://github.com/ThibaultLatrille/NucleotideBias/pulls)

## Licence

The MIT License (MIT)

Copyright (c) 2019 Thibault Latrille

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
