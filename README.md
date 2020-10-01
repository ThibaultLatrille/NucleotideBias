# Nucleotide Bias

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
Install the compiling toolchains:
```
sudo apt install -qq -y make cmake clang
```
Clone and compile the C++ code for *SimuEvol*
```
git clone https://github.com/ThibaultLatrille/SimuEvol && cd SimuEvol && make release && cd ..
```
Install python3 packages
```
sudo apt install -qq -y python3-dev python3-pip
pip3 install jupyterlab snakemake numpy matplotlib pandas ete3 --user
```
Install Hyphy from http://hyphy.org/installation/ using `miniconda`
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
conda config --add channels bioconda
conda install -c bioconda hyphy
```
## 2. Run experiments

### Installation on debian
In both folder `DataEmpirical` and `DataSimulated`, run `snakemake`:
```
cd DataSimulated
snakemake
```

## 3. Add features or debug in the python scripts
You made modifications to one of the python script, a notebook, this README.md, or you added new features.
You wish this work benefits to all (futur) users of this repository?
Please, feel free to open a [pull-request](https://github.com/ThibaultLatrille/NucleotideBias/pulls)

## 4. Add features or debug in *SimuEvol*
You made modifications to the C++ code of the simulation framework *SimuEvol*.
You wish this changes benefit to all users of these software?

Please, feel free to open pull-requests in the respective GitHub repository:
* https://github.com/ThibaultLatrille/SimuEvol 

## Licence

The MIT License (MIT)

Copyright (c) 2019 Thibault Latrille

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
