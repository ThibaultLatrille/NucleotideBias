# Data Simulated

The simulator [SimuEvol](https://github.com/ThibaultLatrille/SimuEvol) is required to generate the alignments on which the MG and MF models are tested against.

## Installing *SimuEvol*

Install the compiling toolchains:
```
sudo apt install -qq -y make cmake clang
```
In the root folder of NucleotideBias, clone and compile the C++ code for *SimuEvol*
```
cd NucleotideBias
git clone https://github.com/ThibaultLatrille/SimuEvol && cd SimuEvol && make release && cd ..
```

## Replicate experiments

To reproduce figures of the manuscript, one can use the config files of the simulations (.yaml) on the folder `DataSimulated`.
```
python3 simulated_experiment.py --config PrimatesNucBiasExons10Mu1.0.yaml --nbr_cpu 4
```
The script _simulated_experiment.py_ also contains options to run the experiments on a SLURM cluster with the option `--sbatch true`, and in screen mode with the option `--screen true`.

To reproduce all the experiments of the manuscript, this loop run all of them.
```
for CONFIG in *.yaml; do
  python3 simulated_experiment.py --config "${CONFIG}" --nbr_cpu 4
done
```

Running the script `simulated_experiment.py` on the file `MySimulation.yaml` will create a folder in `DataSimulated/Experiments/MySimulation` containing the files necessary to initialize Snakemake, then run the command `snakemake`.
Hence, if the command fails or needs to be restarted, one can simply `cd` into the directory and run `snakemake`:
```
cd DataSimulated/Experiments/MySimulation
snakemake -j 8 -k --printshellcmds
```

## Add features or debug in *SimuEvol*
You made modifications to the C++ code of the simulation framework *SimuEvol*.
You wish this changes benefit to all users of these software?

Please, feel free to open pull-requests in the respective GitHub repository:
* https://github.com/ThibaultLatrille/SimuEvol
