#!/usr/bin/env bash
CPU=8
cd DataSimulated
for EXPERIMENT in *.yaml; do
  python3 simulated_experiment.py -c ${EXPERIMENT} -j ${CPU}
done

cd ../DataEmpirical
snakemake --printshellcmds --rerun-incomplete -j 8
