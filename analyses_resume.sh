#!/usr/bin/env bash
for EXPERIMENT in ./DataSimulated/Experiments/Simii*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  rm -rf ${EXPERIMENT}/Snakefile
  ln -s $(pwd)/DataSimulated/Snakefile ${EXPERIMENT}
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake --touch inference
  snakemake --printshellcmds --rerun-incomplete -j 4
  cd ../../..
done
