#!/usr/bin/env bash
for EXPERIMENT in ./DataSimulated/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  rm -rf ${EXPERIMENT}/Snakefile
  ln -s $(pwd)/DataSimulated/Snakefile ${EXPERIMENT}
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake --touch -k simulation
  rm -rf */*_run.bf
  snakemake build -j 8
  snakemake --touch inference
  snakemake --printshellcmds -k --rerun-incomplete -j 8
  cd ../../..
done