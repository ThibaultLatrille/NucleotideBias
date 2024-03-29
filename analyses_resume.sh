#!/usr/bin/env bash
for EXPERIMENT in DataSimulated/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --unlock
  snakemake --printshellcmds --rerun-incomplete -j 8 -k
  cd ../../..
done