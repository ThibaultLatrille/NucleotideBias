#!/usr/bin/env bash
for EXPERIMENT in DataSimulated/Experiments/*; do
  NAME=$(basename "${EXPERIMENT}")
  echo "${NAME}"
  cd ${EXPERIMENT}
  snakemake --printshellcmds -n
  cd ../../..
done