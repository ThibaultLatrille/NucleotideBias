#!/usr/bin/env bash
rm -rf ./DataEmpirical/*/*.stderr
rm -rf ./DataEmpirical/*/*.stdout
rm -rf ./DataSimulated/Experiments/*/slurm
rm -rf ./DataSimulated/Experiments/*/benchmarks
rm -rf ./DataSimulated/Experiments/*/std
rm -rf ./DataSimulated/Experiments/*/*/slurm
rm -rf ./DataSimulated/Experiments/*/*/benchmarks
rm -rf ./DataSimulated/Experiments/*/*/std
rm -rf ./DataSimulated/Experiments/*/*/*.fossils.tsv
rm -rf ./DataSimulated/Experiments/*/*/*.substitutions.tsv
rm -rf ./DataSimulated/Experiments/*/*/*.DFEdistrib.tsv
rm -rf ./DataSimulated/Experiments/*/*/*.sequences.tsv
rm -rf ./DataSimulated/Experiments/*/*/*.traits.tsv