#!/usr/bin/env bash
# Simulations
rm -rf ./manuscript/artworks/simulations/
mkdir ./manuscript/artworks/simulations/
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/at_over_gc.pdf ./manuscript/artworks/simulations/at_over_gc.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/at_over_gc_1.pdf ./manuscript/artworks/simulations/at_over_gc_1.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/at_over_gc_2.pdf ./manuscript/artworks/simulations/at_over_gc_2.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/at_over_gc_3.pdf ./manuscript/artworks/simulations/at_over_gc_3.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/diversity_aa.pdf ./manuscript/artworks/simulations/diversity_aa.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/omega.pdf ./manuscript/artworks/simulations/omega.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/omega_SW.pdf ./manuscript/artworks/simulations/omega_SW.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/omega_WS.pdf ./manuscript/artworks/simulations/omega_WS.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/omega_WS_over_SW.pdf ./manuscript/artworks/simulations/omega_WS_over_SW.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/omega_SW_over_WS.pdf ./manuscript/artworks/simulations/omega_SW_over_WS.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons100Mu1.0/plot_simulation/site_diversity_aa.pdf ./manuscript/artworks/simulations/site_diversity_aa.pdf

# Inference
rm -rf ./manuscript/artworks/inference_simulations/
mkdir ./manuscript/artworks/inference_simulations/
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_simulations/lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_simulations/lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_simulations/lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/PrimatesRepExons10Mu1.0/2.000000/MF_plot/mean.omega.pdf ./manuscript/artworks/inference_simulations/omega_MF.pdf

# Empirical inference
\cp -f ./DataEmpirical/merged_results.T.tex ./manuscript/artworks/empirical_results.T.tex
\cp -f ./DataEmpirical/merged_results.tex ./manuscript/artworks/empirical_results.tex
\cp -f ./DataEmpirical/merged_results.T.tsv ./manuscript/artworks/empirical_results.T.tsv
\cp -f ./DataEmpirical/merged_results.tsv ./manuscript/artworks/empirical_results.tsv


