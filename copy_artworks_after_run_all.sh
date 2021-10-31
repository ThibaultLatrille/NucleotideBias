#!/usr/bin/env bash
# Simulations
rm -rf ./manuscript/artworks/simulations/
mkdir ./manuscript/artworks/simulations/

\cp -f ./DataSimulated/Experiments/PrimatesSimulationsNucBiasExons10Mu1.0/plot_simulation/at_over_gc.pdf ./manuscript/artworks/simulations/at_over_gc.pdf
\cp -f ./DataSimulated/Experiments/PrimatesSimulationsNucBiasExons10Mu1.0/plot_simulation/at_over_gc_1.pdf ./manuscript/artworks/simulations/at_over_gc_1.pdf
\cp -f ./DataSimulated/Experiments/PrimatesSimulationsNucBiasExons10Mu1.0/plot_simulation/at_over_gc_2.pdf ./manuscript/artworks/simulations/at_over_gc_2.pdf
\cp -f ./DataSimulated/Experiments/PrimatesSimulationsNucBiasExons10Mu1.0/plot_simulation/at_over_gc_3.pdf ./manuscript/artworks/simulations/at_over_gc_3.pdf

\cp -f ./DataSimulated/Experiments/PrimatesSimulationsNucBiasExons10Mu1.0/plot_simulation/omega.pdf ./manuscript/artworks/simulations/omega.pdf
\cp -f ./DataSimulated/Experiments/PrimatesSimulationsNucBiasExons10Mu1.0/plot_simulation/omega_SW.pdf ./manuscript/artworks/simulations/omega_SW.pdf
\cp -f ./DataSimulated/Experiments/PrimatesSimulationsNucBiasExons10Mu1.0/plot_simulation/omega_WS.pdf ./manuscript/artworks/simulations/omega_WS.pdf
\cp -f ./DataSimulated/Experiments/PrimatesSimulationsNucBiasExons10Mu1.0/plot_simulation/omega_WS_over_SW.pdf ./manuscript/artworks/simulations/omega_WS_over_SW.pdf
\cp -f ./DataSimulated/Experiments/PrimatesSimulationsNucBiasExons10Mu1.0/plot_simulation/omega_SW_over_WS.pdf ./manuscript/artworks/simulations/omega_SW_over_WS.pdf

# Inference
rm -rf ./manuscript/artworks/inference_simulations/
mkdir ./manuscript/artworks/inference_simulations/

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MF_plot/obs_atgc_1.pdf ./manuscript/artworks/inference_simulations/obs_atgc_1_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MF_plot/obs_atgc_2.pdf ./manuscript/artworks/inference_simulations/obs_atgc_2_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MF_plot/obs_atgc_3.pdf ./manuscript/artworks/inference_simulations/obs_atgc_3_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MG_plot/obs_atgc_1.pdf ./manuscript/artworks/inference_simulations/obs_atgc_1_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MG_plot/obs_atgc_2.pdf ./manuscript/artworks/inference_simulations/obs_atgc_2_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MG_plot/obs_atgc_3.pdf ./manuscript/artworks/inference_simulations/obs_atgc_3_MG.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_simulations/lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_simulations/lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_simulations/lambda_GTR.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_simulations/omega_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_simulations/omega_MF.pdf

# Inference replicates
\cp -f ./DataSimulated/Experiments/PrimatesRepExons10Mu1.0/2.000000/MF_plot/mean.omega.pdf ./manuscript/artworks/inference_simulations/omega_pair_MF.pdf

# Empirical inference
\cp -f ./DataEmpirical/merged_results.T.tex ./manuscript/artworks/empirical_results.T.tex
\cp -f ./DataEmpirical/merged_results.tex ./manuscript/artworks/empirical_results.tex
\cp -f ./DataEmpirical/merged_results.T.tsv ./manuscript/artworks/empirical_results.T.tsv
\cp -f ./DataEmpirical/merged_results.tsv ./manuscript/artworks/empirical_results.tsv

# Supplementary Materials
rm -rf ./manuscript/artworks/inference_supp_mat/
mkdir ./manuscript/artworks/inference_supp_mat/

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons1Mu1.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons1Mu1.0_lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons1Mu1.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons1Mu1.0_lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons1Mu1.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons1Mu1.0_lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons1Mu1.0/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons1Mu1.0_omega_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons1Mu1.0/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons1Mu1.0_omega_MF.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons2Mu1.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons2Mu1.0_lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons2Mu1.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons2Mu1.0_lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons2Mu1.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons2Mu1.0_lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons2Mu1.0/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons2Mu1.0_omega_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons2Mu1.0/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons2Mu1.0_omega_MF.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons5Mu1.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons5Mu1.0_lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons5Mu1.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons5Mu1.0_lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons5Mu1.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons5Mu1.0_lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons5Mu1.0/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons5Mu1.0_omega_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons5Mu1.0/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons5Mu1.0_omega_MF.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons20Mu1.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons20Mu1.0_lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons20Mu1.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons20Mu1.0_lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons20Mu1.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons20Mu1.0_lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons20Mu1.0/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons20Mu1.0_omega_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons20Mu1.0/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons20Mu1.0_omega_MF.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu0.5/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu0.5_lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu0.5/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu0.5_lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu0.5/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu0.5_lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu0.5/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu0.5_omega_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu0.5/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu0.5_omega_MF.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu2.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu2.0_lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu2.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu2.0_lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu2.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu2.0_lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu2.0/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu2.0_omega_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu2.0/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu2.0_omega_MF.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu4.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu4.0_lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu4.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu4.0_lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu4.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu4.0_lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu4.0/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu4.0_omega_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu4.0/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu4.0_omega_MF.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu8.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu8.0_lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu8.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu8.0_lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu8.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu8.0_lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu8.0/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu8.0_omega_MG.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu8.0/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu8.0_omega_MF.pdf

\cp -f ./DataSimulated/Experiments/MammalsNucBiasExons10Mu1.0/2.000000/MF_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/MammalsExons10Mu1.0_lambda_MF.pdf
\cp -f ./DataSimulated/Experiments/MammalsNucBiasExons10Mu1.0/2.000000/MG_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/MammalsExons10Mu1.0_lambda_MG.pdf
\cp -f ./DataSimulated/Experiments/MammalsNucBiasExons10Mu1.0/2.000000/GTR_plot/lambda.pdf ./manuscript/artworks/inference_supp_mat/MammalsExons10Mu1.0_lambda_GTR.pdf
\cp -f ./DataSimulated/Experiments/MammalsNucBiasExons10Mu1.0/2.000000/MG_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/MammalsExons10Mu1.0_omega_MG.pdf
\cp -f ./DataSimulated/Experiments/MammalsNucBiasExons10Mu1.0/2.000000/MF_plot/omega.pdf ./manuscript/artworks/inference_supp_mat/MammalsExons10Mu1.0_omega_MF.pdf

\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons1Mu1.0/2.000000/MF_plot/mean.omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons1Mu1.0_omega_pair_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons2Mu1.0/2.000000/MF_plot/mean.omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons2Mu1.0_omega_pair_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons5Mu1.0/2.000000/MF_plot/mean.omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons5Mu1.0_omega_pair_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons10Mu1.0/2.000000/MF_plot/mean.omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons10Mu1.0_omega_pair_MF.pdf
\cp -f ./DataSimulated/Experiments/PrimatesNucBiasExons20Mu1.0/2.000000/MF_plot/mean.omega.pdf ./manuscript/artworks/inference_supp_mat/PrimatesExons20Mu1.0_omega_pair_MF.pdf


