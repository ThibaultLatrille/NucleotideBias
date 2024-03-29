# GLOBAL IMPORTS
import argparse
import pandas as pd
from hyphy_format import *
from stat_simulated import stats_from_ali, open_fasta_file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--input', required=True, type=str, nargs='+', dest="input")
    args = parser.parse_args()

    species, ali = open_fasta_file("/".join(args.output.split("/")[:-1]) + "/alignment.fasta")
    out_dico = stats_from_ali(ali)
    out_dico["name"] = args.output.split("/")[-2]

    for hyphy_result in args.input:
        model = hyphy_result.split("/")[-1].split("_")[0]
        hyphy_dico = dico_from_file(hyphy_result)
        format_hyphy_dico(hyphy_dico, model)

        out_dico["lambda_obs_" + model], out_dico["lambda_obs_1" + model], out_dico["lambda_obs_2" + model], out_dico[
            "lambda_obs_3" + model] = equilibrium_lambda(hyphy_dico)
        if ("pnG" in hyphy_dico) and ("pnG" in hyphy_dico):
            gc_pct = (hyphy_dico["pnG"] + hyphy_dico["pnC"])
            out_dico["lambda_" + model] = (1.0 - gc_pct) / gc_pct

        out_dico["AIC_" + model] = hyphy_dico["AIC"]
        out_dico["LnL_" + model] = hyphy_dico["Log Likelihood"]
        out_dico["df_" + model] = hyphy_dico["df"]

        if model == "GTR": continue
        out_dico["w_" + model] = hyphy_dico["w"]

        for subset in subset_list:
            omega = "w_" + subset
            if omega in hyphy_dico:
                out_dico[omega + "_" + model] = hyphy_dico[omega]

    pd.DataFrame({k: [v] for k, v in out_dico.items()}).to_csv(args.output, index=False, sep="\t")
