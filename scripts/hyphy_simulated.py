# GLOBAL IMPORTS
import argparse
import pandas as pd
from plot_module import *
from hyphy_format import *
from stat_simulated import stats_from_ali, open_ali_file

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--trace', required=True, type=str, nargs='+', dest="input")
    parser.add_argument('-m', '--model', required=True, type=str, dest="model")

    args = parser.parse_args()
    nested_dict = nested_dict_init()
    for batch in args.input:
        at_gc_pct = float(batch.split("/")[-1].split("_")[0])

        exp = batch.replace(args.model + "_run.bf", "exp")
        ali_dico = stats_from_ali(open_ali_file(exp + ".ali"))
        nested_dict["AT/GC_obs"][at_gc_pct] = ali_dico["at_over_gc"]

        hyphy_result = "{0}_hyout.txt".format(batch)
        hyphy_dico = dico_from_file(hyphy_result)
        format_hyphy_dico(hyphy_dico)

        results_dico = {k: v[0] for k, v in pd.read_csv(exp + ".tsv", sep='\t').items()}
        nested_dict["w_obs"][at_gc_pct] = results_dico["dnd0_event_tot"]
        nested_dict["w_inf"][at_gc_pct] = hyphy_dico["w"]

        if ("pnG" in hyphy_dico) and ("pnG" in hyphy_dico):
            gc_pct = hyphy_dico["pnG"] + hyphy_dico["pnC"]
            nested_dict["lambda_inf"][at_gc_pct] = (1 - gc_pct) / gc_pct

        nested_dict["AT/GC_inf"][at_gc_pct] = equilibrium_lambda(hyphy_dico)

    my_dpi = 96
    fig, ax = plt.subplots()
    index = 0
    list_plot = list()
    list_plot.append(
        {"experiment": "AT/GC_obs", "color": BLUE, "linestyle": '-', "linewidth": 2,
         "label": "$AT/GC$ observed"})
    list_plot.append({"experiment": "AT/GC_inf", "color": YELLOW, "linestyle": '--', "linewidth": 4,
                      "label": "$\\widehat{AT/GC}$ inferred"})
    list_plot.append({"experiment": "lambda_inf", "color": GREEN, "linestyle": '--', "linewidth": 4,
                      "label": "$\\widehat{\\lambda}$ inferred"})

    x_list = sorted(nested_dict["AT/GC_obs"].keys())
    ax.plot(x_list, x_list, color="black", linestyle='-', linewidth=2, label="y=x")

    omega_obs = np.array([nested_dict["w_obs"][k] for k in x_list])
    omega_inf = np.array([nested_dict["w_inf"][k] for k in x_list])
    print("|w_obs-w_inf|/w_obs = {0:.2f}%".format(100 * np.mean(np.abs((omega_inf - omega_obs)) / omega_obs)))

    for param in list_plot:
        x_list = sorted(nested_dict[param["experiment"]].keys())
        y_list = [nested_dict[param["experiment"]][k] for k in x_list]
        ax.plot(x_list, y_list, linestyle=param["linestyle"], label=param["label"],
                color=param["color"], linewidth=param["linewidth"])

    lambda_obs = np.array(x_list)
    lambda_inf = np.array([nested_dict["lambda_inf"][k] for k in x_list])
    print("|lambda_obs-lambda_inf|/lambda_obs = {0:.2f}%".format(
        100 * np.mean(np.abs((lambda_inf - lambda_obs)) / lambda_obs)))

    ax.set_xscale('log')
    ax.set_xlabel('$\\lambda$ used for the simulation')
    ax.set_yscale('log')
    ax.set_ylabel('$\\lambda$ inferred')
    ax.legend()
    ax.set_title(args.model)
    plt.tight_layout()
    plt.savefig("{0}.pdf".format(args.output), format="pdf")
    plt.savefig("{0}.png".format(args.output), format="png")
    plt.clf()
    plt.close('all')
