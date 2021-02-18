# GLOBAL IMPORTS
import argparse
import pandas as pd
from plot_module import *
from hyphy_format import *
import statsmodels.api as sm
from stat_simulated import stats_from_ali, open_ali_file, omega_pairwise_from_profile


def plot_pairwise_matrices(predicted, estimated, output):
    fig, axs = plt.subplots(1, 3, figsize=(16, 6))
    fig.colorbar(axs[0].imshow(predicted), ax=axs[0], orientation='horizontal', fraction=.1)
    axs[0].set_title('$\\omega$ predicted between pairs of amino-acids')
    fig.colorbar(axs[1].imshow(estimated), ax=axs[1], orientation='horizontal', fraction=.1)
    axs[1].set_title('$\\widehat{\\omega}$ estimated between pairs of amino-acids')

    filt = np.logical_and(np.isfinite(predicted) & np.isfinite(estimated), estimated > 0)
    x = predicted[filt].flatten()
    y = estimated[filt].flatten()
    axs[2].scatter(x, y)
    axs[2].set_title('$\\omega$ between pairs of amino-acids')

    model = sm.OLS(y, sm.add_constant(x))
    results = model.fit()
    b, a = results.params[0:2]
    idf = np.linspace(min(x), max(x), 100)
    axs[2].plot(idf, a * idf + b, '-',
                label=r"$y={0:.2g}x {3} {1:.2g}$ ($r^2={2:.2g})$".format(a, abs(b), results.rsquared,
                                                                         "+" if float(b) > 0 else "-"))
    axs[2].set_xlabel('Predicted')
    axs[2].set_ylabel('Estimated')
    axs[2].legend()
    fig.tight_layout()
    plt.savefig(output + ".pdf", format="pdf")
    plt.savefig(output + ".png", format="png")
    plt.clf()
    plt.close('all')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--trace', required=True, type=str, nargs='+', dest="input")
    parser.add_argument('-m', '--model', required=True, type=str, dest="model")

    args = parser.parse_args()
    nested_dict = nested_dict_init()
    predicted_array, estimated_array = [], []
    for batch in args.input:
        at_gc_pct = float(batch.split("/")[-1].split("_")[0])
        alpha = float(batch.split("/")[-2])
        exp = batch.replace(args.model + "_run.bf_hyout.txt", "exp")
        species, alignment = open_ali_file(exp + ".ali")
        ali_dico = stats_from_ali(alignment)
        nested_dict["AT/GC_obs"][at_gc_pct] = ali_dico["at_over_gc"]

        hyphy_dico = dico_from_file(batch)
        format_hyphy_dico(hyphy_dico, args.model)

        if args.model == "MF":
            profile_path = "/".join(batch.split("/")[:-1]) + "_profile.prefs"
            predicted_array.append(omega_pairwise_from_profile(profile_path, at_gc_pct))
            estimated_array.append(omega_pairwise_from_hyphy(hyphy_dico))
            plot_pairwise_matrices(predicted_array[-1], estimated_array[-1],
                                   "{0}/omega.{1}".format(args.output, at_gc_pct))

        results_dico = {k: v[0] for k, v in pd.read_csv(exp + ".tsv", sep='\t').items()}
        nested_dict["w_obs"][at_gc_pct] = results_dico["dnd0_event_tot"]
        nested_dict["w_inf"][at_gc_pct] = hyphy_dico["w"]

        if ("pnG" in hyphy_dico) and ("pnG" in hyphy_dico):
            gc_pct = hyphy_dico["pnG"] + hyphy_dico["pnC"]
            nested_dict["lambda_inf"][at_gc_pct] = (1 - gc_pct) / gc_pct

        nested_dict["AT/GC_inf"][at_gc_pct] = equilibrium_lambda(hyphy_dico)

    if args.model == "MF":
        estimated_mean = np.mean(estimated_array, axis=0)
        predicted_mean = np.mean(predicted_array, axis=0)
        plot_pairwise_matrices(predicted_mean, estimated_mean, "{0}/mean.omega".format(args.output))

    fig, ax = plt.subplots()
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
    plt.savefig("{0}/lambda.pdf".format(args.output), format="pdf")
    plt.savefig("{0}/lambda.png".format(args.output), format="png")
    plt.clf()
    plt.close('all')
