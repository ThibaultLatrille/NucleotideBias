# GLOBAL IMPORTS
import argparse
import pandas as pd
from plot_module import *
from hyphy_format import *
import statsmodels.api as sm
from stat_simulated import stats_from_ali, open_fasta_file, omega_pairwise_from_profile
from scipy.linalg import null_space


def plot_pairwise_matrices(predicted, estimated, output):
    fig, axs = plt.subplots(1, 3, figsize=(16, 6))
    fig.colorbar(axs[0].imshow(predicted), ax=axs[0], orientation='horizontal', fraction=.1)
    axs[0].set_title('$\\phi$ predicted between pairs of amino-acids')
    fig.colorbar(axs[1].imshow(estimated), ax=axs[1], orientation='horizontal', fraction=.1)
    axs[1].set_title('$\\widehat{\\omega}$ estimated between pairs of amino-acids')

    for index in [0, 1]:
        # We want to show all ticks...
        axs[index].set_xticks(np.arange(len(amino_acids)))
        axs[index].set_yticks(np.arange(len(amino_acids)))
        # ... and label them with the respective list entries
        axs[index].set_xticklabels(amino_acids)
        axs[index].set_yticklabels(amino_acids)
        # Loop over data dimensions and create text annotations.
        '''
        for i in range(len(amino_acids)):
            for j in range(len(amino_acids)):
                text = axs[index].text(j, i, "{0:.2f}".format(predicted[i, j] if index == 0 else estimated[i, j]),
                                       ha="center", va="center", color="w")
        '''

    filt = np.isfinite(predicted) & np.isfinite(estimated)
    x = predicted[filt].flatten()
    y = estimated[filt].flatten()
    axs[2].scatter(x, y)
    axs[2].set_title('Fixation probabilities', fontsize=font_size)

    model = sm.OLS(y, sm.add_constant(x))
    results = model.fit()
    b, a = results.params[0:2]
    idf = np.linspace(min(x), max(x), 100)
    axs[2].plot(idf, a * idf + b, '-',
                label=r"$y={0:.2g}x {3} {1:.2g}$ ($r^2={2:.2g})$".format(a, abs(b), results.rsquared,
                                                                         "+" if float(b) > 0 else "-"))
    axs[2].set_xlabel('Predicted ($\\phi$)', fontsize=font_size)
    axs[2].set_ylabel('Estimated ($\\widehat{\\omega}$)', fontsize=font_size)
    axs[2].legend(fontsize=font_size)
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
    parser.add_argument('-n', '--mutation_matrix', required=False, type=str, default="False", dest="mutation_matrix")

    args = parser.parse_args()

    if args.mutation_matrix != "False":
        df = pd.read_csv(args.mutation_matrix, sep="\t")
        nuc_matrix = np.zeros((len(nucleotides), len(nucleotides)))

        for id_s, nuc_source in enumerate(nucleotides):
            for id_t, nuc_target in enumerate(nucleotides):
                if nuc_source == nuc_target: continue
                sub = "q_" + nuc_source + nuc_target
                nuc_matrix[id_s, id_t] = df[sub].values[0]
                nuc_matrix[id_s, id_s] -= nuc_matrix[id_s, id_t]

        freqs = null_space(nuc_matrix.T).T[0]
        freqs /= freqs.sum()
        assert (np.sum(np.abs(np.dot(freqs, nuc_matrix))) < 1e-5)

    nested_dict = nested_dict_init()
    predicted_array, estimated_array = [], []
    for batch in args.input:
        at_gc_pct = float(batch.split("/")[-1].split("_")[0])

        if args.mutation_matrix == "False":
            freqs = [(at_gc_pct if nuc in weak_nucleotides else 1.0) for nuc in nucleotides]
            freqs = np.array(freqs) / sum(freqs)
        alpha = float(batch.split("/")[-2])
        exp = batch.replace(args.model + "_run.bf_hyout.txt", "exp")
        species, alignment = open_fasta_file(exp + (".ThirdPos.fasta" if args.model == "GTR" else ".fasta"))
        ali_dico = stats_from_ali(alignment)
        nested_dict["AT/GC_obs"][at_gc_pct] = ali_dico["at_over_gc"]

        hyphy_dico = dico_from_file(batch)
        format_hyphy_dico(hyphy_dico, args.model)

        if args.model == "MF":
            profile_path = "/".join(batch.split("/")[:-1]) + "_profile.prefs"
            predicted_array.append(omega_pairwise_from_profile(profile_path, freqs))
            estimated_array.append(omega_pairwise_from_hyphy(hyphy_dico))
            plot_pairwise_matrices(predicted_array[-1], estimated_array[-1],
                                   "{0}/omega.{1}".format(args.output, at_gc_pct))

        results_dico = {k: v[0] for k, v in pd.read_csv(exp + ".tsv", sep='\t').items()}
        nested_dict["w_obs"][at_gc_pct] = results_dico["dnd0_event_tot"]
        if args.model != "GTR":
            nested_dict["w_inf"][at_gc_pct] = hyphy_dico["w"]

        if ("pnG" in hyphy_dico) and ("pnG" in hyphy_dico):
            gc_pct = hyphy_dico["pnG"] + hyphy_dico["pnC"]
            nested_dict["lambda_inf"][at_gc_pct] = (1 - gc_pct) / gc_pct

        nested_dict["AT/GC_inf"][at_gc_pct] = equilibrium_lambda(hyphy_dico)

    if args.model == "MF":
        estimated_mean = np.mean(estimated_array, axis=0)
        predicted_mean = np.mean(predicted_array, axis=0)
        plot_pairwise_matrices(predicted_mean, estimated_mean, "{0}/mean.omega".format(args.output))

    x_list = sorted(nested_dict["AT/GC_obs"].keys())
    if args.model != "GTR":
        omega_obs = np.array([nested_dict["w_obs"][k] for k in x_list])
        omega_inf = np.array([nested_dict["w_inf"][k] for k in x_list])
        f = open("{0}/omega.txt".format(args.output), "w")
        f.write("|w_obs-w_inf|/w_obs = {0:.2f}%".format(100 * np.mean(np.abs((omega_inf - omega_obs)) / omega_obs)))
        f.close()

    if args.mutation_matrix != "False":
        exit(0)

    fig, ax = plt.subplots()
    ax.plot(x_list, x_list, color="black", linestyle='-', linewidth=2, label="y=x")

    list_plot = list()
    list_plot.append({"experiment": "lambda_inf", "color": GREEN, "linestyle": '--', "linewidth": 4,
                      "label": "$\\widehat{\\lambda}$ inferred"})
    list_plot.append(
        {"experiment": "AT/GC_obs", "color": BLUE, "linestyle": '-', "linewidth": 2, "label": "$AT/GC$ observed"})
    # list_plot.append({"experiment": "AT/GC_inf", "color": YELLOW, "linestyle": '--', "linewidth": 4, "label": "$\\widehat{AT/GC}$ inferred"})

    for param in list_plot:
        x_list = sorted(nested_dict[param["experiment"]].keys())
        y_list = [nested_dict[param["experiment"]][k] for k in x_list]
        ax.plot(x_list, y_list, linestyle=param["linestyle"], label=param["label"],
                color=param["color"], linewidth=param["linewidth"])

    lambda_obs = np.array(x_list)
    lambda_inf = np.array([nested_dict["lambda_inf"][k] for k in x_list])
    f = open("{0}/lambda.txt".format(args.output), "w")
    f.write("|lambda_obs-lambda_inf|/lambda_obs = {0:.2f}%".format(
        100 * np.mean(np.abs((lambda_inf - lambda_obs)) / lambda_obs)))
    f.close()
    ax.set_xscale('log')
    ax.set_xlabel('$\\lambda$ used for the simulation', fontsize=font_size)
    ax.set_yscale('log')
    ax.set_ylabel('$\\lambda$ estimated', fontsize=font_size)
    ax.get_xaxis().get_major_formatter().labelOnlyBase = False
    ax.legend(fontsize=font_size)
    model_name = {"GTR": "General time-reversible (GTR) on third positions", "MG": "Muse & Gaut codon model", "MF": "Mean-field codon model"}
    ax.set_title(model_name[args.model], fontsize=font_size)
    ax.set_xticks([0.2, 1, 5])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_yticks([0.2, 1, 5])
    ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt.xticks(fontsize=legend_size)
    plt.yticks(fontsize=legend_size)
    plt.tight_layout()
    plt.savefig("{0}/lambda.pdf".format(args.output), format="pdf")
    plt.savefig("{0}/lambda.png".format(args.output), format="png")
    plt.clf()
    plt.close('all')
