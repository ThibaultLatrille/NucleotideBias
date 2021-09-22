# GLOBAL IMPORTS
import argparse
import pandas as pd
from plot_module import *
from hyphy_format import *
import statsmodels.api as sm
from stat_simulated import stats_from_ali, open_fasta_file, omega_pairwise_from_profile
from scipy.linalg import null_space


def plot_pairwise_matrices(predicted, estimated, output, estimated_list=None):
    if estimated_list is None: estimated_list = []
    fig, axs = plt.subplots(1, 3, figsize=(16, 6))
    cbar = fig.colorbar(axs[0].imshow(predicted), ax=axs[0], orientation='horizontal', fraction=.05)
    cbar.ax.tick_params(labelsize=font_size)
    axs[0].set_title('$\\left\\langle 2 N_{\\mathrm{e}} \\mathbb{P}_{\\mathrm{fix}} \\right\\rangle $ '
                     'predicted between\npairs of amino-acids', fontsize=font_size * 1.2)
    cbar = fig.colorbar(axs[1].imshow(estimated), ax=axs[1], orientation='horizontal', fraction=.05)
    cbar.ax.tick_params(labelsize=font_size)
    axs[1].set_title('$\\widehat{\\omega}$ estimated between\npairs of amino-acids', fontsize=font_size * 1.2)

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
    axs[2].set_title('Fixation probabilities', fontsize=font_size * 1.2)

    model = sm.OLS(y, sm.add_constant(x))
    results = model.fit()
    b, a = results.params[0:2]
    idf = np.linspace(min(x), max(x), 100)
    axs[2].plot(idf, a * idf + b, '-', color=BLUE,
                label=r"$y={0:.2g}x {3} {1:.2g}$ ($r^2={2:.2g})$".format(a, abs(b), results.rsquared,
                                                                         "+" if float(b) > 0 else "-"))

    if estimated_list and len(estimated_list) > 5:
        # yerr = np.array([y - np.percentile(estimated_list, 5, axis=0)[filt].flatten(), np.percentile(estimated_list, 95, axis=0)[filt].flatten() - y])
        yerr = 1.96 * np.std(estimated_list, axis=0)[filt].flatten() / np.sqrt(len(estimated_list))
        axs[2].errorbar(x, y, yerr=yerr, fmt='o', marker=None, mew=0, ecolor=BLUE, lw=0.5,
                        zorder=-1)
    axs[2].set_xlabel('Predicted $\\left\\langle 2 N_{\\mathrm{e}} \\mathbb{P}_{\\mathrm{fix}} \\right\\rangle $',
                      fontsize=font_size * 1.2)
    axs[2].set_ylabel('Estimated $\\widehat{\\omega}$', fontsize=font_size * 1.2)
    axs[2].legend(fontsize=font_size * 0.8, loc="lower right")
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

    freqs = np.zeros(4) * 0.25
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

    nested_dict = defaultdict(lambda: defaultdict(lambda: list()))
    predicted_dico = {}
    estimated_array = []
    for batch in args.input:
        at_gc_pct = float(batch.split("/")[-1].split("_")[0])
        replicate = int(batch.split("/")[-1].split("_")[1])

        if args.mutation_matrix == "False":
            freqs = [(at_gc_pct if nuc in weak_nucleotides else 1.0) for nuc in nucleotides]
            freqs = np.array(freqs) / sum(freqs)
        alpha = float(batch.split("/")[-2])
        exp = batch.replace(args.model + "_run.bf_hyout.txt", "exp")
        species, alignment = open_fasta_file(exp + (".ThirdPos.fasta" if args.model == "GTR" else ".fasta"))
        ali_dico = stats_from_ali(alignment)
        nested_dict[at_gc_pct]["AT/GC_obs"].append(ali_dico["at_over_gc"])

        hyphy_dico = dico_from_file(batch)
        format_hyphy_dico(hyphy_dico, args.model)

        if args.model == "MF":
            profile_path = "/".join(batch.split("/")[:-1]) + "_profile.prefs"
            key = profile_path + str(at_gc_pct)
            if key not in predicted_dico:
                predicted_dico[key] = omega_pairwise_from_profile(profile_path, freqs)
            estimated_array.append(omega_pairwise_from_hyphy(hyphy_dico))
            plot_pairwise_matrices(predicted_dico[key], estimated_array[-1],
                                   "{0}/omega.{1}_{2}".format(args.output, at_gc_pct, replicate))

        results_dico = {k: v[0] for k, v in pd.read_csv(exp + ".tsv", sep='\t').items()}
        nested_dict[at_gc_pct]["w_obs"].append(results_dico["dnd0_event_tot"])
        if args.model != "GTR":
            nested_dict[at_gc_pct]["w_inf"].append(hyphy_dico["w"])

        if ("pnG" in hyphy_dico) and ("pnG" in hyphy_dico):
            gc_pct = hyphy_dico["pnG"] + hyphy_dico["pnC"]
            nested_dict[at_gc_pct]["lambda_inf"].append((1 - gc_pct) / gc_pct)

        nested_dict[at_gc_pct]["AT/GC_inf"].append(equilibrium_lambda(hyphy_dico))

    if args.model == "MF":
        predicted_mean = np.mean(list(predicted_dico.values()), axis=0)
        estimated_mean = np.mean(estimated_array, axis=0)
        plot_pairwise_matrices(predicted_mean, estimated_mean, "{0}/mean.omega".format(args.output), estimated_array)

    if args.model != "GTR":
        f = open("{0}/omega.tsv".format(args.output), "w")
        f.write("ω (precision in %)" + ("\tλ\n" if len(nested_dict) > 1 else "\n"))
        for lambda_mut in nested_dict.keys():
            omega_obs = np.array(nested_dict[lambda_mut]["w_obs"])
            omega_inf = np.array(nested_dict[lambda_mut]["w_inf"])
            f.write("{0:.2f}".format(100 * np.mean(np.abs((omega_inf - omega_obs)) / omega_obs)) + (
                "\t{0}".format(lambda_mut) if len(nested_dict) > 1 else ""))
        f.close()

    lambda_mut = list(nested_dict.keys())
    if len(lambda_mut) < 1: exit(0)

    fig, ax = plt.subplots()
    ax.plot(lambda_mut, lambda_mut, color="black", linestyle='-', linewidth=2, label="y=x")

    list_plot = list()
    list_plot.append({"experiment": "lambda_inf", "color": GREEN, "linestyle": '--', "linewidth": 4,
                      "label": "$\\widehat{\\lambda}$ inferred"})
    list_plot.append(
        {"experiment": "AT/GC_obs", "color": BLUE, "linestyle": '-', "linewidth": 2, "label": "AT/GC observed"})
    list_plot.append({"experiment": "AT/GC_inf", "color": YELLOW, "linestyle": '--', "linewidth": 4, "label": "$\\widehat{AT/GC}$ predicted"})

    for param in list_plot:
        y_list = np.array([np.mean(k[param["experiment"]]) for k in nested_dict.values()])
        ax.plot(lambda_mut, y_list, linestyle=param["linestyle"], label=param["label"],
                color=param["color"], linewidth=param["linewidth"])

        reps_set_len = set([len(k[param["experiment"]]) for k in nested_dict.values()])
        assert (len(reps_set_len) == 1)
        if reps_set_len.pop() > 5:
            lower = [np.percentile(k[param["experiment"]], 5) for k in nested_dict.values()]
            upper = [np.percentile(k[param["experiment"]], 95) for k in nested_dict.values()]
            ax.fill_between(lambda_mut, lower, upper, alpha=0.3, color=param["color"], facecolor=param["color"])

    lambda_inf = np.array([np.mean(k["lambda_inf"]) for k in nested_dict.values()])
    f = open("{0}/lambda.tsv".format(args.output), "w")
    f.write("λ (precision in %)\n{0:.2f}%".format(100 * np.mean(np.abs((lambda_inf - lambda_mut)) / lambda_mut)))
    f.close()
    ax.set_xscale('log')
    ax.set_xlabel('$\\lambda$ used for the simulation', fontsize=font_size)
    ax.set_yscale('log')
    ax.set_ylabel('$\\lambda$ estimated', fontsize=font_size)
    ax.get_xaxis().get_major_formatter().labelOnlyBase = False
    ax.legend(fontsize=font_size, loc="lower right")
    model_name = {"GTR": "General time-reversible (GTR) on third positions", "MG": "Muse & Gaut codon model",
                  "MF": "Mean-field codon model"}
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
