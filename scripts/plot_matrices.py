#!python3
import argparse
from codons import *
from os.path import basename
import numpy as np
import pandas as pd
from plot_module import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--input', required=True, type=str, nargs='+', dest="input")
    args = parser.parse_args()

    fig, axs = plt.subplots(1, len(args.input), figsize=(len(args.input) * 4, 4))

    model_nuc = {}
    for file in args.input:
        df = pd.read_csv(file, sep="\t")
        nuc_matrix = np.ones((len(nucleotides), len(nucleotides)))
        nuc_matrix[:] = np.NaN

        for id_nuc_source, nuc_source in enumerate(nucleotides):
            for id_nuc_target, nuc_target in enumerate(nucleotides):
                if nuc_source == nuc_target: continue
                sub = "q_" + nuc_source + nuc_target
                nuc_matrix[id_nuc_source, id_nuc_target] = df[sub].values[0]

        model_nuc[basename(file).split("_")[0]] = nuc_matrix

    vmin = min([np.min(n[np.isfinite(n)]) for n in model_nuc.values()])
    vmax = max([np.max(n[np.isfinite(n)]) for n in model_nuc.values()])

    for index, (model, nuc_matrix) in enumerate(model_nuc.items()):

        fig.colorbar(axs[index].imshow(nuc_matrix, vmin=vmin, vmax=vmax), ax=axs[index], orientation='horizontal',
                     fraction=.1)
        axs[index].set_title(label=model)

        # We want to show all ticks...
        axs[index].set_xticks(np.arange(len(nucleotides)))
        axs[index].set_yticks(np.arange(len(nucleotides)))
        # ... and label them with the respective list entries
        axs[index].set_xticklabels(nucleotides)
        axs[index].set_yticklabels(nucleotides)

        # Loop over data dimensions and create text annotations.
        for i in range(len(nucleotides)):
            for j in range(len(nucleotides)):
                text = axs[index].text(j, i, "{0:.2f}".format(nuc_matrix[i, j]), ha="center", va="center", color="w")

    fig.tight_layout()
    plt.savefig(args.output, format="png")
    plt.savefig(args.output.replace(".png", ".pdf"), format="pdf")
    plt.clf()
    plt.close('all')
