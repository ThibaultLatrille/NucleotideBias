# GLOBAL IMPORTS
import argparse
import pandas as pd
import os
from collections import Counter
from hyphy_format import *
from plot_module import *


def simpson_diversity(frequencies):
    return 1.0 / np.sum(frequencies * frequencies)


def shanon_diversity(frequencies):
    f = [i for i in frequencies if i != 0.0]
    return np.exp(-np.sum(f * np.log(f)))


def hill_diversity(frequencies, q):
    if q == 1:
        return shanon_diversity(frequencies)
    else:
        return np.power(np.sum(np.power(frequencies, q)), 1.0 / (1.0 - q))


def compute_diversity(frequencies_seq):
    n = frequencies_seq.shape[0]
    site_diversity = [shanon_diversity(frequencies_seq[site, :]) for site in range(n)]
    diversity = shanon_diversity(np.mean(frequencies_seq, axis=0))
    return np.mean(site_diversity), diversity


def open_ali_file(ali_path):
    ali_list = []
    with open(ali_path, 'r') as ali_file:
        next(ali_file)
        for line in ali_file:
            if line != "\n":
                _, seq = line.replace("  ", " ").replace("\n", "").split(" ")
                assert len(seq) % 3 == 0
                ali_list.append(seq)
    return ali_list


def open_fasta_file(fasta_path):
    ali_list = []
    with open(fasta_path, 'r') as fasta_file:
        for seq_unstriped in fasta_file:
            if seq_unstriped[0] != ">":
                seq = seq_unstriped.strip().replace("\n", "")
                assert len(seq) % 3 == 0
                ali_list.append(seq)
    return ali_list


def stats_from_ali(alignment):
    out_dico = {}
    assert (len(set([len(s) for s in alignment])) == 1)
    nbr_sites = int(len(alignment[0]) / 3)
    concat = "".join(alignment)
    count = Counter(concat)
    out_dico["at_over_gc"] = (count["A"] + count["T"]) / (count["C"] + count["G"])
    for pos in range(0, 3):
        count = Counter(concat[pos::3])
        out_dico["at_over_gc_{0}".format(pos + 1)] = (count["A"] + count["T"]) / (count["C"] + count["G"])
    out_dico["NbrSites"] = nbr_sites
    out_dico["NbrTaxa"] = len(alignment)

    aa_freqs = np.zeros((nbr_sites, len(amino_acids)))
    for site in range(nbr_sites):
        count = Counter("".join(codon_table[s[site * 3:(site + 1) * 3]] for s in alignment))
        assert ("X" not in count)
        for aa in amino_acids:
            aa_freqs[site, aa_char_to_int[aa]] = count[aa]
    aa_freqs /= len(alignment)

    out_dico["site_diversity_aa"], out_dico["diversity_aa"] = compute_diversity(aa_freqs)
    return out_dico


def y_label(p):
    if "at_over_gc" in p:
        return 'AT/GC'
    elif "diversity" in p:
        return "D"
    elif "omega_WS" in p:
        return "$\\omega_{\\mathrm{AT} \\rightarrow \\mathrm{GC}}$"
    elif "omega_SW" in p:
        return "$\\omega_{\\mathrm{GC} \\rightarrow \\mathrm{AT}}$"
    elif "omega" in p:
        return "$\\omega$"
    else:
        return p


def title(p):
    if "at_over_gc_" in p:
        return "Codon position " + p[-1]
    elif "site_diversity_aa" in p:
        return "Site-specific amino-acid diversity"
    elif "diversity_aa" in p:
        return "Sequence amino-acid diversity"
    elif "omega_WS" in p:
        return "Restricted to non-synonymous\nmutations from AT to GC"
    elif "omega_SW" in p:
        return "Restricted to non-synonymous\nmutations from GC to AT"
    elif "omega" in p:
        return "All non-synonymous mutations"
    else:
        return p


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--input', required=True, type=str, nargs='+', dest="input")
    args = parser.parse_args()

    dico_params = {"omega": {}, "omega_WS": {}, "omega_SW": {}}
    x_range, y_range = set(), set()
    for filepath in args.input:
        if not os.path.isfile(filepath): continue
        f_slash = filepath.split("/")
        x, y = float(f_slash[-1].split("_")[-2]), float(f_slash[-2])
        x_range.add(x)
        y_range.add(y)
        for param, val in stats_from_ali(open_ali_file(filepath + ".ali")).items():
            if param not in dico_params: dico_params[param] = dict()
            dico_params[param][(x, y)] = val
        results_dico = {k: v[0] for k, v in pd.read_csv(filepath + ".tsv", sep='\t').items()}
        dico_params["omega"][(x, y)] = results_dico["dnd0_event_tot"]
        dico_params["omega_WS"][(x, y)] = results_dico["dnd0_WS_event_tot"]
        dico_params["omega_SW"][(x, y)] = results_dico["dnd0_SW_event_tot"]

    x_range, y_range = sorted(x_range), sorted(y_range)
    dico_params.pop("NbrSites")
    dico_params.pop("NbrTaxa")
    plt.figure()

    for param, dico_values in dico_params.items():
        for y in y_range:
            plt.plot(x_range, [dico_values[(x, y)] for x in x_range], linewidth=2, label="$\\alpha={0:.3g}$".format(y))
        plt.xscale("log")
        plt.xlabel("$\\lambda$", size=font_size)
        plt.ylabel(y_label(param), size=font_size)
        plt.title(title(param), size=font_size)
        if "_aa" in param:
            plt.ylim((1, 20))
        elif "omega" in param:
            plt.ylim((0, 1))
        elif "at_over_gc" in param:
            plt.yscale("log")
        plt.legend(fontsize=font_size)
        plt.xticks(fontsize=font_size)
        plt.yticks(fontsize=font_size)
        plt.tight_layout()
        plt.savefig("{0}/{1}.pdf".format(args.output, param), format="pdf")
        plt.savefig("{0}/{1}.png".format(args.output, param), format="png")
        plt.clf()
    plt.close('all')
