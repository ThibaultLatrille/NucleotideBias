# GLOBAL IMPORTS
import argparse
import pandas as pd
import os
from collections import Counter
from hyphy_format import *
from plot_module import *
from codons import amino_acids


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
    assert (os.path.isfile(ali_path))
    species_list, ali_list = [], []
    with open(ali_path, 'r') as ali_file:
        next(ali_file)
        for line in ali_file:
            if line != "\n":
                spe, seq = line.replace("  ", " ").replace("\t", " ").replace("\n", "").split(" ")
                ali_list.append(seq)
                species_list.append(spe)
    return species_list, ali_list


def open_fasta_file(fasta_path):
    species_list, ali_list = [], []
    with open(fasta_path, 'r') as fasta_file:
        for seq_unstriped in fasta_file:
            if seq_unstriped[0] != ">":
                seq = seq_unstriped.strip().replace("\n", "")
                ali_list.append(seq)
            else:
                species_list.append(seq_unstriped.replace(">", ""))
    return species_list, ali_list


def omega_pairwise_from_profile(profile_file, nuc_freqs):
    assert (os.path.isfile(profile_file))
    df = pd.read_csv(profile_file, sep=",")
    profiles = df.drop('site', axis=1).values
    fitnesses = np.log(profiles)
    inv_profiles = 1.0 / profiles
    nbr_sites = len(profiles)
    codon_frequencies = np.ones((nbr_sites, len(codons)))
    for site in range(nbr_sites):
        for codon_index, codon in enumerate(codons):
            for nuc in codon:
                codon_frequencies[site, codon_index] *= nuc_freqs[nucindex[nuc]]
            codon_frequencies[site, codon_index] *= profiles[site, aa_char_to_int[codon_table[codon]]]

    z_sites = 1.0 / codon_frequencies.sum(axis=1)

    omega_array = np.zeros((len(amino_acids), len(amino_acids)))
    omega_array[:] = np.NaN

    for source, aa_source in enumerate(amino_acids):
        for target, aa_target in enumerate(amino_acids):
            if not aa_neighbors[source, target]: continue

            omega_array[source, target] = 0.0
            for site in range(nbr_sites):
                if abs(fitnesses[site, target] - fitnesses[site, source]) < 1e-3:
                    omega_array[source, target] += z_sites[site] * (
                                1 + fitnesses[site, target] - fitnesses[site, source])
                else:
                    omega_array[source, target] += z_sites[site] * (
                                fitnesses[site, target] - fitnesses[site, source]) / (
                                                           inv_profiles[site, source] - inv_profiles[site, target])

            omega_array[source, target] /= np.sum([z_sites[site] * profiles[site, source] for site in range(nbr_sites)])
    return omega_array


def stats_from_ali(alignment):
    out_dico = {}
    assert (len(set([len(s) for s in alignment])) == 1)

    concat = "".join(alignment)
    count = Counter(concat)
    out_dico["at_over_gc"] = (count["A"] + count["T"]) / (count["C"] + count["G"])

    nbr_sites = int(len(alignment[0]) / 3)
    if nbr_sites % 3 != 0: return out_dico
    for pos in range(0, 3):
        count = Counter(concat[pos::3])
        out_dico["at_over_gc_{0}".format(pos + 1)] = (count["A"] + count["T"]) / (count["C"] + count["G"])
    out_dico["NbrSites"] = nbr_sites
    out_dico["NbrTaxa"] = len(alignment)

    aa_freqs = np.zeros((nbr_sites, len(amino_acids)))
    for site in range(nbr_sites):
        count = Counter("".join(codon_table[s[site * 3:(site + 1) * 3]] for s in alignment))
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
    elif "omega_WS_over_SW" in p:
        return "$\\phi_{\\mathrm{AT} \\rightarrow \\mathrm{GC}} / \\phi_{\\mathrm{GC} \\rightarrow \\mathrm{AT}}$"
    elif "omega_SW_over_WS" in p:
        return "$\\phi_{\\mathrm{GC} \\rightarrow \\mathrm{AT}} / \\phi_{\\mathrm{AT} \\rightarrow \\mathrm{GC}}$"
    elif "omega_WS" in p:
        return "$\\phi_{\\mathrm{AT} \\rightarrow \\mathrm{GC}}$"
    elif "omega_SW" in p:
        return "$\\phi_{\\mathrm{GC} \\rightarrow \\mathrm{AT}}$"
    elif "omega" in p:
        return "$\\phi$"
    else:
        return p


def title(p):
    if "at_over_gc_" in p:
        return "Codon position " + p[-1]
    elif "site_diversity_aa" in p:
        return "Site-specific amino-acid diversity"
    elif "diversity_aa" in p:
        return "Sequence amino-acid diversity"
    elif "omega_" in p and "over" in p:
        return "Restricted to non-synonymous mutations\nfrom AT to GC and GC to AT"
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

    dico_params = {"omega": {}, "omega_WS": {}, "omega_SW": {}, "omega_SW_over_WS": {}, "omega_WS_over_SW": {}}
    x_range, y_range = set(), set()
    for filepath in args.input:
        if not os.path.isfile(filepath): continue
        f_slash = filepath.split("/")
        x, y = float(f_slash[-1].split("_")[-2]), float(f_slash[-2])
        x_range.add(x)
        y_range.add(y)
        species, alignment = open_ali_file(filepath + ".ali")
        for param, val in stats_from_ali(alignment).items():
            if param not in dico_params: dico_params[param] = dict()
            dico_params[param][(x, y)] = val
        results_dico = {k: v[0] for k, v in pd.read_csv(filepath + ".tsv", sep='\t').items()}
        dico_params["omega"][(x, y)] = results_dico["dnd0_event_tot"]
        dico_params["omega_WS"][(x, y)] = results_dico["dnd0_WS_event_tot"]
        dico_params["omega_SW"][(x, y)] = results_dico["dnd0_SW_event_tot"]
        dico_params["omega_SW_over_WS"][(x, y)] = results_dico["dnd0_SW_event_tot"] / results_dico["dnd0_WS_event_tot"]
        dico_params["omega_WS_over_SW"][(x, y)] = results_dico["dnd0_WS_event_tot"] / results_dico["dnd0_SW_event_tot"]

    x_range, y_range = sorted(x_range), sorted(y_range)
    dico_params.pop("NbrSites")
    dico_params.pop("NbrTaxa")

    for param, dico_values in dico_params.items():
        fig, ax = plt.subplots()
        for y in y_range:
            label = "$N_{\\mathrm{r}}" + "={0:.3g}".format(y) + "$"
            ax.plot(x_range, [dico_values[(x, y)] for x in x_range], linewidth=2, label=label)
        ax.set_xscale("log")
        ax.set_xticks([0.2, 1, 5])
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.set_xlabel("$\\lambda$", size=font_size)
        ax.set_ylabel(y_label(param), size=font_size)
        ax.set_title(title(param), size=font_size * 1.2)
        if "_aa" in param:
            ax.set_ylim((1, 20))
        elif ("omega" in param) and ("over" not in param):
            ax.set_ylim((0, 1))
        elif "at_over_gc" in param:
            ax.set_yscale("log")
            ax.set_yticks([0.2, 1, 5])
            ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.legend(loc='upper left', fontsize=legend_size)
        plt.xticks(fontsize=legend_size)
        plt.yticks(fontsize=legend_size)
        plt.tight_layout()
        plt.savefig("{0}/{1}.pdf".format(args.output, param), format="pdf")
        plt.savefig("{0}/{1}.png".format(args.output, param), format="png")
        plt.clf()
    plt.close('all')
