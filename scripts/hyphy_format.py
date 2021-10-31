# GLOBAL IMPORTS
from codons import *


def dico_from_file(filename):
    tmp_dico = {}
    tmp_file = open(filename, "r")
    for line in tmp_file:
        for sub_line in line.split(";"):
            split_line = sub_line.split("=")
            if len(split_line) <= 1: continue
            key = split_line[0].replace("global ", "").strip()
            if key in tmp_dico: continue
            value = split_line[-1].strip()
            try:
                tmp_dico[key] = float(value)
            except:
                pass

    tmp_file.close()
    return tmp_dico


def format_hyphy_dico(hyphy_dico, model):
    df = {"GTR": 6 + 3, "MG": 6 + 3 + 1, "MF": 6 + 3 + 75 + 19}
    if model in df.keys() and "Log Likelihood" in hyphy_dico:
        hyphy_dico["df"] = df[model]
        hyphy_dico["AIC"] = 2 * df[model] - 2 * hyphy_dico["Log Likelihood"]

    if "pnCG" in hyphy_dico:
        hyphy_dico["pnC"] = hyphy_dico["pnCG"]
        hyphy_dico["pnG"] = hyphy_dico["pnCG"]
        hyphy_dico["pnA"] = hyphy_dico["pnAT"]
    hyphy_dico["pnT"] = 1 - (hyphy_dico["pnC"] + hyphy_dico["pnG"] + hyphy_dico["pnA"])

    if model in ["GTR", "MG"]:
        return hyphy_dico

    if "epsA" in hyphy_dico and "epsM" not in hyphy_dico:
        hyphy_dico["epsM"] = 20 - sum([hyphy_dico["eps" + aa] for aa in amino_acids_set if aa != "M"])

    codon_frequencies = np.ones(len(codons))
    for codon_index, codon in enumerate(codons):
        for nuc in codon:
            codon_frequencies[codon_index] *= hyphy_dico["pn" + nuc]
        epsilon = "eps" + codon_table[codon]
        if epsilon in hyphy_dico:
            codon_frequencies[codon_index] *= hyphy_dico[epsilon]

    codon_frequencies /= np.sum(codon_frequencies)

    d_dict, d0_dict = defaultdict(float), defaultdict(float)
    d, d0 = 0.0, 0.0

    for x, codon_origin in enumerate(codons):
        for y, a, b in non_syn_neighbors[x]:
            nuc_origin, nuc_target = nucleotides[a], nucleotides[b]
            codon_target = codons[y]

            mut_flow_tmp = codon_frequencies[x] * hyphy_dico["pn" + nuc_target]

            exchan = "exch" + "".join(sorted((nuc_target, nuc_origin)))
            if exchan in hyphy_dico:
                mut_flow_tmp *= hyphy_dico[exchan]

            p_fix = 1.0
            beta = "b_" + "".join(sorted((codon_table[codon_origin], codon_table[codon_target])))
            if beta in hyphy_dico:
                assert (
                    aa_neighbors[aa_char_to_int[codon_table[codon_origin]], aa_char_to_int[codon_table[codon_target]]])
                if hyphy_dico[beta] > 100.0:
                    print("{0}={1}".format(beta, hyphy_dico[beta]))
                p_fix *= hyphy_dico[beta]

            epsilon = "eps" + codon_table[codon_target]
            if epsilon in hyphy_dico:
                p_fix *= hyphy_dico[epsilon]

            d += mut_flow_tmp * p_fix
            d0 += mut_flow_tmp

            omega_subset = "w_" + weak_strong(nuc_origin) + weak_strong(nuc_target)
            d_dict[omega_subset] += mut_flow_tmp * p_fix
            d0_dict[omega_subset] += mut_flow_tmp

    hyphy_dico["w"] = d / d0
    for omega_subset in d_dict.keys():
        if d_dict[omega_subset] != 0.0:
            hyphy_dico[omega_subset] = d_dict[omega_subset] / d0_dict[omega_subset]

    return hyphy_dico


def equilibrium_lambda(hyphy_dico):
    at_pct_tot, at_pct_1, at_pct_2, at_pct_3 = equilibrium_at_pct(hyphy_dico)
    return at_pct_tot / (1 - at_pct_tot), at_pct_1 / (1 - at_pct_1), at_pct_2 / (1 - at_pct_2), at_pct_3 / (
                1 - at_pct_3)


def equilibrium_at_pct(hyphy_dico):
    nbr_weak_tot = np.array([c.count("A") + c.count("T") for c in codons])
    nbr_weak_1 = np.array([c[0].count("A") + c[0].count("T") for c in codons])
    nbr_weak_2 = np.array([c[1].count("A") + c[1].count("T") for c in codons])
    nbr_weak_3 = np.array([c[2].count("A") + c[2].count("T") for c in codons])
    codon_frequencies = np.ones(len(codons))

    for index, codon in enumerate(codons):
        for nuc in codon:
            codon_frequencies[index] *= hyphy_dico["pn" + nuc]
        epsilon = "eps" + codon_table[codon]
        if epsilon in hyphy_dico:
            codon_frequencies[index] *= hyphy_dico[epsilon]

    codon_frequencies /= np.sum(codon_frequencies)

    return np.sum(codon_frequencies * nbr_weak_tot) / 3, np.sum(codon_frequencies * nbr_weak_1), np.sum(
        codon_frequencies * nbr_weak_2), np.sum(codon_frequencies * nbr_weak_3)


def omega_pairwise_from_hyphy(hyphy_dico):
    omega_array = np.ones((len(amino_acids), len(amino_acids)))
    omega_array[:] = np.NaN

    for aa_source, source in aa_char_to_int.items():
        for aa_target, target in aa_char_to_int.items():
            beta = "b_" + "".join(sorted((aa_source, aa_target)))
            eps = "eps" + aa_target
            if (beta in hyphy_dico) and (eps in hyphy_dico):
                omega_array[source, target] = hyphy_dico[beta] * hyphy_dico[eps]
    return omega_array
