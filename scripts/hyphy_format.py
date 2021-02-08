# GLOBAL IMPORTS
from codons import *


def dico_from_file(filename):
    tmp_dico = {}
    tmp_file = open(filename, "r")
    for line in tmp_file:
        for sub_line in line.split(";"):
            split_line = sub_line.split("=")
            if len(split_line) > 1:
                value = split_line[-1].strip()
                try:
                    tmp_dico[split_line[0]] = float(value)
                except:
                    pass
    tmp_file.close()
    return tmp_dico


def format_hyphy_dico(hyphy_dico):
    if "pnCG" in hyphy_dico:
        hyphy_dico["pnC"] = hyphy_dico["pnCG"]
        hyphy_dico["pnG"] = hyphy_dico["pnCG"]
        hyphy_dico["pnA"] = hyphy_dico["pnAT"]
        hyphy_dico["pnT"] = hyphy_dico["pnAT"]

    if 'w' not in hyphy_dico:
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
                    if hyphy_dico[beta] > 100.0:
                        print("{0}={1}".format(beta, hyphy_dico[beta]))
                        hyphy_dico[beta] = 1.0
                    p_fix *= hyphy_dico[beta]

                omega_subset = "w_" + weak_strong(nuc_origin) + weak_strong(nuc_target)
                if omega_subset in hyphy_dico:
                    p_fix *= hyphy_dico[omega_subset]

                epsilon = "eps" + codon_table[codon_target]
                if epsilon in hyphy_dico:
                    p_fix *= hyphy_dico[epsilon]

                d += mut_flow_tmp * p_fix
                d0 += mut_flow_tmp

                if omega_subset not in hyphy_dico:
                    d_dict[omega_subset] += mut_flow_tmp * p_fix
                    d0_dict[omega_subset] += mut_flow_tmp

        hyphy_dico["w"] = d / d0
        for omega_subset in d_dict.keys():
            if d_dict[omega_subset] != 0.0 and (omega_subset not in hyphy_dico):
                hyphy_dico[omega_subset] = d_dict[omega_subset] / d0_dict[omega_subset]

    return hyphy_dico


def equilibrium_lambda(hyphy_dico):
    at_pct = equilibrium_at_pct(hyphy_dico)
    return at_pct / (1 - at_pct)


def equilibrium_at_pct(hyphy_dico):
    nbr_weak = np.array([c.count("A") + c.count("T") for c in codons])
    codon_frequencies = np.ones(len(codons))

    for index, codon in enumerate(codons):
        for nuc in codon:
            codon_frequencies[index] *= hyphy_dico["pn" + nuc]
        epsilon = "eps" + codon_table[codon]
        if epsilon in hyphy_dico:
            codon_frequencies[index] *= hyphy_dico[epsilon]

    codon_frequencies /= np.sum(codon_frequencies)

    return np.sum(codon_frequencies * nbr_weak) / 3
