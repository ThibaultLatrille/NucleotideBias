from pandas import DataFrame
from codons import nucleotides, weak_nucleotides
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-l', '--lambda', required=True, type=float, dest="mut_bias")
    args = parser.parse_args()
    out = {}

    freqs = np.array([args.mut_bias if t in weak_nucleotides else 1.0 for t in nucleotides])
    freqs /= freqs.sum()
    nuc_matrix = np.zeros((len(nucleotides), len(nucleotides)))
    events = 0
    for id_s, s in enumerate(nucleotides):
        for id_t, t in enumerate(nucleotides):
            if s == t: continue
            sub = "q_" + s + t
            nuc_matrix[id_s, id_t] = args.mut_bias if t in weak_nucleotides else 1.0

        nuc_matrix[id_s, id_s] = -nuc_matrix[id_s, :].sum()

    assert (np.sum(np.abs(np.dot(freqs, nuc_matrix))) < 1e-5)
    for id_s, s in enumerate(nucleotides):
        for id_t, t in enumerate(nucleotides):
            if s == t: continue
            out["q_" + s + t] = [nuc_matrix[id_s, id_t]]
    DataFrame(out).to_csv(args.output, index=False, sep="\t")
