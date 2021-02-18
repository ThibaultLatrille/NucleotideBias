from pandas import DataFrame
from hyphy_format import *
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--input', required=True, type=str, dest="input")
    args = parser.parse_args()

    hyphy_dico = dico_from_file(args.input)
    format_hyphy_dico(hyphy_dico, model="GTR")
    out = {}
    for s in nucleotides:
        for t in nucleotides:
            if s == t: continue
            out["q_" + s + t] = hyphy_dico["pn" + t]
            ex = "exch" + "".join(sorted((s, t)))
            if ex in hyphy_dico:
                out["q_" + s + t] *= hyphy_dico[ex]
    DataFrame({k: [v] for k, v in out.items()}).to_csv(args.output, index=False, sep="\t")
