from pandas import DataFrame
from codons import nucleotides, weak_nucleotides
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-l', '--lambda', required=True, type=float, dest="mut_bias")
    args = parser.parse_args()

    out = {}
    for s in nucleotides:
        for t in nucleotides:
            if s == t: continue
            out["q_" + s + t] = [args.mut_bias if t in weak_nucleotides else 1.0]
    DataFrame(out).to_csv(args.output, index=False, sep="\t")
