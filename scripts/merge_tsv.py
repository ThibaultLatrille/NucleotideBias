# GLOBAL IMPORTS
import argparse
import pandas as pd


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--input', required=True, type=str, nargs='+', dest="input")
    args = parser.parse_args()

    concat = pd.concat([pd.read_csv(filepath, sep='\t') for filepath in args.input], axis=0, ignore_index=True)
    concat.to_csv(args.output, sep="\t", index=False)
    concat.to_latex(args.output.replace(".tsv", '.tex'), index=False)

    header = list(concat["name"].values)
    concat.pop("name")
    concat_t = concat.transpose()
    concat_t.to_csv(args.output.replace(".tsv", '.T.tsv'), sep="\t", header=header)
    concat_t.to_latex(args.output.replace(".tsv", '.T.tex'), header=header)
