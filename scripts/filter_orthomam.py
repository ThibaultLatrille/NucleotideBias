#!python3
import argparse
import pandas as pd
import numpy as np
from collections import Counter
from stat_simulated import open_ali_file
import os
from codons import codon_table


def at_gc(filepath, third_pos):
    species, alignment = open_ali_file(filepath)
    size = len(alignment[0])
    if third_pos:
        ali_array = np.array([list(s) for s in alignment])
        alignment = []
        for site in range(int(ali_array.shape[1] / 3)):
            set_aa = set([codon_table["".join(c)] for c in ali_array[:, 3 * site:3 * site + 3]])
            set_aa.discard("-")
            if len(set_aa) == 1:
                alignment.append("".join(ali_array[:, (3 * site) + 2]))

    at = 0.0
    gc = 0.0
    for seq in alignment:
        count = Counter(seq)
        at += count["A"] + count["T"]
        gc += count["C"] + count["G"]
    return at / gc, size


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, default='../OrthoMam/singlegene_alignments', type=str,
                        dest="input")
    parser.add_argument('-o', '--only_third_pos', required=False, default=False, type=bool, dest="only_third_pos")
    args = parser.parse_args()
    files = os.listdir(args.input)
    ali_at_gc = {}
    for n, file in enumerate(files):
        ali_at_gc[file] = at_gc(args.input + "/" + file, args.only_third_pos)
        print("{0:.2f}% completed ({1})".format(100 * n / len(files), file))
    sorted_ali = list(sorted(ali_at_gc, key=lambda k: ali_at_gc[k][0], reverse=True))
    df = pd.DataFrame({"CDS": sorted_ali, "AT/GC": [ali_at_gc[k][0] for k in sorted_ali],
                       "Sites": [ali_at_gc[k][1] for k in sorted_ali]})
    df.to_csv("../OrthoMam/AT_GC{0}.tsv".format("_4F" if args.only_third_pos else ""), index=False, sep="\t")
