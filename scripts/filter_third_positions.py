# GLOBAL IMPORTS
import argparse
from codons import codon_table
import numpy as np
from stat_simulated import open_ali_file, open_fasta_file
from ete3 import Tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--input', required=True, type=str, dest="input")
    parser.add_argument('-t', '--tree', required=False, default="", type=str, dest="tree")
    parser.add_argument('-a', '--ancestral', required=False, default="", type=str, dest="ancestral")
    args = parser.parse_args()

    if ".ali" in args.input:
        species, alignment = open_ali_file(args.input)
        filt_fasta = open(args.input.replace(".ali", ".fasta"), 'w')
        for id_sp, sp in enumerate(species):
            filt_fasta.write(">{0}\n{1}\n".format(sp, alignment[id_sp]))
        filt_fasta.close()
    else:
        species, alignment = open_fasta_file(args.input)
    dict_ali = {sp: alignment[i] for i, sp in enumerate(species)}
    if args.tree:
        t = Tree(args.tree, format=1)
        anc_nodes = [n for n in t.iter_descendants() if n.name == args.ancestral]
        assert (len(anc_nodes) == 1)
        dict_ali = {sp: alignment[i] for i, sp in enumerate(species)}
        species = anc_nodes[0].get_leaf_names()
        alignment = [dict_ali[sp] for sp in species]
        anc_nodes[0].write(outfile=args.tree.replace(".nhx", "." + args.ancestral + ".nhx"), format=1)

        filt_fasta = open(args.output.replace(".fasta", ".codons.fasta"), 'w')
        for id_sp, sp in enumerate(species):
            filt_fasta.write(">{0}\n{1}\n".format(sp, dict_ali[sp]))
        filt_fasta.close()

    alignment = np.array([list(s) for s in alignment])
    filtered_sites = []
    for site in range(int(alignment.shape[1] / 3)):
        set_aa = set([codon_table["".join(c)] for c in alignment[:, 3 * site:3 * site + 3]])
        set_aa.discard("-")
        if len(set_aa) == 1:
            filtered_sites.append(alignment[:, (3 * site) + 2])

    filtered_sites = np.array(filtered_sites).T

    filt_fasta = open(args.output, 'w')
    for id_sp, sp in enumerate(species):
        filt_fasta.write(">{0}\n{1}\n".format(sp.strip(), "".join(filtered_sites[id_sp, :]).replace('?', '-')))
    filt_fasta.close()
