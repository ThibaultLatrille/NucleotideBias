#!python3
import os
import argparse
import pandas as pd
from stat_simulated import open_ali_file
from ete3 import Tree

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--ali_folder', required=False, default='OrthoMam/singlegene_alignments', type=str, dest="ali_folder")
    parser.add_argument('-i', '--input', required=False, default='OrthoMam/AT_GC_4F.tsv', type=str, dest="input")
    parser.add_argument('-t', '--tree', required=False, default="OrthoMam/rootedtree.nhx", type=str, dest="tree")
    parser.add_argument('-o', '--output', required=False, default='DataEmpirical/OrthoMamPrimatesHighGC-4F/', type=str,
                        dest="output")
    parser.add_argument('-a', '--ancestral', required=False, default="GOPMCNPPHGRPCCMPMCSCCA", type=str, dest="ancestral")
    parser.add_argument('-s', '--nbr_sites', required=False, default=15000, type=int, dest="nbr_sites")
    parser.add_argument('-r', '--reversed', required=False, default=True, type=bool, dest="reversed")

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    ancestral = ''
    t = Tree(args.tree, format=1)
    anc_nodes = [n for n in t.iter_descendants() if n.name == args.ancestral]
    assert (len(anc_nodes) == 1)
    merge_alignment = {k: "" for k in anc_nodes[0].get_leaf_names()}
    anc_nodes[0].write(outfile="{0}/tree.newick".format(args.output), format=1)

    args = parser.parse_args()
    df = pd.read_csv(args.input, sep="\t")
    sites = 0

    for i, row in reversed(list(df.iterrows())) if args.reversed else df.iterrows():
        sites += row["Sites"]
        if sites > args.nbr_sites: break
        species, alignment = open_ali_file(args.ali_folder + "/" + row["CDS"])
        dict_ali = {sp: alignment[i] for i, sp in enumerate(species) if sp in merge_alignment}
        for taxon, seq_merge in merge_alignment.items():
            merge_alignment[taxon] += dict_ali[taxon] if taxon in dict_ali else "-" * row["Sites"]

    filt_fasta = open("{0}/alignment.fasta".format(args.output), 'w')
    for sp, seq in merge_alignment.items():
        filt_fasta.write(">{0}\n{1}\n".format(sp, seq.replace('?', '-')))
    filt_fasta.close()
