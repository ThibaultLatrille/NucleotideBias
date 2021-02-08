#!python3
from ete3 import Tree
import argparse


def is_float(x):
    try:
        float(x.replace('\'', ""))
        return True
    except ValueError:
        return False


def tree_plot(input_tree):
    t = Tree(input_tree, format=1)
    names = set()
    for node in t.traverse():
        if not node.name or is_float(node.name):
            if node.is_root():
                name = "Root"
            else:
                leaves = node.get_leaf_names()
                name = "SEQ" + "".join([i.replace("SEQ", "") for i in leaves])
                while name in names:
                    name += "Bis"
            names.add(name)
            node.name = name
        print(node.name)
    t.write(format=1, outfile=input_tree + ".annotated")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--tree', required=False, type=str,
                        default='../DataEmpirical/ha.newick', dest="t", metavar="<tree>",
                        help="The tree to be re-written")
    args = parser.parse_args()
    tree_plot(args.t)
