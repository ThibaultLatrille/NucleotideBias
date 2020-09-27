import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-a', '--alpha', required=True, type=float, dest="alpha")
    parser.add_argument('-n', '--nbr_sites', required=True, type=int, dest="nbr_sites")
    parser.add_argument('-s', '--seed', required=False, type=int, default=42, dest="seed")
    args = parser.parse_args()
    preferences = np.random.dirichlet(args.alpha * np.ones(20))
    prefs_file = open(args.output, 'w')
    np.random.seed(args.seed)
    prefs_file.write("site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y\n")
    for i in range(1, args.nbr_sites + 1):
        preferences = np.random.dirichlet(args.alpha * np.ones(20))
        prefs_file.write("{0},".format(i) + ",".join([str(i) for i in preferences]) + "\n")
    prefs_file.close()
