import numpy as np
import argparse
import pandas as pd

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--input', required=True, type=str, dest="input")
    parser.add_argument('-n', '--nbr_exons', required=True, type=int, dest="nbr_exons")
    parser.add_argument('-p', '--relative_pop_size', required=False, type=float, default=False,
                        dest="relative_pop_size")
    args = parser.parse_args()
    preferences = pd.read_csv(args.input, sep=",")
    profiles = np.power(preferences.drop('site', axis=1).values, args.relative_pop_size)
    for site in range(profiles.shape[0]):
        profiles[site, :] /= np.sum(profiles[site, :])
        assert(abs(np.sum(profiles[site, :]) - 1.0) < 1e-5)
    df = pd.DataFrame(np.concatenate([profiles] * args.nbr_exons, axis=0), columns=preferences.columns[1:])
    a = range(1, profiles.shape[0] + 1)
    df["site"] = range(1, profiles.shape[0] * args.nbr_exons + 1)
    df.to_csv(args.output, sep=',', encoding='utf-8', index=False)
