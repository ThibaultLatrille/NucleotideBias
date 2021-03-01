# GLOBAL IMPORTS
import argparse
import pandas as pd
import numpy as np
from scipy.stats.distributions import chi2


def format_float(x):
    try:
        return [format_float_v(y) for y in x]
    except:
        return format_float_v(x)


def format_float_v(x):
    if 0.001 < abs(x) < 10:
        return "{:6.3f}".format(x)
    elif 0.001 < abs(x) < 10**6:
        return "{:6.1f}".format(x)
    else:
        s = "{:6.2g}".format(x)
        if "e" in s:
            mantissa, exp = s.split('e')
            s = mantissa + 'e$^{' + str(int(exp)) + '}$'
            s = " " * (5 + 6 - len(s)) + s
        return s


def tex_float(x):
    s = "{0:.3g}".format(x)
    if "e" in s:
        mantissa, exp = s.split('e')
        return mantissa + 'e^{' + exp + '}'
    else:
        return s


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, dest="output")
    parser.add_argument('-i', '--input', required=True, type=str, nargs='+', dest="input")
    args = parser.parse_args()

    concat = pd.concat([pd.read_csv(filepath, sep='\t') for filepath in args.input], axis=0, ignore_index=True)
    if "AIC_MG" in concat and "AIC_MF" in concat:
        concat["Δ_AIC"] = concat["AIC_MG"] - concat["AIC_MF"]
        concat["weight_MF"] = 1 / (np.exp(-concat["Δ_AIC"] / 2) + 1)
        concat["weight_MG"] = np.exp(-concat["Δ_AIC"] / 2) / (np.exp(-concat["Δ_AIC"] / 2) + 1)
        assert (concat["df_MF"] - concat["df_MG"] > 0).all()
        concat["LRT_MF_MG"] = -2 * (concat["LnL_MG"] - concat["LnL_MF"])
        assert (concat["LRT_MF_MG"] > 0).all()
        concat["Chi2_LRT_MF_MG"] = chi2.sf(concat["LRT_MF_MG"], concat["df_MF"] - concat["df_MG"])
    concat.to_csv(args.output, float_format=format_float, sep="\t", index=False, columns=list(sorted(concat)))
    concat.to_latex(args.output.replace(".tsv", '.tex'), float_format=format_float, index=False,
                    columns=list(sorted(concat)))

    header = list(concat["name"].values)
    concat.pop("name")
    concat_t = concat.transpose()
    concat_t.to_csv(args.output.replace(".tsv", '.T.tsv'), float_format=format_float, sep="\t", header=header)
    concat_t.to_latex(args.output.replace(".tsv", '.T.tex'), float_format=format_float, header=header)
