import matplotlib

matplotlib.rcParams["font.family"] = ["Latin Modern Sans"]
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler

plt.rc('axes', prop_cycle=cycler(color=["#5D80B4", "#E29D26", "#8FB03E", "#EB6231", "#857BA1"]))
font_size = 18
legend_size = 14