import numpy as np
import os

preferences = np.random.dirichlet(alpha*np.ones(20))
prefs_file = open(prefs_path, 'w')
prefs_file.write("# POSITION WT SITE_ENTROPY PI_A PI_C PI_D PI_E PI_F PI_G PI_H PI_I PI_K PI_L PI_M PI_N PI_P PI_Q PI_R PI_S PI_T PI_V PI_W PI_Y\n")
for i in range(1, nbr_sites + 1):
	preferences = np.random.dirichlet(alpha * np.ones(20))
	prefs_file.write("{0} A {0} ".format(i, alpha) + " ".join([str(i) for i in preferences]) + "\n")
prefs_file.close()