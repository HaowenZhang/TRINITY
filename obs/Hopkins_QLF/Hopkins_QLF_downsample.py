import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from scipy.interpolate import UnivariateSpline


zs = [0.1, 0.5, 1.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

for z in zs:
	QLF_H = np.loadtxt('./Hopkins_z%.1f.qlf' % z)
	Lbol_grid = QLF_H[:,0]
	ND_grid = QLF_H[:,4]
	interp_H = UnivariateSpline(Lbol_grid, ND_grid, s=0)

	QLF_U = np.loadtxt('../Ueda_QLF/ueda_z%.1f.qlf' % z, comments='#')
	Lbol_down = (90 - QLF_U[:,0]) / 2.5
	ND_down = interp_H(Lbol_down)
	err_h = np.array([0.3] * len(ND_down))
	err_l = np.array([0.3] * len(ND_down))
	data = np.array([QLF_U[:,0], np.log10(ND_down), err_h, err_l]).transpose()

	np.savetxt('./Hopkins_z%.1f_downsample.qlf' % z, data, header='zlow: %.1f\nzhigh: %.1f\ntype: qlf\nref: Hopkins et al. 2007\nLog(Lx) Log(Nd)' % (z, z), comments='#')
