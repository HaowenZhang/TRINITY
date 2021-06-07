import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

Lsun = np.log10(3.9) + 33

zs = ['1.0', '2.0', '3.0', '4.0', '5.0']
for z in zs:
	try:
		df = pd.read_csv('./hopkins_z' + z + '_EL.qlf', names=['Mi', 'nd', 'err_h', 'err_l'], comment='#',\
				delim_whitespace=True, dtype=float, skiprows=5)
		# Mi = np.array(df['Mi'])
		# nd = np.array
		plt.errorbar((90 - df['Mi']) / 2.5 - Lsun, df['nd'], yerr=[df['err_h'], df['err_l']], capsize=2, linewidth=0, elinewidth=0.5, marker='o',\
				label='Hopkins, emission line', color='g')
	except IOError:
		pass

	try:
		df = pd.read_csv('./hopkins_z' + z + '_softX.qlf', names=['Mi', 'nd', 'err_h', 'err_l'], comment='#',\
				delim_whitespace=True, dtype=float, skiprows=5)
		plt.errorbar((90 - df['Mi']) / 2.5 - Lsun, df['nd'], yerr=[df['err_h'], df['err_l']], capsize=2, linewidth=0, elinewidth=0.5, marker='+',\
				label='Hopkins, soft X-ray', color='b')
	except IOError:
		pass

	try:
		df = pd.read_csv('./hopkins_z' + z + '_hardX.qlf', names=['Mi', 'nd', 'err_h', 'err_l'], comment='#',\
				delim_whitespace=True, dtype=float, skiprows=5)
		plt.errorbar((90 - df['Mi']) / 2.5 - Lsun, df['nd'], yerr=[df['err_h'], df['err_l']], capsize=2, linewidth=0, elinewidth=0.5, marker='*',\
				label='Hopkins, hard X-ray', color='orange')
	except IOError:
		pass

	try:
		df = pd.read_csv('./hopkins_z' + z + '_optical.qlf', names=['Mi', 'nd', 'err_h', 'err_l'], comment='#',\
				delim_whitespace=True, dtype=float, skiprows=5)
		plt.errorbar((90 - df['Mi']) / 2.5 - Lsun, df['nd'], yerr=[df['err_h'], df['err_l']], capsize=2, linewidth=0, elinewidth=0.5, marker='d',\
				label='Hopkins, optical', color='magenta')
	except IOError:
		pass

	try:
		df = pd.read_csv('./hopkins_z' + z + '_IR.qlf', names=['Mi', 'nd', 'err_h', 'err_l'], comment='#',\
				delim_whitespace=True, dtype=float, skiprows=5)
		plt.errorbar((90 - df['Mi']) / 2.5 - Lsun, df['nd'], yerr=[df['err_h'], df['err_l']], capsize=2, linewidth=0, elinewidth=0.5, marker='s',\
				label='Hopkins, infrared', color='k')
	except IOError:
		pass


	df = pd.read_csv('../Ueda_QLF/ueda_z' + z + '.qlf', names=['Mi', 'nd', 'err_h', 'err_l'], comment='#',\
				delim_whitespace=True, dtype=float, skiprows=5)
	plt.errorbar((90 - df['Mi']) / 2.5 - Lsun, df['nd'], yerr=[df['err_h'], df['err_l']], capsize=2, linewidth=0, elinewidth=0.5, marker='o',\
			label='Ueda, X-ray', color='r')
	plt.xlabel(r'$\log\ L_{\mathrm{bol}} [\mathrm{ergs\cdot s^{-1}}]$')
	plt.ylabel(r'$\log\ \mathrm{Number\ Density}$')
	plt.ylim(-9.5, -3.3)
	plt.xlim((10.1, 15))
	plt.legend()
	plt.show()
	plt.close()
