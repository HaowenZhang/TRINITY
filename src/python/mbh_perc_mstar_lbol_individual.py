import os
import sys
import subprocess

# This is the python wrapper of the C binary file, TRINITY/bin/mbh_perc_mstar_lbol_individual,
# which calculates how many sigma the SMBH masses (given by the quasar_catalog) are away from
# the SMBH mass--galaxy mass relations for quasars at their luminosities (also given by the catalog).

if len(sys.argv) < 3:
	print('Usage: %s mass_cache quasar_catalog (mcmc output).' % sys.argv[0])
	exit(1)

dir_bin = "../../bin/"
name_bin = dir_bin + sys.argv[0].split('.')[0]
print(name_bin)

cps = subprocess.run([name_bin, sys.argv[1], sys.argv[2]])

print(cps.stdout, file=sys.stdout)
print(cps.stderr, file=sys.stderr)
