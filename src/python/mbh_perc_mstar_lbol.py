import os
import sys
import subprocess

# This is the python wrapper of the C binary file, TRINITY/bin/mbh_perc_mstar_lbol,
# which generates 16th, 50th, and 84th percentiles in SMBH mass at an input redshift,
# on a grid of galaxy mass, assuming an input random scatter in observed SMBH mass
# around the intrinsic mass.

if len(sys.argv) < 5:
	print('Usage: %s z lbol_lim sigma_mbh mass_cache (mcmc output).' % sys.argv[0])
	exit(1)

dir_bin = "../../bin/"
name_bin = dir_bin + sys.argv[0].split('.')[0]
print(name_bin)

cps = subprocess.run([name_bin, sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]])

print(cps.stdout, file=sys.stdout)
print(cps.stderr, file=sys.stderr)
