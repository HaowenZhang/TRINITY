import os
import sys
import subprocess

# This is the python wrapper of the C binary file, TRINITY/bin/gen_edd_r_color,
# which generates various galaxy&SMBH properties, as functions of redshift and
# halo mass

if len(sys.argv) < 2:
	print('Usage: %s mass_cache (mcmc output).' % sys.argv[0])
	exit(1)

dir_bin = "../../bin/"
name_bin = dir_bin + sys.argv[0].split('.')[0]
print(name_bin)

cps = subprocess.run([name_bin, sys.argv[1]])

print(cps.stdout, file=sys.stdout)
print(cps.stderr, file=sys.stderr)
