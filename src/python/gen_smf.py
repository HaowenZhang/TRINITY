import os
import sys
import subprocess

# This is the python wrapper of the C binary file, TRINITY/bin/gen_smf,
# which generates galaxy stellar mass functions in the redshift bin
# (z_low, z_high).

if len(sys.argv) < 4:
	print('Usage: %s z_low z_high mass_cache (mcmc output).' % sys.argv[0])
	exit(1)

dir_bin = "../../bin/"
name_bin = dir_bin + sys.argv[0].split('.')[0]
print(name_bin)

cps = subprocess.run([name_bin, sys.argv[1], sys.argv[2], sys.argv[3]])

print(cps.stdout, file=sys.stdout)
print(cps.stderr, file=sys.stderr)
