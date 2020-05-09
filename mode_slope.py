#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Determine Sigma from pi times the slope of the mode number
# as a function of lambda, with a multi-measurement-elimination jackknife

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than auto-correlation time)
# We ignore any partial blocks at the end
# Third and fourth are limits of the fit
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <la_min> <la_max>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
la_min = float(sys.argv[3])
la_max = float(sys.argv[4])
runtime = -time.time()

# Make sure results directory exists
if not os.path.isdir('results'):
  print "ERROR: results/ does not exist"
  sys.exit(1)

# Extract lattice volume from path
# For now assume L and Nt are both two-digit numbers
path = os.getcwd()
path = path.replace('mnt', '')    # Accommodate Barkla filesystem
temp = path.split('nt')
L = int(temp[0][-2:])    # Last two digits before 'nt'
Nt = int(temp[1][:2])    # First two digits after 'nt'
vol = L**3 * Nt
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Accumulate eigenvalues from output files
# Construct jackknife samples by skipping over a block of files
# These files may include various numbers of (squared) eigenvalues
# We only consider files with at least 100, and ignore all past 100
# Finally, we use the minimum 100th eigenvalue as a maximum cutoff
files = 'eig.*'
allFiles = glob.glob(files)
cfgs = []
for filename in allFiles:
  cfg = int((filename.split('.'))[-1])  # Number after last .
  if cfg >= cut:
    cfgs.append(cfg)
cfgs.sort()

nblocks = int(np.floor((cfgs[-1] - cfgs[0]) / float(block_size)))
jkslopes = np.empty(nblocks, dtype = np.float)
j = 0

begin = cfgs[0]   # Where each block begins, to be incremented
while (begin + block_size) < cfgs[-1]:
  count = 0;      eig = [];     max_eig = []
  for i in cfgs:            # Make jackknife sample
    if i >= begin and i < (begin + block_size):
      continue
    count += 1
    toOpen = 'eig.' + str(i)
    for line in open(toOpen):
      if line.startswith('EIG'):
        temp = line.split()
        la = np.sqrt(float(temp[2]))
        if la_min <= la and la <= la_max:
          eig.append(la);

  # If we don't have anything to print, then we're done
  if len(eig) == 0:
    print "ERROR: no data in jackknife sample %d", i
    sys.exit(1)

  # Extract slope from linear fit
  # Normalizing nu by the volume and number of eigenvalue measurements
  eig.sort()
  x = np.array(eig)
  norm = 4 * count * vol    # Normalize for single continuum flavor
  y = np.arange(len(eig)) / float(norm)
  out = np.polyfit(x, y, 1)
  jkslopes[j] = out[0]
  begin += block_size
  j += 1

# Check that we went through the expected number of samples
if j != nblocks:
  print "ERROR: no data in jackknife sample %d" % j
  sys.exit(1)

# Now we can average over jackknife samples and print out results
ave = np.average(jkslopes)
var = (nblocks - 1.) * np.sum((jkslopes - ave)**2) / float(nblocks)
outfilename = 'results/Sigma.dat'
outfile = open(outfilename, 'w')
print >> outfile, "%.6g %.4g # %d" % (np.pi * ave, np.pi * np.sqrt(var), nblocks)
outfile.close()

runtime += time.time()
print "# Runtime: %.2g seconds" % runtime
# ------------------------------------------------------------------
