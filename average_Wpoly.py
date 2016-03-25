#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse ./Wpoly.dat to construct averages and standard errors
# with given thermalization cut and block size

# Assume one ensemble per directory (overwrite results files)

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than auto-correlation time)
# We discard any partial blocks at the end
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<cut> <block>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])

# For now, let's print out both c=0.2 and 0.3,
# the second and third data on each line
count = 0
ave2 = 0.0        # Accumulate within each block
ave3 = 0.0
datList2 = []
datList3 = []
begin = cut       # Where each block begins, to be incremented
flowfile = 'Wpoly.dat'
for line in open(flowfile):
  if line.startswith('M'):
    continue
  temp = line.split()
  MDTU = float(temp[0])
  if MDTU < cut:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
    ave2 += float(temp[1])
    ave3 += float(temp[2])
    count += 1
  elif MDTU >= (begin + block_size):  # Move on to next bloc
    datList2.append(ave2 / float(count))
    datList3.append(ave3 / float(count))
    begin += block_size
    count = 1                     # Next block begins with this line
    ave2 = float(temp[1])
    ave3 = float(temp[2])

# Now print mean and standard error, assuming N>1
dat2 = np.array(datList2)
dat3 = np.array(datList3)
N = np.size(dat2)
ave2 = np.mean(dat2, dtype = np.float64)
err2 = np.std(dat2, dtype = np.float64) / np.sqrt(N - 1.)
ave3 = np.mean(dat3, dtype = np.float64)
err3 = np.std(dat3, dtype = np.float64) / np.sqrt(N - 1.)
outfilename = 'results/Wpoly.dat'
outfile = open(outfilename, 'w')
print >> outfile, "# c=0.2 err c=0.3 err # Nblocks"
print >> outfile, "%.8g %.4g %.8g %.4g # %d" % (ave2, err2, ave3, err3, N)
outfile.close()
# ------------------------------------------------------------------
