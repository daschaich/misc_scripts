#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse given file construct averages and standard errors
# with given thermalization cut and block size

# Assume one ensemble per directory (overwrite results files)

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than auto-correlation time)
# We discard any partial blocks at the end
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <file>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
infile = str(sys.argv[3])

MDTU = 0.0
count = 0
ave = 0.0        # Accumulate within each block
datList = []
begin = cut       # Where each block begins, to be incremented
for line in open(infile):
  if line.startswith('M'):
    continue
#  temp = line.split()
  MDTU += 1.0
  if MDTU < cut:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
    ave += float(line)
    count += 1
  elif MDTU >= (begin + block_size):  # Move on to next bloc
    datList.append(ave / float(count))
    begin += block_size
    count = 1                     # Next block begins with this line
    ave = float(line)

# Now print mean and standard error, assuming N>1
dat = np.array(datList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.)
print "%.8g %.4g # %d" % (ave, err, N)
# ------------------------------------------------------------------
