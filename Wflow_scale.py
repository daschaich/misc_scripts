#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# This script determines the Wilson flow scale t0,
# defined as the t for which t^2 E(t) = target (traditionally 0.3)

# First make sure directory for output file exists
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than auto-correlation time)
# We discard any partial blocks at the end
# Third is target <t^2E> for which we find the corresponding sqrt(8t0)
# Fourth is base name of files to analyze, including directory
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <target> <dir/tag>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
target = float(sys.argv[3])
tag = str(sys.argv[4])
runtime = -time.time()

# Construct and sort list of output files
cfgs = []
temp = tag + '*'
for filename in glob.glob(temp):
  cfgs.append(int((filename.split('.'))[-1]))  # Number after last .
cfgs.sort()

# Check to make sure the arguments are appropriate
if len(cfgs) == 0:
  print "ERROR: no files named", temp
  sys.exit(1)

# Extract lattice volume from first output file
firstFile = tag + str(cfgs[0])
for line in open(firstFile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  elif line.startswith('iseed '):
    break   # Done scanning through file
if L > Nt:
  L = Nt    # Take minimum
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Find t0 in each file -- might as well save to data/ directory
outfilename = 'data/scale.csv'
outfile = open(outfilename, 'w')
print >> outfile, "MDTU,sqrt(8t0)"
count = 0
missing = 0
for i in cfgs:
  toOpen = tag + str(i)
  for line in open(toOpen):
    if line.startswith('WFLOW '):
      temp = line.split()
      tSqE = np.fabs(float(temp[4]))    # For 64nt128...
      if tSqE > target:
        count += 1
        print >> outfile, "%d,%.8g" % (i, np.sqrt(8.0 * float(temp[1])))
        break             # Don't bother to check file for completion
    elif line.startswith('RUNNING COMPLETED'):
      if i > cut:
        print "WARNING: Measurement %d never reached target %.2g" % (i, target)
      missing += 1
outfile.close()

if not count == len(cfgs) - missing:
  print "ERROR: Have %d files but only %d measurements" % (len(cfgs), count)
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now analyze newly-created data/t0 file
count = 0
ave = 0.0         # Accumulate within each block
datList = []
begin = cut       # Where each block begins, to be incremented
for line in open(outfilename):
  if line.startswith('MDTU'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU < cut:
    continue
  elif MDTU >= begin and float(temp[0]) < (begin + block_size):
    ave += float(temp[1])
    count += 1
  elif MDTU >= (begin + block_size):  # Move on to next block
    datList.append(ave / float(count))
    begin += block_size
    count = 1                     # Next block begins with this line
    ave = float(temp[1])

# Check special case that the final file fills the last block
if MDTU >= (begin + block_size - cfgs[-1] + cfgs[-2]):
  datList.append(ave / float(count))
  begin += block_size

# Now print mean and standard error, requiring N>1
if len(datList) < 2:
  print "ERROR: Need multiple blocks to average"
  sys.exit(1)

dat = np.array(datList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
print "%.2g %.8g %.4g # %d" % (target, ave, err, N)

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
