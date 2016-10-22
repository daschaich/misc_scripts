#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Determine V \chi_t = <Q^2> - <Q>^2 from the given data file,
# running a block-elimination jackknife

# Parse arguments: first is file to analyze,
# second is thermalization cut,
# third is block size (should be larger than auto-correlation time)
# We discard any partial blocks at the end
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<file> <cut> <block>"
  sys.exit(1)
filename = str(sys.argv[1])
cut = int(sys.argv[2])
block_size = int(sys.argv[3])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)

# Quickly extract last MDTU to figure out how many measurements we have
for line in open(filename):
  temp = line.split()
  maxTU = float(temp[0])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct array of blocked measurements
dat = []
datSq = []

# Monitor block lengths, starting and ending MDTU
block_data = [[], [], []]
count = 0
begin = cut   # Where each block begins, to be incremented
toAve = 0.
toAveSq = 0.

# File format: MDTU   c=0.2   c=0.3   c=0.4   c=0.5
# We're interested in c=0.5
for line in open(filename):
  temp = line.split()
  MDTU = float(temp[0])
  if MDTU < cut:
    continue
  # If we're done with this block, record it and reset for the next
  if MDTU >= (begin + block_size):
    dat.append(toAve / float(count))
    datSq.append(toAveSq / float(count))

    # Record and reset block data
    block_data[0].append(count)
    count = 0
    block_data[1].append(begin)
    begin += block_size
    block_data[2].append(begin)
    toAve = 0.
    toAveSq = 0.

  # Running averages
  tr = float(temp[4])
  toAve += tr
  toAveSq += tr**2
  count += 1
Nblocks = len(dat)

# If only one block, we should fit the histogram to a gaussian
# But for now we don't have that in this script
if Nblocks == 1:
  Vchi = datSq[0] - dat[0]**2
  print "%.6g from one block" % Vchi
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination
tot = sum(dat)
totSq = sum(datSq)
jkVchi = np.empty(Nblocks, dtype = np.float)
for i in range(Nblocks):
  Q = (tot - dat[i]) / (Nblocks - 1.)
  QSq = (totSq - datSq[i]) / (Nblocks - 1.)
  jkVchi[i] = QSq - Q * Q

# Now we can average over jackknife samples and print out results
ave = np.average(jkVchi)
var = (Nblocks - 1.) * np.sum((jkVchi - ave)**2) / float(Nblocks)
print "%.6g %.4g # %d" % (ave, np.sqrt(var), Nblocks)

# Quick check of blocks
#for i in range(Nblocks):
#  print "# Block %2d has %d measurements from MDTU in [%d, %d)" \
#        % (i + 1, block_data[0][i], block_data[1][i], block_data[2][i])
# ------------------------------------------------------------------

