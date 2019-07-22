#!/usr/bin/python
import os
import sys
import glob
import numpy as np
from scipy.misc import comb   # For N choose k
# ------------------------------------------------------------------
# Parse dygraph time-series data file to compute cumulants ka_n
# for the (Wilson-flowed) Polyakov loop
# (Other targets may be added in the future)

# Following arXiv:1411.7461, this should give us the susceptibility ka_2/V,
# the skewness ka_3/ka_2^{3/2} and the kurtosis ka_4/ka_2
# We can check the first against direct computations from ts_suscept.py
# (accounting for offline normalization by volume)

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<cut> <block>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
Norder = 4                          # How many cumulants to consider
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# Check that we actually have data to average
MDTUfile = 'data/TU.csv'
good = -1
for line in open(MDTUfile):
  if line.startswith('t'):
    continue
  temp = line.split(',')
  if float(temp[1]) > cut:
    good = 1
    break

if good == -1:
  print "Error: no data to analyze",
  print "since cut=%d but we only have %d MDTU" % (cut, float(temp[1]))
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For Wpoly(_mod), grab the last (fourth) number for c=0.5
# For poly* and pbp, just grab the single number after the MDTU label
for obs in ['Wpoly', 'Wpoly_mod', 'poly_r', 'poly_mod']:
  count = 0
  ave = 0.0         # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= cut:
      continue
    elif MDTU > begin and MDTU < (begin + block_size):
      if obs == 'Wpoly' or obs == 'Wpoly_mod':
        tr = float(temp[-1])
      elif obs == 'poly_mod' or obs == 'poly_r':
        tr = float(temp[1])
      ave += tr
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      if count == 0:
        print "WARNING: no %s data to average at %d MDTU" % (obs, int(MDTU))
        skip = 1
        break
      datList.append(ave / count)
      begin += block_size
      count = 1                         # Next block begins here
      if obs == 'Wpoly' or obs == 'Wpoly_mod':
        tr = float(temp[-1])
      elif obs == 'poly_mod' or obs == 'poly_r':
        tr = float(temp[1])
      ave = tr

  # Require multiple blocks, N>1
  N = len(datList)
  if N < 2:
    print "ERROR: need multiple blocks for jackknifing"
    sys.exit(1)

  # Now construct jackknife samples through single-block elimination
  
  # Jackknife first Norder cumulants ka (and moments mu)
  # through single-block elimination
  dat = np.array(datList)
  muJK = np.empty((Norder, N), dtype = np.float)
  kaJK = np.empty((Norder, N), dtype = np.float)
  for i in range(N):
    JKdat = np.delete(dat, i)             # Remove ith data point
    for n in range(Norder):
      muJK[n][i] = np.mean(JKdat**(n + 1), dtype = np.float64)
      kaJK[n][i] = muJK[n][i]
      for m in range(n):
        # comb(N, k) is N-choose-k
        kaJK[n][i] -= comb(n, m) * kaJK[m][i] * muJK[n - m - 1][i]

  # Average over jackknife samples for cumulants through ka_N
  # Checked ka_1 against average and ka_2 against susceptibility...TODO
  outfilename = 'results/' + obs + '.cumulants'
  outfile = open(outfilename, 'w')
  for n in range(Norder):
    ave = np.average(kaJK[n])
    var = (N - 1.0) * np.sum((kaJK[n] - ave)**2) / float(N)
    print >> outfile, "kappa_%d %.8g %.4g # %d" % (n+1, ave, np.sqrt(var), N)

  # Average over jackknife samples for skewness and kurtosis
  if (Norder > 3):
    # ka_3 / ka_2^{1.5} (indexing from 1)
    skewJK = kaJK[2] / np.power(kaJK[1], 1.5)
    ave = np.average(skewJK)
    var = (N - 1.0) * np.sum((skewJK - ave)**2) / float(N)
    print >> outfile, "skewness %.8g %.4g # %d" % (ave, np.sqrt(var), N)

    # ka_4 / ka_2^2 (indexing from 1)
    kurtJK = kaJK[3] / kaJK[1]**2
    ave = np.average(kurtJK)
    var = (N - 1.0) * np.sum((kurtJK - ave)**2) / float(N)
    print >> outfile, "kurtosis %.8g %.4g # %d" % (ave, np.sqrt(var), N)
  outfile.close()
# ------------------------------------------------------------------
