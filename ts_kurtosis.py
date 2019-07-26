#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph time-series data file to compute (excess) kurtosis
# for the (Wilson-flowed) Polyakov loop
# (Other targets may be added in the future)

# Also check susceptibility and skewness
# TODO: The former can be compared with the result from ts_suscept.py...

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
TODO: I THINK I NEED TO BLOCK ALL OF tr^2, tr^3 and tr^4!!!...
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

    # Accumulate within block
    elif MDTU > begin and MDTU < (begin + block_size):
      if obs == 'Wpoly' or obs == 'Wpoly_mod':
        tr = float(temp[-1])
      elif obs == 'poly_mod' or obs == 'poly_r':
        tr = float(temp[1])
      ave += tr
      count += 1

    # Done with this block
    elif MDTU == (begin + block_size):    # Done with this block
      if obs == 'Wpoly' or obs == 'Wpoly_mod':
        tr = float(temp[-1])
      elif obs == 'poly_mod' or obs == 'poly_r':
        tr = float(temp[1])
      ave += tr
      count += 1

      # Record this block
      datList.append(ave / count)
      begin += block_size

      # Re-initialize for next block
      count = 0
      ave = 0.0

    # This should never happen
    elif MDTU > (begin + block_size):
      print "ERROR: Unexpected behavior in %s, aborting" % obsfile
      sys.exit(1)

  # Require multiple blocks, N>1
  N = len(datList)
  if N < 2:
    print "ERROR: need multiple blocks for jackknifing"
    sys.exit(1)

  # Now construct jackknife samples through single-block elimination
  # Check susceptibility and skewness in addition to kurtosis
  #   chi = (1/N) sum_i (dat_i - vev)^2
  #   S = [(1/N) sum_i (dat_i - vev)^3] / chi^{3/2}
  #   ka = [(1/N) sum_i (dat_i - vev)^4] / chi^2
  chi = np.zeros(N, dtype = np.float64)
  S = np.zeros(N, dtype = np.float64)
  ka = np.zeros(N, dtype = np.float64)
  dat = np.array(datList, dtype = np.float64)
  for i in range(N):
    JKdat = np.delete(dat, i)             # Remove ith data point
    vev = np.mean(JKdat)
    chi[i] = np.mean((JKdat - vev)**2)
    S[i] = np.mean((JKdat - vev)**3) / (np.power(chi[i], 1.5))
    ka[i] = np.mean((JKdat - vev)**4) / (chi[i]**2) - 3.0

  # Average over jackknife samples and print
  outfilename = 'results/' + obs + '.kurtosis'
  outfile = open(outfilename, 'w')

  ave = np.mean(chi)
  var = (N - 1.0) * np.mean((chi - ave)**2)
  print >> outfile, "suscept %.8g %.4g # %d" % (ave, np.sqrt(var), N)

  ave = np.mean(S)
  var = (N - 1.0) * np.mean((S - ave)**2)
  print >> outfile, "skewness %.8g %.4g # %d" % (ave, np.sqrt(var), N)

  ave = np.mean(ka)
  var = (N - 1.0) * np.mean((ka - ave)**2)
  print >> outfile, "kurtosis %.8g %.4g # %d" % (ave, np.sqrt(var), N)
  outfile.close()
# ------------------------------------------------------------------
