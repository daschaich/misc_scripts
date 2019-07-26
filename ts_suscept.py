#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph time-series data file to compute susceptibility,
# skewness and (excess) kurtosis
# for the plaquette, (Wilson-flowed) Polyakov loop and chiral condensate
# (Other targets may be added in the future)

# When there are multiple measurements in a block,
# we can't use the naive sum_i(dat_i - vev)^n expressions!
# Need to accumulate dat^n within each block instead

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<cut> <block>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
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
# For plaquette, average two data per line
# For Wpoly(_mod), grab the last (fourth) number for c=0.5
# For poly* and pbp, just grab the single number after the MDTU label
for obs in ['plaq', 'Wpoly', 'Wpoly_mod', 'poly_r', 'poly_mod', 'pbp']:
  count = 0
  ave = 0.0         # Accumulate within each block
  aveSq = 0.0
  datList = []
  sqList = []
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
      if obs == 'plaq':
        tr = 0.5 * (float(temp[1]) + float(temp[2]))
      elif obs == 'Wpoly' or obs == 'Wpoly_mod':
        tr = float(temp[-1])
      elif obs == 'poly_mod' or obs == 'pbp' or obs == 'poly_r':
        tr = float(temp[1])
      ave += tr
      aveSq += tr * tr
      count += 1

    # Done with this block
    elif MDTU == (begin + block_size):
      if obs == 'plaq':
        tr = 0.5 * (float(temp[1]) + float(temp[2]))
      elif obs == 'Wpoly' or obs == 'Wpoly_mod':
        tr = float(temp[-1])
      elif obs == 'poly_mod' or obs == 'pbp' or obs == 'poly_r':
        tr = float(temp[1])
      ave += tr
      aveSq += tr * tr
      count += 1

      # Record this block
      datList.append(ave / count)
      sqList.append(aveSq / count)
      begin += block_size

      # Re-initialize for next block
      count = 0
      ave = 0.0
      aveSq = 0.0

    # This should never happen
    elif MDTU > (begin + block_size):
      print "ERROR: Unexpected behavior in %s, aborting" % obsfile
      sys.exit(1)

  # Require multiple blocks, N>1
  N = len(datList)
  if N < 2:
    print "ERROR: need multiple blocks to take average"
    sys.exit(1)

  # Now construct jackknife samples through single-block elimination
  # Do normalization offline, so here just have
  #   chi = <obs^2> - <obs>^2
  dat = np.array(datList, dtype = np.float64)
  sq = np.array(sqList, dtype = np.float64)
  chi = np.zeros(N, dtype = np.float64)
  for i in range(N):    # Jackknife samples
    vev = np.mean(np.delete(dat, i))
    sq_vev = np.mean(np.delete(sq, i))
    chi[i] = sq_vev - vev * vev

  # Sanity check -- compare against averages computed separately
#  print obs, "ave = %.8g" % np.mean(dat, dtype = np.float64)

  # Now we can average over jackknife samples and print out results
  ave = np.mean(chi)
  var = (N - 1.0) * np.mean((chi - ave)**2)
  outfilename = 'results/' + obs + '.suscept'
  outfile = open(outfilename, 'w')
  print >> outfile, "%.8g %.4g # %d" % (ave, np.sqrt(var), N)
  outfile.close()
# ------------------------------------------------------------------
