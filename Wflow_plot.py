#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# Plot the Wilson flow t^2<E> as a function of flow time

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than auto-correlation time)
# Third is base name of files to analyze, including directory
# We discard any partial blocks at the end
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <dir/tag>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
tag = str(sys.argv[3])
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

# Figure out maximum safe flow time t=(L/2)^2/8from lattice volume
# Extract lattice volume from first output file
# Also record epsilon, to be checked in each file
firstFile = tag + str(cfgs[0])
for line in open(firstFile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  elif line.startswith('epsilon '):
    epsilon = float((line.split())[1])
  elif line.startswith('alpha_hyp0'):
    break   # Done scanning through file
if L > Nt:
  L = Nt    # Take minimum
tmax = L * L / 32.0
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Go through output files
# Sanity check: monitor t to make sure we combine the right data
# Accumulate tSqE and its square to be averaged and stdev'd
# Discard any partial bin at the end
Nmeas = 0
binsize = 0
t = []
dat = []
tSqE = []
tSqE_err = []
for i in cfgs:
  toOpen = tag + str(i)
  with open(toOpen) as infile:
    iter = 0
    for line in infile:
      # Make sure all files use the same epsilon
      if line.startswith('epsilon '):
        temp = line.split()   # Convert line into list
        if float(temp[1]) != epsilon:
          print('ERROR: EPSILON MISMATCH:'),
          print('{0:.2g} vs. {1:.2g}'.format(float(temp[1]), epsilon)),
          print('in file {0:d}'.format(i))
          sys.exit(1)

      # Go.  Initialize list from first file
      # Maintain running averages afterward
      # Sanity checks: compare flow time in each file
      #                ignore flow times that look too large
      elif line.startswith('WFLOW '):
        temp = line.split()   # Convert line into list
        toCheck = float(temp[1])
        if toCheck > tmax:    # Ignore flow times that look too large
          continue
        if Nmeas == 0 and binsize == 0:
          t.append(toCheck)
          dat.append(float(temp[4]))
        else:
          dat[iter] += np.fabs(float(temp[4]))    # For 64nt128...

          # Sanity check
          if toCheck != t[iter]:
            print('ERROR: FLOW TIME MISMATCH:'),
            print('{0:.2g} vs. {1:.2g}'.format(toCheck, t[iter])),
            print('in file {0:d}'.format(i))
            sys.exit(1)
        iter += 1

  # Accumulate results every bin
  binsize += 1
  if binsize == block_size:
    for iter in range(0, len(t)):
      temp = dat[iter] / binsize
      if Nmeas == 0:
        tSqE.append(temp)
        tSqE_err.append(temp * temp)
      else:
        tSqE[iter] += temp
        tSqE_err[iter] += temp * temp
      dat[iter] = 0

    # Prepare for next bin
    binsize = 0
    Nmeas += 1
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print averages and standard deviations
print('# Thermalization cut at configuration {0:d} of {1:d}'
      .format(cut, cfgs[-1]))
print('# {0:d} files --> {1:d} measurements'.format(len(cfgs), Nmeas))
for iter in range(0, len(t)):
  tSqE[iter] /= Nmeas
  tSqE_err[iter] /= Nmeas
  tSqE_err[iter] -= tSqE[iter] * tSqE[iter]
  tSqE_err[iter] = np.sqrt(tSqE_err[iter] / (Nmeas - 1))
  print('{0:.5g} {1:.5g} {2:.5g}'.format(t[iter], tSqE[iter], \
                                         tSqE_err[iter]))

runtime += time.time()
print "# Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
