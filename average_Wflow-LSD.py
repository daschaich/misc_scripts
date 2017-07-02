#!/usr/bin/python
import os
import sys
import time
import glob
import numpy as np
# ------------------------------------------------------------------
# Construct averages and standard errors of the gradient flow coupling
# Process all t up to given target <= L^2 / 32
# Now constructing blocks based on number of measurements, not MDTU

# Parse arguments: first is initial measurement,
# second is number of measurements per block
# We discard any partial blocks at the end
# Third is base name of files to analyze, including directory
# Fourth argument is t-shift improvement parameter tau
# Fifth argument is maximum t to consider, at most L^2 / 32
# Optional sixth argument tells us which observable to consider
if len(sys.argv) < 6:
  print "Usage:", str(sys.argv[0]),
  print "<first> <num> <dir/tag> <tau> <target> <obs>"
  sys.exit(1)
first = int(sys.argv[1])
num = int(sys.argv[2])
tag = str(sys.argv[3])
tau = float(sys.argv[4])
target = float(sys.argv[5])
runtime = -time.time()

# Choose which observable to use -- require 'plaq' as specific argument
# Will probably have tau=0, so don't include in output file name
# TODO: Set up plaquette log-derivative
plaq = -1
deriv = -1
outfilename = 'results/Wflow_coupling.dat'
if len(sys.argv) > 6:
  if str(sys.argv[6]) == 'plaq':
    plaq = 1
    outfilename = 'results/Wplaq_coupling.dat'
  elif str(sys.argv[6]) == 'deriv':
    deriv = 1
    outfilename = 'results/Wflow_deriv.dat'
    print "Warning: Not yet sure about factors of (t + tau) in derivative"
  else:
    print "Warning: Implicitly using clover observable"
# ------------------------------------------------------------------



# ------------------------------------------------------------------
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
# Generally will have Nt = 2L, but let's not assume this
dt = 0.01
lookup = []
firstFile = tag + str(cfgs[0])
if not os.path.isfile(firstFile):
  print "ERROR:", firstFile, "does not exist"
  sys.exit(1)
for line in open(firstFile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
    # Set L to minimum of nx and nt
    if Nt < L:
      L = Nt
    if L < 0 or Nt < 0:
      print "ERROR: couldn't extract lattice size from", firstFile
      sys.exit(1)

  # Count number of points of t >= tau
  # Perturbative correction would bring in a factor of 1/(t - tau)^2
  # so we need non-zero (t - tau)...
  # We also need to stop at target even if the file has more data
  elif line.startswith('WFLOW '):
    temp = line.split()
    t = float(temp[1])
    if t > target:              # Target considers unshifted t
      break
    t -= tau                    # Include t-shift
    if t < 0.001:               # Start when t > tau
      continue
    lookup.append(t)
Npt = len(lookup)

if target > float(L**2) / 32.0:         # Sanity check
  print "ERROR: can't reach target %.2g > %.2g" % (target, float(L**2) / 32.0)
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct jackknife blocks from Wflow files
ave = np.zeros(Npt, dtype = np.float)       # Accumulate within each block
count = np.zeros(Npt, dtype = np.int)
index = 0
datList = [[] for x in range(Npt)]
begin = first     # Where each block begins, to be incremented
for i in cfgs:
  if i < first:
    continue
  elif count[0] < num:
    WflowFile = tag + str(i)
    if not os.path.isfile(WflowFile):
      print "WARNING:", WflowFile, "does not exist"
      continue
    index = 0
    for line in open(WflowFile):
      if line.startswith('epsilon'):    # Check epsilon
        if not dt == float((line.split())[1]):
          print "ERROR:", WflowFile, "uses wrong epsilon",
          print (line.split())[1]
          sys.exit(1)
      elif line.startswith('tmax '):    # Check tmax
        if float((line.split())[1]) < target:
          print "ERROR:", WflowFile, "uses too-small tmax",
          print (line.split())[1]
          sys.exit(1)
      # Format: WFLOW t plaq E t^2*E t*d(t^2*E)/dt 12t^2*(3-plaq) top.charge
      elif line.startswith('WFLOW '):
        temp = line.split()
        t = float(temp[1])
        if t > target:
          break
        t -= tau                                # Include t-shift
        if t < 0.001:                           # Start when t > tau
          continue
        if plaq < 0 and deriv < 0:              # t**2 * E(t + tau)
          ave[index] += t**2 * np.fabs(float(temp[3]))   # fabs for 64nt128...
        elif plaq > 0:
          # Data is (t + tau)**2 * E(t + tau)...
          rescale = t / float(temp[1])
          ave[index] += rescale**2 * float(temp[6])
        elif deriv > 0:
          # Data is (t + tau) * d[(t + tau)**2 E(t + tau)]/d[t + tau]...
          rescale = t / float(temp[1])
          if float(temp[3]) < 0:
            rescale *= -1.0                     # For 64nt128...
          ave[index] += rescale * float(temp[5])
        count[index] += 1
        index += 1
    # Check that correct number of points were read
    if not index == Npt:
      print "ERROR: Only", index, "of", Npt, "points counted"
      sys.exit(1)

  elif count[0] == num:                         # Move on to next block
    for index in range(Npt):
      if count[index] >= 1:
        datList[index].append(ave[index] / float(count[index]))
      if count[index] != count[0]:
        print "WARNING: mismatch %d vs %d for index %d in block %d" \
          % (count[index], count[0], index, len(datList[0]))
    ave = np.zeros(Npt, dtype = np.float)
    count = np.zeros(Npt, dtype = np.int)
    begin = i

    WflowFile = tag + str(i)
    if not os.path.isfile(WflowFile):
      print "WARNING:", WflowFile, "does not exist"
      continue
    index = 0
    for line in open(WflowFile):
      if line.startswith('epsilon'):    # Check epsilon
        if not dt == float((line.split())[1]):
          print "ERROR:", WflowFile, "uses wrong epsilon",
          print (line.split())[1]
          sys.exit(1)
      elif line.startswith('tmax '):    # Check tmax
        if float((line.split())[1]) < target:
          print "ERROR:", WflowFile, "uses too-small tmax",
          print (line.split())[1]
          sys.exit(1)
      elif line.startswith('WFLOW '):
        temp = line.split()
        t = float(temp[1])
        if t > target:
          break
        t -= tau                                # Include t-shift
        if t < 0.001:                           # Start when t > tau
          continue
        if plaq < 0:
          ave[index] += t**2 * float(temp[3])   # t**2 * E(t + tau)
        elif plaq > 0:
          # Data is (t + tau)**2 * E(t + tau)...
          rescale = t / float(temp[1])
          ave[index] += rescale**2 * float(temp[6])
        elif deriv > 0:
          # Data is (t + tau) * d[(t + tau)**2 E(t + tau)]/d[t + tau]...
          rescale = t / float(temp[1])
          ave[index] += rescale * float(temp[5])
        count[index] += 1
        index += 1
    # Check that correct number of points were read
    if not index == Npt:
      print "ERROR: Only", index, "of", Npt, "points counted"
      sys.exit(1)
  else: # count[0] > num:                       # Sanity check
    print "ERROR: Something funny is going on..."
    sys.exit(1)

# The file_tag check in the loop above means that
# full blocks aren't recorded until we reach the next file
# Check special case that the final file fills the last block
if count[0] == num:
  for index in range(Npt):
    if count[index] >= 1:
      datList[index].append(ave[index] / float(count[index]))
    if count[index] != count[0]:
      print "WARNING: mismatch %d vs %d for index %d in block %d" \
            % (count[index], count[0], index, len(datList[0]))
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
N = float(len(datList[0]))
for index in range(Npt):
  if len(datList[index]) != len(datList[0]):
    print "ERROR: mismatch %d vs %d for index %d" \
          % (len(datList[index]), len(datList[0]), index)
    sys.exit(1)
if N <= 1:
  print "ERROR: Not enough blocks (%d) to average" % N
  sys.exit(1)

outfile = open(outfilename, 'w')
print >> outfile, "# tau=%g with %d blocks" % (tau, int(N))

# Separately compute perturbative finite-volume + zero-mode corrections
# and save to lookup table that can be `paste`d and `awk`d when plotting
if deriv > 0:
  gprop = 1.0
else:
  gprop = 128.0 * 3.14159**2 / 24.0
for i in range(Npt):
  dat = np.array(datList[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1)
  print >> outfile, "%.4g %.8g %.4g" % (lookup[i], gprop * ave, gprop * err)

runtime += time.time()
print >> outfile, "# Runtime: %d seconds" % int(runtime)
outfile.close()
# ------------------------------------------------------------------
