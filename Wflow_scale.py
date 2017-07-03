#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# This script determines the Wilson flow scales t0 and w0
# t0 is defined as the t for which t^2 E(t) = target (traditionally 0.3)
# w0 is defined similarly for the log-derivative of t^2 E(t)
# For the former we consider the clover discretization
# with and without perturbative finite-volume corrections
# For now we don't consider the plaquette discretization
# Now constructing blocks based on number of measurements, not MDTU
# TODO: Incorporate perturbative corrections into log-derivative
# TODO: Reconstruct log-derivative for plaquette discretization
# TODO: Could do different targets for all observables

# First make sure directory for output file exists
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# Parse arguments: first is initial measurement,
# second is number of measurements per block
# We discard any partial blocks at the end
# Third is target <t^2E> for which we find the corresponding sqrt(8t0)
# Fourth is base name of files to analyze, including directory
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<first> <num> <target> <dir/tag>"
  sys.exit(1)
first = int(sys.argv[1])
num = int(sys.argv[2])
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
pertFile = '../pert_corr.dat'
if not os.path.isfile(pertFile):
  print "ERROR:", pertFile, "does not exist"
  sys.exit(1)

# The perturbative corrections require epsilon = 0.01
dt = 0.01

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
# Consider clover tSqE with and without perturbative improvement
# Skip plaquette discretization but do (clover) w0 at the same time
# Include bins with a few missing measurements, for consistent analyses
outfilename = 'data/scale.csv'
outfile = open(outfilename, 'w')
print >> outfile, "# t^2 E = %g" % target
print >> outfile, "# meas,sqrt(8t0)_raw,sqrt(8t0)_pert,w0_raw"
t_pert, pert = np.loadtxt(pertFile, unpack=True)
count = 0
missing = 0
for i in cfgs:
  toOpen = tag + str(i)
  done = [0, 0, 0]        # Stop when all three are done
  t0_raw = -1.0
  t0_pert = -1.0
  w0_raw = -1.0
  index = 0
  for line in open(toOpen):
    if done == [1, 1, 1]:
      count += 1
      print >> outfile, "%d,%.8g,%.8g,%.8g" % (i, t0_raw, t0_pert, w0_raw)
      break               # Don't bother to check file for completion
    if line.startswith('epsilon'):    # Check epsilon
      if not dt == float((line.split())[1]):
        print "ERROR:", toOpen, "uses wrong epsilon",
        print (line.split())[1]
        sys.exit(1)
    # Format: WFLOW t plaq E t^2*E t*d(t^2*E)/dt 12t^2*(3-plaq) top.charge
    elif line.startswith('WFLOW '):
      temp = line.split()
      t = float(temp[1])
      tSqEraw = np.fabs(float(temp[4]))    # fabs for 64nt128...
      deriv_raw = np.fabs(float(temp[5]))
      # Extract and check perturbative correction for this t
      if t_pert[index] != t:
        print "ERROR: t=%.2g in %s doesn't match t=%.2g in  %s" \
              % (t_pert[index], pertFile, t, toOpen)
        sys.exit(1)
      tSqEpert = tSqEraw / pert[index]
      index += 1
      # Check if any of the three have just crossed the target
      if done[0] == 0 and tSqEraw > target:
        t0_raw = np.sqrt(8.0 * t)
        done[0] = 1
      if done[1] == 0 and tSqEpert > target:
        t0_pert = np.sqrt(8.0 * t)
        done[1] = 1
      if done[2] == 0 and deriv_raw > target:
        w0_raw = np.sqrt(t)
        done[2] = 1
    elif line.startswith('RUNNING COMPLETED'):
      if i > first:
        print "WARNING: Measurement %d never reached target %.2g" % (i, target)
      # TODO: This inserts extraneous spaces
      # May need to import print() from future to fix
      print >> outfile, "%d," % i,
      if done[0] == 0:
        print >> outfile, "null,",
      else:
        print >> outfile, "%.8g," % t0_raw,
      if done[1] == 0:
        print >> outfile, "null,",
      else:
        print >> outfile, "%.8g," % t0_pert,
      if done[2] == 0:
        print >> outfile, "null"
      else:
        print >> outfile, "%.8g" % w0_raw
      missing += 1
outfile.close()

if not count == len(cfgs) - missing:
  print "ERROR: Have %d files but only %d measurements" % (len(cfgs), count)
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now analyze newly-created data/t0 file
# Print each jackknife bin for correlated fitting
binfilename = 'data/scale_bins.csv'
binfile = open(binfilename, 'w')
print >> binfile, "# t^2 E = %g" % target
print >> binfile, "# bin,start,stop,sqrt(8t0)_raw,sqrt(8t0)_pert,w0_raw"
missing = 0
count = 0
ave_t0raw = 0.0         # Accumulate within each block
ave_t0pert = 0.0
ave_w0raw = 0.0
t0rawList = []
t0pertList = []
w0rawList = []
begin = first     # Where each block begins, to be incremented
for line in open(outfilename):
  if line.startswith('#'):
    continue
  temp = line.split(',')
  i = float(temp[0])
  if i < first:
    continue
  elif count < num:
    if 'null' in line:    # TODO: Could refine to handle different obs
      missing += 1
    else:
      ave_t0raw += float(temp[1])
      ave_t0pert += float(temp[2])
      ave_w0raw += float(temp[3])
    count += 1
  elif count == num:                        # Move on to next block
    if count == missing:
      print "ERROR: All measurements in this block don't reach target"
      print "       Need to do something smarter..."
      sys.exit(1)
    t0rawList.append(ave_t0raw / float(count - missing))
    t0pertList.append(ave_t0pert / float(count - missing))
    w0rawList.append(ave_w0raw / float(count - missing))
    print >> binfile, "%d,%d,%d,%.8g,%.8g,%.8g" \
                      % (len(t0rawList), begin, previous, \
                         t0rawList[-1], t0pertList[-1], w0rawList[-1])
    begin = i
    missing = 0
    count = 1                     # Next block begins with this line
    if 'null' in line:
      ave_t0raw = 0.0
      ave_t0pert = 0.0
      ave_w0raw = 0.0
      missing = 1
    else:
      ave_t0raw = float(temp[1])
      ave_t0pert = float(temp[2])
      ave_w0raw = float(temp[3])
  else: # count[0] > num:                       # Sanity check
    print "ERROR: Something funny is going on..."
    sys.exit(1)
  previous = i

# Check special case that the final file fills the last block
if count == num:
  t0rawList.append(ave_t0raw / float(count - missing))
  t0pertList.append(ave_t0pert / float(count - missing))
  w0rawList.append(ave_w0raw / float(count - missing))
  print >> binfile, "%d,%d,%d,%.8g,%.8g,%.8g" \
                    % (len(t0rawList), begin, previous, \
                       t0rawList[-1], t0pertList[-1], w0rawList[-1])
binfile.close()

# Now print mean and standard error, requiring N>1
if len(t0rawList) < 2:
  print "ERROR: Need multiple blocks to average"
  sys.exit(1)

print "%.2g" % target,

dat = np.array(t0rawList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
print "%.8g %.4g" % (ave, err),

dat = np.array(t0pertList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
print "%.8g %.4g" % (ave, err),

dat = np.array(w0rawList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
print "%.8g %.4g # %d" % (ave, err, N)

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
