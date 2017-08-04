#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# This script determines the Wilson flow scale t0 from QHMC output
# t0 is defined as the t for which t^2 E(t) = target (traditionally 0.3)
# Consider clover discretization with and without perturbative correction
# For now we don't consider the plaquette discretization or w0 scale
# Constructing blocks based on number of measurements, not MDTU

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
# Number is hiding between tag and first '.'
cfgs = []
temp = tag + '*'
for filename in glob.glob(temp):
  start = (filename.split('.'))[0]          # Everything before first '.'
  cfgs.append(int((start.split(tag))[-1]))  # Number after tag
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
firstFile = glob.glob(tag + str(cfgs[0]) + '.*')
for line in open(firstFile[0]):
  if line.startswith('fname: '):
    temp = (line.split())[-1]       # f#l#t#b#m#
    end = (temp.split('l'))[-1]     # Stuff after 'l'
    L = int((end.split('t'))[0])
    end = (temp.split('t'))[-1]     # Stuff after 't'
    Nt = int((end.split('b'))[0])
  elif line.startswith('jname: '):
    break   # Done scanning through file
if L > Nt:
  L = Nt    # Take minimum
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Find t0 in each file -- might as well save to data/ directory
# Consider clover tSqE with and without perturbative improvement
# Skip plaquette discretization and w0
# Include bins with a few missing measurements, for consistent analyses
outfilename = 'data/scale.csv'
outfile = open(outfilename, 'w')
print >> outfile, "# t^2 E = %g" % target
print >> outfile, "# meas,sqrt(8t0)_raw,sqrt(8t0)_pert"
t_pert, pert = np.loadtxt(pertFile, unpack=True)
max_t = max(t_pert)       # Some of these measurements go to c>0.5
count = 0
missing = 0
for i in cfgs:
  toOpen = glob.glob(tag + str(i) + '.*')
  done = [0, 0]           # Stop when both are done
  t0_raw = -1.0
  t0_pert = -1.0
  index = 0
  for line in open(toOpen[0]):
    if done == [1, 1]:
      count += 1
      print >> outfile, "%d,%.8g,%.8g" % (i, t0_raw, t0_pert)
      break               # Don't bother to check file for completion
    # Format: Parameters: eps = #, NWilsonFlow = #, freq = #
    # So need to strip comma
    if line.startswith('Parameters'):    # Check epsilon
      if not dt == float(((line.split())[3]).rstrip(',')):
        print "ERROR:", toOpen[0], "uses wrong epsilon",
        print (line.split())[1]
        sys.exit(1)
    # Format: Step #: Mean Plaquette: # symmE: # symmQ: #
    # Need to strip colon from step number
    elif line.startswith('Step '):
      temp = line.split()
      t = dt * float((temp[1]).rstrip(':'))
      if t - max_t > 1e-4:          # This measurement goes beyond c=0.5
        break
      tSqEraw = t * t * np.fabs(float(temp[6]))
      # Extract and check perturbative correction for this t
      # Converting from steps to t requires floating-point comparison
      if np.fabs(t_pert[index] - t) > 1e-4:
        print "ERROR: t=%.2g in %s doesn't match t=%.2g in  %s" \
              % (t_pert[index], pertFile, t, toOpen[0])
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
print >> binfile, "# bin,start,stop,sqrt(8t0)_raw,sqrt(8t0)_pert"
missing = 0
count = 0
ave_t0raw = 0.0         # Accumulate within each block
ave_t0pert = 0.0
t0rawList = []
t0pertList = []
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
    count += 1
  elif count == num:                        # Move on to next block
    if count == missing:
      print "ERROR: All measurements in this block don't reach target"
      print "       Need to do something smarter..."
      sys.exit(1)
    t0rawList.append(ave_t0raw / float(count - missing))
    t0pertList.append(ave_t0pert / float(count - missing))
    print >> binfile, "%d,%d,%d,%.8g,%.8g" \
                      % (len(t0rawList), begin, previous, \
                         t0rawList[-1], t0pertList[-1])
    begin = i
    missing = 0
    count = 1                     # Next block begins with this line
    if 'null' in line:
      ave_t0raw = 0.0
      ave_t0pert = 0.0
      missing = 1
    else:
      ave_t0raw = float(temp[1])
      ave_t0pert = float(temp[2])
  else: # count[0] > num:                       # Sanity check
    print "ERROR: Something funny is going on..."
    sys.exit(1)
  previous = i

# Check special case that the final file fills the last block
if count == num:
  t0rawList.append(ave_t0raw / float(count - missing))
  t0pertList.append(ave_t0pert / float(count - missing))
  print >> binfile, "%d,%d,%d,%.8g,%.8g" \
                    % (len(t0rawList), begin, previous, \
                       t0rawList[-1], t0pertList[-1])
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
print "%.8g %.4g # %d" % (ave, err, N)

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
