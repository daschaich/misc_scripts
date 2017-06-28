#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
# ------------------------------------------------------------------
# This script determines the Wilson flow scales t0 and w0
# defined as the t for which t^2 E(t) = target (traditionally 0.3)
# For the former we compare both clover and plaquette discretizations
# Now constructing blocks based on number of measurements, not MDTU
# TODO: Reconstruct w0 for plaquette discretization
# TODO: Could do different targets for all three

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
# Also do plaquette discretization and (clover) w0 at the same time
# Include bins with a few missing measurements, for consistent analyses
outfilename = 'data/scale.csv'
outfile = open(outfilename, 'w')
print >> outfile, "meas,sqrt(8t0),sqrt(8plaq),w0"
count = 0
missing = 0
for i in cfgs:
  toOpen = tag + str(i)
  done = [0, 0, 0]        # Stop when all three are done
  out_t0 = -1.0
  out_plaq = -1.0
  out_w0 = -1.0
  for line in open(toOpen):
    if done == [1, 1, 1]:
      count += 1
      print >> outfile, "%d,%.8g,%.8g,%.8g" % (i, out_t0, out_plaq, out_w0)
      break               # Don't bother to check file for completion
    # Format: WFLOW t plaq E t^2*E t*d(t^2*E)/dt 12t^2*(3-plaq) top.charge
    if line.startswith('WFLOW '):
      temp = line.split()
      tSqE = np.fabs(float(temp[4]))    # fabs for 64nt128...
      deriv = np.fabs(float(temp[5]))
      plaq = np.fabs(float(temp[6]))
      if done[0] == 0 and tSqE > target:
        out_t0 = np.sqrt(8.0 * float(temp[1]))
        done[0] = 1
      if done[1] == 0 and plaq > target:
        out_plaq = np.sqrt(8.0 * float(temp[1]))
        done[1] = 1
      if done[2] == 0 and deriv > target:
        out_w0 = np.sqrt(float(temp[1]))
        done[2] = 1
    elif line.startswith('RUNNING COMPLETED'):
      if i > first:
        print "WARNING: Measurement %d never reached target %.2g" % (i, target)
      print >> outfile, i,
      if done[0] == 0:
        print >> outfile, ",null,"
      else:
        print >> outfile, ",%.8g," % out_t0
      if done[1] == 0:
        print >> outfile, "null,"
      else:
        print >> outfile, "%.8g," % out_plaq
      if done[2] == 0:
        print >> outfile, "null"
      else:
        print >> outfile, "%.8g" % out_w0
      missing += 1
outfile.close()

if not count == len(cfgs) - missing:
  print "ERROR: Have %d files but only %d measurements" % (len(cfgs), count)
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now analyze newly-created data/t0 file
missing = 0
count = 0
ave_t0 = 0.0         # Accumulate within each block
ave_plaq = 0.0
ave_w0 = 0.0
t0List = []
plaqList = []
w0List = []
begin = first     # Where each block begins, to be incremented
for line in open(outfilename):
  if line.startswith('meas'):
    continue
  temp = line.split(',')
  i = float(temp[0])
  if i < first:
    continue
  elif count < num:
    if 'null' in line:    # TODO: Could refine to handle different obs
      missing += 1
    else:
      ave_t0 += float(temp[1])
      ave_plaq += float(temp[2])
      ave_w0 += float(temp[3])
    count += 1
  elif count == num:                        # Move on to next block
    if count == missing:
      print "ERROR: All measurements in this block don't reach target"
      print "       Need to do something smarter..."
      sys.exit(1)
    t0List.append(ave_t0 / float(count - missing))
    plaqList.append(ave_plaq / float(count - missing))
    w0List.append(ave_w0 / float(count - missing))
    begin = i
    missing = 0
    count = 1                     # Next block begins with this line
    if 'null' in line:
      ave_t0 = 0.0
      ave_plaq = 0.0
      ave_w0 = 0.0
      missing = 1
    else:
      ave_t0 = float(temp[1])
      ave_plaq = float(temp[2])
      ave_w0 = float(temp[3])
  else: # count[0] > num:                       # Sanity check
    print "ERROR: Something funny is going on..."
    sys.exit(1)

# Check special case that the final file fills the last block
if count == num:
  t0List.append(ave_t0 / float(count - missing))
  plaqList.append(ave_plaq / float(count - missing))
  w0List.append(ave_w0 / float(count - missing))

# Now print mean and standard error, requiring N>1
if len(t0List) < 2:
  print "ERROR: Need multiple blocks to average"
  sys.exit(1)

print "%.2g" % target,

dat = np.array(t0List)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
print "%.8g %.4g" % (ave, err),

dat = np.array(plaqList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
print "%.8g %.4g" % (ave, err),

dat = np.array(w0List)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
print "%.8g %.4g # %d" % (ave, err, N)

runtime += time.time()
print "Runtime: %0.1f seconds" % runtime
# ------------------------------------------------------------------
