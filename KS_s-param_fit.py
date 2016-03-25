#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Run jackknifed fits of V-A polarization function Pi(Q^2)
# to (1, 2) rational function

# Staggered vs. domain wall fermions affect renormalization factor
# (here Z=1) and file names

# Parse arguments: first is number of directories,
# each of which much then be listed along with a thermalization cut
# The last argument is the fixed block size
# Follow Meifeng by omitting the partial block at the end
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<# of dirs>",
  print "<name> and <cut> for each dir",
  print "<block>"
  sys.exit(1)
Ndirs = int(sys.argv[1])
dirnames = []
cuts = []
for i in range(Ndirs):
  dirnames.append(str(sys.argv[2 * i + 2]))
  cuts.append(int(sys.argv[2 * i + 3]))
block_size = int(sys.argv[len(sys.argv) - 1])
runtime = -time.time()

# !!! Number of points to use in fit
Npts = 21     # Get up to QSq=0.4 for 32nt64
#Npts = 11     # Get up to QSq=0.4 for 24nt48

# errfunc will be minimized via least-squares optimization
pade12 = lambda p, x: (p[0] + p[1] * x) / (1 + x * (p[2] + x * p[3]))
errfunc = lambda p, x, y, err: (pade12(p, x) - y) / err
p_in = [-0.01, -0.01, -0.1, -0.1]   # Order-of-magnitude initial guesses

# !!! This could make it easier to check for bad poles...
#pade12 = lambda p, x: (p[0] + p[1] * x) / ((x + p[2]) * (x + p[3]))
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
for dirname in dirnames:
  if not os.path.isdir(dirname):
    print "ERROR:", dirname, "does not exist"
    sys.exit(1)

# Construct lists of which sources have been used on which configurations
cfgs = []
meas = [[] for x in range(Ndirs)]
tot_meas = 0
for i in range(Ndirs):
  cfgs.append([])
  files = dirnames[i] + '/KSdecomp.*'
  for filename in glob.glob(files):
    cfg = int((filename.split('.'))[-1])  # Number after last .
    if cfg not in cfgs[i] and cfg >= cuts[i]:
      cfgs[i].append(cfg)
  cfgs[i].sort()

  if len(cfgs[i]) == 0:
    print "ERROR: no files", files, "found"
    sys.exit(1)

  # Now that we have the configurations sorted,
  # construct the corresponding lists of sources
  check = -1
  for j in cfgs[i]:
    meas[i].append([])
    files = dirnames[i] + '/KSdecomp.t*.' + str(j)
    for filename in glob.glob(files):
      temp = filename.split('KSdecomp.t')[1] # Everything after KSdecomp.t
      t_src = int((temp.split('.'))[0])    # Number after '.t'
      meas[i][-1].append(t_src)
    if check > 0:
      # Check whether any configurations appear to be missing
      if (j - prev) != check:
        print "WARNING: spacing between %d and %d isn't %d" \
              % (prev, j, check)
    else:
      meas_check = len(meas[i][0])
      check = cfgs[i][1] - cfgs[i][0]
    prev = j

    # Check whether any configurations seem to be missing measurements
    temp = len(meas[i][-1])
    if temp != meas_check:
      print "WARNING: %d measurements for configuration %d instead of %d" \
            % (temp, i, meas_check)
    tot_meas += temp

tot_cfgs = sum([len(cfgs[i]) for i in range(Ndirs)])
print "# In total %d measurements on %d configurations:" \
      % (tot_meas, tot_cfgs),
if Ndirs == 1:
  print "(%d--%d by %d)" % (min(cfgs[i]), max(cfgs[i]), check)
else:
  print "(%d--%d," % (min(cfgs[0]), max(cfgs[0])),
  for i in range(1, Ndirs - 1):
    print "%d--%d," % (min(cfgs[i]), max(cfgs[i])),
  print "%d--%d by %d)" % (min(cfgs[Ndirs - i]), max(cfgs[Ndirs - 1]), check)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct array of blocked measurements
dat = [[] for x in range(Npts)]
datSq = [[] for x in range(Npts)]
lengths = []      # Number of measurements in each jackknife block

# First save a reference list of the QSq we should find below
filename = dirnames[0] + '/KSdecomp.t' + str(meas[0][0][0]) + '.' + str(cfgs[0][0])
temp = []
for line in open(filename):
  temp.append(float((line.split('\t'))[0]))
  if len(temp) == Npts:
    break
QSq = np.array(temp)

for i in range(Ndirs):  # Reset block for each directory
  count = 0
  begin = cuts[i]       # Where each block begins, to be incremented
  toAve = [0 for x in range(Npts)]
  toAveSq = [0 for x in range(Npts)]
  j = -1
  for MDTU in cfgs[i]:
    j += 1
    for src in meas[i][j]:
      # If we're done with this block, record it and reset for the next
      if MDTU >= (begin + block_size):
        for Q in range(Npts):
          dat[Q].append(toAve[Q] / float(count))
          datSq[Q].append(toAveSq[Q] / float(count))

        lengths.append(count)
        # Print summary of blocks
        print "# Block %d holds %d measurements from MDTU in [%d, %d)" \
              % (len(lengths), lengths[-1], begin, begin + block_size)
        count = 0
        begin += block_size
        toAve = [0 for x in range(Npts)]
        toAveSq = [0 for x in range(Npts)]
        # Handle missing measurements -- skip to next block at MDTU
        if MDTU >= (begin + block_size):
          begin = MDTU

      # Running averages
      # Real part of transverse V-A vacuum polarization is 12th datum of 13
      # !!! Scale by -0.5 to roughly match DWF results
      x = -1
      filename = dirnames[i] + '/KSdecomp.t' + str(src) + '.' + str(MDTU)
      for line in open(filename):
        x += 1
        if x == Npts:
          break
        temp = line.split('\t')
        if float(temp[0]) != QSq[x]:
          print "ERROR: QSq mismatch in file %s: %g vs. %g" \
                % (filename, float(temp[0]), QSq[x])
          sys.exit(1)
        toAdd = float(temp[11]) / 2.0
        toAdd *= -0.5           # !!!
        toAve[x] += toAdd
        toAveSq[x] += toAdd**2
      count += 1

  # Check special case that the loop above ended
  # before the last full block was recorded
  if cfgs[i][-1] >= (begin + block_size - cfgs[i][-1] + cfgs[i][-2]):
    for Q in range(Npts):
      dat[Q].append(toAve[Q] / float(count))
      datSq[Q].append(toAveSq[Q] / float(count))
    lengths.append(count)

Nblocks = len(dat[0])

# If we're doing a single fit with propagated uncertainties,
# make sure the data is actually set up to fit
if Nblocks == 0:
  for Q in range(Npts):
    dat[Q].append(toAve[Q] / float(count))
    datSq[Q].append(toAveSq[Q] / float(count))
  lengths.append(count)
  Nblocks = 1
  # Print summary of blocks
  print "# Block %d holds %d measurements from MDTU in [%d, %d)" \
      % (len(lengths), lengths[-1], begin, begin + block_size)

print "# Total of %d measurements in %d block(s) of length %d" \
    % (sum(lengths), Nblocks, block_size)
print "# Using %d points, up to Q^2 = %.4g" % (Npts, QSq[Npts - 1])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# If only one block, just do a single fit with propagated uncertainties
# In this case only count is the total number of measurements
# optimize.leastsq does not scale by chiSq_dof, which we desire
if Nblocks == 1:
  Pi = np.array([dat[i][0] for i in range(Npts)])
  PiErr = np.empty(Npts, dtype = np.float)
  for i in range(Npts):
    err = datSq[i][0] - dat[i][0]**2
    PiErr[i] = np.sqrt(err / (float(count) - 1.))
#    print QSq[i], Pi[i], PiErr[0]
  all_out = optimize.leastsq(errfunc, p_in[:], args=(QSq, Pi, PiErr),
                             full_output = 1)
  p_out = all_out[0]
  covar = all_out[1]

  # Sanity check
#  print p_out
  if p_out[0] > 0:
    print "Ack!  Fitter found wrong-sign pole"
    sys.exit(1)

  # Use covariance matrix to propagate uncertainties
  # derivs are derivatives of slope with each of p_out[:]
  slope = 4 * np.pi * (p_out[1] - p_out[0] * p_out[2])
  derivs = 4 * np.pi * np.array([-p_out[2], 1, -p_out[0], 0])
  err = np.sqrt(np.dot(derivs, np.dot(covar, derivs)))
  print "slope %.6g %.4g" % (slope, err)

  FP = np.sqrt(-1. * p_out[0])
  derivs = np.array([-1  / (2 * FP), 0, 0, 0])
  err = np.sqrt(np.dot(derivs, np.dot(covar, derivs)))
  print "FP %.6g %.4g" % (FP, err)
  runtime += time.time()
  print "# Runtime: %.2g seconds" % runtime
  sys.exit(0)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can construct jackknife samples through single-block elimination
tot = np.array([sum(dat[i]) for i in range(Npts)])
totSq = np.array([sum(datSq[i]) for i in range(Npts)])

# Fit results for all jk samples
p_out = np.empty((len(p_in), Nblocks), dtype = np.float)
jkslopes = np.empty(Nblocks, dtype = np.float)
jkFP = np.empty(Nblocks, dtype = np.float)
for i in range(Nblocks):  # Jackknife samples
  Pi = np.empty(Npts, dtype = np.float)
  PiErr = np.empty(Npts, dtype = np.float)
  for Q in range(Npts):
    Pi[Q] = (tot[Q] - dat[Q][i]) / (Nblocks - 1.)
    temp = (totSq[Q] - datSq[Q][i]) / (Nblocks - 1.)
    PiErr[Q] = np.sqrt(temp - Pi[Q]**2)
  temp, success = optimize.leastsq(errfunc, p_in[:], args=(QSq, Pi, PiErr))
  for j in range(len(p_in)):
    p_out[j][i] = temp[j]
  jkslopes[i] = 4. * np.pi * (p_out[1][i] - p_out[0][i] * p_out[2][i])
  jkFP[i] = np.sqrt(-1. * p_out[0][i])

  # Sanity check
  if p_out[0][i] > 0:
    print "Ack!  Fitter found wrong-sign pole"
    sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
# Try estimating the full error matrix
p_ave = np.empty(len(p_in), dtype = np.float)
cov = np.empty((len(p_in), len(p_in)), dtype = np.float)
for i in range(len(p_in)):
  p_ave[i] = np.average(p_out[i])
for i in range(len(p_in)):
  for j in range(i, len(p_in)):
    temp = np.sum((p_out[i] - p_ave[i]) * (p_out[j] - p_ave[j]))
    cov[i][j] = (Nblocks - 1.) * temp / float(Nblocks)
    if j > i:
      cov[j][i] = cov[i][j]

# Slope
slope = 4. * np.pi * (p_ave[1] - p_ave[0] * p_ave[2])
derivs = 4. * np.pi * np.array([-p_ave[2], 1, -p_ave[0], 0])
err = np.sqrt(np.dot(derivs, np.dot(cov, derivs)))
print "slope %.6g %.4g" % (slope, err)

ave = np.average(jkslopes)
var = (Nblocks - 1.) * np.sum((jkslopes - ave)**2) / float(Nblocks)
print "check %.6g %.4g" % (ave, np.sqrt(var))

# Intercept
FP = np.sqrt(-1. * p_ave[0])
derivs = np.array([-1  / (2 * FP), 0, 0, 0])
err = np.sqrt(np.dot(derivs, np.dot(cov, derivs)))
print "FP %.6g %.4g" % (FP, err)

ave = np.average(jkFP)
var = (Nblocks - 1.) * np.sum((jkFP - ave)**2) / float(Nblocks)
print "check %.6g %.4g" % (ave, np.sqrt(var))

runtime += time.time()
print "# Runtime: %.2g seconds" % runtime
# ------------------------------------------------------------------
