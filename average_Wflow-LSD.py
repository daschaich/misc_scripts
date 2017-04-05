#!/usr/bin/python
import os
import sys
import time
import glob
import numpy as np
# ------------------------------------------------------------------
# Construct averages and standard errors of the gradient flow coupling
# for given ensemble, processing all t <= L^2 / 32
# specified by the given thermalization cut and block size
# Drop any partial blocks at the end

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than auto-correlation time)
# We discard any partial blocks at the end
# Third is base name of files to analyze, including directory
# Fourth argument is t-shift optimization parameter tau
# Optional fifth argument tells us to use the plaquette observable
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <dir/tag> <tau> <obs>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
tag = str(sys.argv[3])
tau = float(sys.argv[4])
runtime = -time.time()

# Choose which observable to use -- require 'plaq' as specific argument
# Will probably have tau=0, so don't include in output file name
plaq = -1
outfilename = 'results/Wflow_coupling.dat'
if len(sys.argv) > 5:
  if str(sys.argv[5]) == 'plaq':
    plaq = 1
    outfilename = 'results/Wplaq_coupling.dat'
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
    Ns = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
    # Set L to minimum of nx and nt
    if Nt < Ns:
      L = Nt
    else:
      L = Ns
    if L < 0 or Nt < 0:
      print "ERROR: couldn't extract lattice size from", firstFile
      sys.exit(1)
    target = float(L**2) / 32.0

  # Count number of points of t >= tau
  # The perturbative correction brings in a factor of 1/(t - tau)^2
  # so we need non-zero (t - tau)...
  # We also need to stop at target even if the file has more data
  elif line.startswith('WFLOW '):
    temp = line.split()
    t = float(temp[1])
    if t > target:
      break
    t -= tau                    # Include t-shift
    if t < 0.001:               # Start when t > tau
      continue
    lookup.append(t)
Npt = len(lookup)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Compute tree-level finite-volume perturbative correction to t^2 E
# Copied from tSqE_tree.py for clover (stripping most inline documentation)
def tree(Ns, Nt, t):
  one_ov_L = np.array([1.0 / float(Ns), 1.0 / float(Ns),
                       1.0 / float(Ns), 1.0 / float(Nt)], dtype = np.float)
  twopiOvL = 2.0 * np.pi * one_ov_L
  tSqE = 0.0
  for n1 in range(Ns):
    for n2 in range(Ns):
      for n3 in range(Ns):
        for n4 in range(Nt):
          if n1 + n2 + n3 + n4 == 0:    # Zero-mode contribution is separate
            continue

          n_mu = np.array([n1, n2, n3, n4], dtype = np.float)
          p_mu = twopiOvL * n_mu
          phat_mu = 2.0 * np.sin(p_mu / 2.0)
          ptw_mu = np.sin(p_mu)

          phatSq = (phat_mu**2).sum()
          ptwSq = (ptw_mu**2).sum()

          TrS = ((ptwSq - ptw_mu**2) * (np.cos(p_mu / 2))**2).sum()

          tSqE += np.exp(-2.0 * t * phatSq) * TrS / phatSq

  print twopiOvL, n_mu, p_mu # !!!TODO
  sys.exit(0)
  tSqE += 2.0                           # Zero-mode contribution
  tSqE *= (64.0 * np.pi**2 * t**2) / (3.0 * float(Ns**3 * Nt))
  return tSqE
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Compute tree-level finite-volume perturbative correction to t^2 E
# Copied from tSqE_plaq.py for plaq (stripping most inline documentation)
def tree_plaq(Ns, Nt, t):
  one_ov_L = np.array([1.0 / float(Ns), 1.0 / float(Ns),
                       1.0 / float(Ns), 1.0 / float(Nt)], dtype = np.float)
  twopiOvL = 2.0 * np.pi * one_ov_L
  tSqE = 0.0
  for n1 in range(Ns):
    for n2 in range(Ns):
      for n3 in range(Ns):
        for n4 in range(Ns):
          if n1 + n2 + n3 + n4 == 0:    # Zero-mode contribution is separate
            continue

          n_mu = np.array([n1, n2, n3, n4], dtype = np.float)
          p_mu = twopiOvL * n_mu
          phat_mu = 2.0 * np.sin(p_mu / 2.0)
          phatSq = (phat_mu**2).sum()

          TrS = (phatSq - phat_mu**2).sum()

          tSqE += np.exp(-2.0 * t * phatSq) * TrS / phatSq

  tSqE += 2.0                           # Zero-mode contribution
  tSqE *= (64.0 * np.pi**2 * t**2) / (3.0 * float(Ns**3 * Nt))
  return tSqE
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct jackknife blocks from Wflow files
ave = np.zeros(Npt, dtype = np.float)       # Accumulate within each block
count = np.zeros(Npt, dtype = np.int)
index = 0
datList = [[] for x in range(Npt)]
begin = cut       # Where each block begins, to be incremented
for MDTU in cfgs:
  if MDTU <= cut:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
    WflowFile = tag + str(MDTU)
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
        else:
          # Data is (t + tau)**2 * E(t + tau)...
          rescale = t / float(temp[1])
          ave[index] += rescale**2 * float(temp[6])
        count[index] += 1
        index += 1
    # Check that correct number of points were read
    if not index == Npt:
      print "ERROR: Only", index, "of", Npt, "points counted"
      sys.exit(1)

  elif MDTU >= (begin + block_size):    # Move on to next block
    for i in range(Npt):
      if count[i] >= 1:
        datList[i].append(ave[i] / float(count[i]))
      if count[i] != count[0]:
        print "WARNING: measurement mismatch %d vs %d for i=%d in block %d" \
          % (count[i], count[0], i, len(datList[0]))
    ave = np.zeros(Npt, dtype = np.float)
    count = np.zeros(Npt, dtype = np.int)
    begin += block_size

    WflowFile = tag + str(MDTU)
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
        else:
          # Data is (t + tau)**2 * E(t + tau)...
          rescale = t / float(temp[1])
          ave[index] += rescale**2 * float(temp[6])
        count[index] += 1
        index += 1
    # Check that correct number of points were read
    if not index == Npt:
      print "ERROR: Only", index, "of", Npt, "points counted"
      sys.exit(1)

# The file_tag check in the loop above means that
# full blocks aren't recorded until we reach the next file
# Check special case that the final file fills the last block
if MDTU >= (begin + block_size) and count[0] > 1:
  for i in range(Npt):
    if count[i] >= 1:
      datList[i].append(ave[i] / float(count[i]))
    if count[i] != count[0]:
      print "WARNING: measurement mismatch %d vs %d for i=%d in block %d" \
          % (count[i], count[0], i, len(datList[0]))
  begin += block_size
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now print mean and standard error, requiring N>1
N = float(len(datList[0]))
for i in range(Npt):
  if len(datList[i]) != len(datList[0]):
    print "ERROR: measurement mismatch %d vs %d for i=%d" \
          % (len(datList[i]), len(datList[0]), i)
    sys.exit(1)
if N <= 1:
  print "ERROR: Not enough blocks (%d) to average" % N
  sys.exit(1)

outfile = open(outfilename, 'w')
print >> outfile, "# tau=%g with %d blocks" % (tau, int(N))

# Include perturbative finite-volume + zero-mode corrections
for i in range(Npt):
  if plaq < 0:
    pert_corr = tree(Ns, Nt, lookup[i])
  else:
    pert_corr = tree_plaq(Ns, Nt, lookup[i])
  gprop = 128.0 * 3.14159**2 / (24.0 * pert_corr)

  dat = np.array(datList[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1)
  print >> outfile, "%.4g %.8g %.4g" % (lookup[i], gprop * ave, gprop * err)

runtime += time.time()
print >> outfile, "# Runtime: %d seconds" % int(runtime)
outfile.close()
# ------------------------------------------------------------------
