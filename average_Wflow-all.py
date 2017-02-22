#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Construct averages and standard errors of the gradient flow coupling
# Consider only one ensemble at a time, processing all t <= L^2 / 32
# Read data/key.$tag.csv to relate Wflow files to the jackknife blocks
# specified by the given thermalization cut and block size
# Drop any partial blocks at the end

# Parse arguments: first is ensemble tag, second is thermalization cut,
# third is block size (should be larger than auto-correlation time)
# We discard any partial blocks at the end
# Fourth argument is t-shift optimization parameter tau
# Optional fifth argument tells us to use the plaquette observable
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<tag> <cut> <block> <tau> <obs>"
  sys.exit(1)
tag = str(sys.argv[1])
cut = int(sys.argv[2])
block_size = int(sys.argv[3])
tau = float(sys.argv[4])

# Choose which observable to use -- require 'plaq' as specific argument
plaq = -1
outfilename = 'results/Wflow_coupling.%s' % tag
if len(sys.argv) > 5:
  if str(sys.argv[5]) == 'plaq':
    plaq = 1
    outfilename = 'results/Wplaq_coupling.%s' % tag
  else:
    print "Warning: Implicitly using clover observable"
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
keyfile = 'data/key.' + tag + '.csv'
if not os.path.isfile(keyfile):
  print "ERROR:", keyfile, "does not exist"
  sys.exit(1)

# Extract lattice volume from a Wilson flow file that should exist
firstfile = 'Out/Wflow_' + tag + '.10'
dt = 0.01
lookup = []
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])

  # Count number of points of t >= tau
  # The perturbative correction brings in a factor of 1/(t - tau)^2
  # so we need non-zero (t - tau)...
  elif line.startswith('WFLOW '):
    temp = line.split()
    t = float(temp[1]) - tau                # Include t-shift
    if t < 0.001:                           # Start when t > tau
      continue
    lookup.append(t)
Npt = len(lookup)

# Set L to minimum of nx and nt
if Nt < L:
  L = Nt
if L < 0 or Nt < 0:
  print "ERROR: couldn't extract lattice size from", firstfile
  sys.exit(1)
target = float(L**2) / 32.0
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Compute tree-level finite-volume perturbative correction to t^2 E
# Copied from tSqE_tree.py for clover (stripping most inline documentation)
def tree(L, t):
  twopiOvL = 2.0 * np.pi / float(L)
  tSqE = 0.0
  for n1 in range(L):
    for n2 in range(L):
      for n3 in range(L):
        for n4 in range(L):
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

  tSqE += 2.0                           # Zero-mode contribution
  tSqE *= (64.0 * np.pi**2 * t**2) / (3.0 * float(L**4))
  return tSqE
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Compute tree-level finite-volume perturbative correction to t^2 E
# Copied from tSqE_plaq.py for plaq (stripping most inline documentation)
def tree_plaq(L, t):
  twopiOvL = 2.0 * np.pi / float(L)
  tSqE = 0.0
  for n1 in range(L):
    for n2 in range(L):
      for n3 in range(L):
        for n4 in range(L):
          if n1 + n2 + n3 + n4 == 0:    # Zero-mode contribution is separate
            continue

          n_mu = np.array([n1, n2, n3, n4], dtype = np.float)
          p_mu = twopiOvL * n_mu
          phat_mu = 2.0 * np.sin(p_mu / 2.0)
          phatSq = (phat_mu**2).sum()

          TrS = (phatSq - phat_mu**2).sum()

          tSqE += np.exp(-2.0 * t * phatSq) * TrS / phatSq

  tSqE += 2.0                           # Zero-mode contribution
  tSqE *= (64.0 * np.pi**2 * t**2) / (3.0 * float(L**4))
  return tSqE
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Use key file to construct jackknife blocks from Wflow files
file_tag = '-1'
ave = np.zeros(Npt, dtype = np.float)       # Accumulate within each block
count = np.zeros(Npt, dtype = np.int)
index = 0
datList = [[] for x in range(Npt)]
begin = cut       # Where each block begins, to be incremented
for traj in open(keyfile):
  if traj.startswith('t'):
    continue
  keytemp = traj.split(',')
  MDTU = float(keytemp[0])
  if MDTU <= cut or keytemp[1].rstrip() == file_tag:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
    file_tag = keytemp[1].rstrip()
    Wflowfile = 'Out/Wflow_' + tag + '.' + file_tag
    if not os.path.isfile(Wflowfile):
      print "WARNING:", Wflowfile, "does not exist"
      continue
    index = 0
    for line in open(Wflowfile):
      if line.startswith('epsilon'):    # Check epsilon
        if not dt == float((line.split())[1]):
          print "ERROR:", Wflowfile, "uses wrong epsilon",
          print (line.split())[1]
          sys.exit(1)
      elif line.startswith('tmax '):    # Check tmax
        if float((line.split())[1]) < target:
          print "ERROR:", Wflowfile, "uses too-small tmax",
          print (line.split())[1]
          sys.exit(1)
      elif line.startswith('WFLOW '):
        temp = line.split()
        t = float(temp[1]) - tau                # Include t-shift
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

    file_tag = keytemp[1].rstrip()
    Wflowfile = 'Out/Wflow_' + tag + '.' + file_tag
    if not os.path.isfile(Wflowfile):
      print "WARNING:", Wflowfile, "does not exist"
      continue
    index = 0
    for line in open(Wflowfile):
      if line.startswith('epsilon'):    # Check epsilon
        if not dt == float((line.split())[1]):
          print "ERROR:", Wflowfile, "uses wrong epsilon",
          print (line.split())[1]
          sys.exit(1)
      elif line.startswith('tmax '):    # Check tmax
        if float((line.split())[1]) < target:
          print "ERROR:", Wflowfile, "uses too-small tmax",
          print (line.split())[1]
          sys.exit(1)
      elif line.startswith('WFLOW '):
        temp = line.split()
        t = float(temp[1]) - tau                # Include t-shift
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
    pert_corr = tree(L, lookup[i])
  else:
    pert_corr = tree_plaq(L, lookup[i])
  gprop = 128.0 * 3.14159**2 / (24.0 * pert_corr)

  dat = np.array(datList[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1)
  print >> outfile, "%.4g %.8g %.4g" % (lookup[i], gprop * ave, gprop * err)
outfile.close()
# ------------------------------------------------------------------
