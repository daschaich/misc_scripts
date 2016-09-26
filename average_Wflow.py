#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Construct averages and standard errors of the gradient flow coupling
# Consider only one ensemble at a time, consider c=0.2, 0.25, 0.3 and 0.35
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
outfilename = 'results/Wflow-tau%g.%s' % (tau, tag)
if len(sys.argv) > 5:
  if str(sys.argv[5]) == 'plaq':
    plaq = 1
    col = 6
    outfilename = 'results/Wplaq-tau%g.%s' % (tau, tag)
  else:
    print "Warning: Implicitly using clover observable"
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
keyfile = 'data/key.' + tag + '.csv'
if not os.path.isfile(keyfile):
  print "ERROR:", keyfile, "does not exist"
  sys.exit(1)

# Extract lattice volume from first output file
firstfile = 'Out/out_' + tag + '.1'
if not os.path.isfile(firstfile):
  print "ERROR:", firstfile, "does not exist"
  sys.exit(1)
for line in open(firstfile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  elif line.startswith('iseed '):
    break   # Done scanning through file

# Set L to minimum of nx and nt
if Nt < L:
  L = Nt
if L < 0 or Nt < 0:
  print "ERROR: couldn't extract lattice size from", firstfile
  sys.exit(1)

# Convert c=0.2, 0.25, 0.3 and 0.35 to corresponding t
# This allows us to include the t-shift optimization
# Note that we need to include it here to keep c fixed!
magic_c = np.array([0.2, 0.25, 0.3, 0.35], dtype = np.float)
magic_t = (L * magic_c)**2 / 8.0
#print magic_t

# This should be a reasonably common epsilon
# We will just rewrite this for each measurement
dt = 0.01
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Use key file to construct jackknife blocks from Wflow files
file_tag = '-1'
ave = np.zeros_like(magic_c)  # Accumulate within each block
count = np.zeros_like(magic_c, dtype = np.int)
datList = [[] for x in range(len(magic_c))]
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
    for line in open(Wflowfile):
      if line.startswith('epsilon'):
        dt = float((line.split())[1])
      if line.startswith('WFLOW '):
        temp = line.split()
        t = float(temp[1]) - tau    # Include t-shift
        for i in range(len(magic_c)):
          if t - dt + 1.0e-4 < magic_t[i] and t + 1.0e-4 >= magic_t[i]:
            ave[i] += t**2 * float(temp[3])    # t**2 * E(t + tau)
            count[i] += 1

  elif MDTU >= (begin + block_size):    # Move on to next block
    for i in range(len(magic_c)):
      if count[i] >= 1:
        datList[i].append(ave[i] / float(count[i]))
      if count[i] != count[0]:
        print "WARNING: measurement mismatch %d vs %d for i=%d in block %d" \
          % (count[i], count[0], i, len(datList[0]))
    ave = np.zeros_like(magic_c)
    count = np.zeros_like(magic_c, dtype = np.int)
    begin += block_size

    file_tag = keytemp[1].rstrip()
    Wflowfile = 'Out/Wflow_' + tag + '.' + file_tag
    if not os.path.isfile(Wflowfile):
      print "WARNING:", Wflowfile, "does not exist"
      continue
    for line in open(Wflowfile):
      if line.startswith('epsilon'):
        dt = float((line.split())[1])
      elif line.startswith('WFLOW '):
        temp = line.split()
        t = float(temp[1]) - tau    # Include t-shift
        for i in range(len(magic_c)):
          if t - dt + 1.0e-4 < magic_t[i] and t + 1.0e-4 >= magic_t[i]:
            if plaq < 0:
              ave[i] += t**2 * float(temp[3])    # t**2 * E(t + tau)
            else:
              # Data is (t + tau)**2 * E(t + tau)...
              rescale = t / float(temp[1])
              ave[i] += rescale**2 * float(temp[col])
            count[i] += 1

# The file_tag check in the loop above means that
# full blocks aren't recorded until we reach the next file
# Check special case that the final file fills the last block
if MDTU >= (begin + block_size) and count[0] > 1:
  for i in range(len(magic_c)):
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
for i in range(len(magic_c)):
  if len(datList[i]) != len(datList[0]):
    print "ERROR: measurement mismatch %d vs %d for i=%d" \
          % (len(datList[i]), len(datList[0]), i)
    sys.exit(1)
if N <= 1:
  print "ERROR: Not enough blocks (%d) to average" % N
  sys.exit(1)

# Include perturbative finite-volume + zero-mode corrections
gprop = np.zeros_like(magic_c)
for i in range(len(magic_c)):
  pert_file = "../scripts/Ct_pert_c" + str(magic_c[i])
  if not os.path.isfile(pert_file):
    print "WARNING:", pert_file, "not found"
    sys.exit(1)
  for line in open(pert_file):
    temp = line.split()
    if temp[1] == str(L):
      gprop[i] = 128.0 * 3.14159**2 / (24.0 * float(temp[7]))
      break
outfile = open(outfilename, 'w')
print >> outfile, "# tau=%g" % tau
print >> outfile, "# c=%g err c=%g err c=%g err c=%g err" \
                  % (magic_c[0], magic_c[1], magic_c[2], magic_c[3])
for i in range(len(magic_c)):
  dat = np.array(datList[i])
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1)
  print >> outfile, "%.8g %.4g" % (gprop[i] * ave, gprop[i] * err),
print >> outfile, "#", int(N)
outfile.close()
# ------------------------------------------------------------------
