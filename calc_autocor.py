#!/usr/bin/python
import glob
import math
import os
import sys
# ------------------------------------------------------------------
# Estimate integrated autocorrelation time in MDTU
# from given dygraph data files

# tau = 0.5 + sum_{t = 1}^{t = W} Gamma(t)
# Gamma(t) = C(t) / C(0)
# C(t) = <(A_i - <A>) (A_{i + t) - <A>)>

# For now check maximum up to W=250 MDTU
# and require constant trajectory lengths and measurement frequency
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Un-normalized correlation function C(t) for data list dat
# Pre-compute mean, which isn't going to change
def C(t, dat, mean):
  N = len(dat)
  C = 0
  count = 0
  for i in range(0, N - 1 - t):
    C += (dat[i] - mean) * (dat[i + t] - mean)
    count += 1
  return (C / count)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: first is ensemble tag, second is thermalization cut
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<tag> <cut>"
  sys.exit(1)
tag = str(sys.argv[1])
cut = int(sys.argv[2])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)
if not 'Run_' in os.getcwd():
  print "ERROR: need to call from a Run_* directory"
  sys.exit(1)

# Translate MDTU cut into corresponding output file
file_cut = -1
keyfile = 'data/key.' + tag + '.csv'
for line in open(keyfile):
  if line.startswith('t'):
    continue
  temp = line.split(',')
  if float(temp[0]) > cut:
    file_cut = int(temp[1])
    break
if file_cut == -1:
  print "Skipping", tag, "since cut =", cut
  sys.exit(0)
else:
  print tag, cut

# Use that file to set the maximum number of RG blockings
# And determine the separation between measurements
blMax = 0
firstFile = 'Out/out_' + tag + '.' + str(file_cut)
for line in open(firstFile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  elif line.startswith('trajecs '):
    meas_gap = int((line.split())[1])
  elif line.startswith('traj_length '):
    t_length = float((line.split())[1])
    meas_gap *= t_length    # Convert to MDTU
  elif line.startswith('reload_serial '):
    break   # Done scanning through file

if L <= Nt:
  temp = L
else:
  temp = Nt
while temp % 2 == 0 and temp > 2:
  temp /= 2
  blMax += 1
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Plaquette is special -- two data (to be averaged) every t_length
oldt = -1
dat = []
plaqfile = 'data/plaq.' + tag + '.csv'
for line in open(plaqfile):
  temp = line.split(',')
  if temp[0].startswith('M') or float(temp[0]) < cut:
    continue
  else:
    if oldt > 0:
      sep = float(temp[0]) - oldt
      if not sep == t_length:   # Data strided by t_length
        print "ERROR: trajectory length seems non-constant:"
        print "      ", sep, "vs.", t_length
        sys.exit(1)
    oldt = float(temp[0])
    dat.append((float(temp[1]) + float(temp[2])) / 2)

# Construct C(t) -- save it in case it might be useful
N = len(dat)
ave = sum(dat) / float(N)
C0 = C(0, dat, ave)
gamma = [1]
tau = 0.5 * sep
maxtau = tau
maxW = 0
for t in range(1, int(250 / sep)):
  gamma.append(C(t, dat, ave) / C0)
#  print "Gamma(%d) = %.2g" % (t, gamma[-1]) # Check
  tau += gamma[-1] * sep
  if tau > maxtau:
    maxtau = tau
    maxW = t * sep

# Print maximum for W up to 250, normalized per MDTU
print "plaq: %.4g (W = %d)" % (maxtau, maxW)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Real part of Polyakov loop only has one datum every t_length
dat = []
polyfile = 'data/poly_r.' + tag + '.csv'
for line in open(polyfile):
  temp = line.split(',')
  if temp[0].startswith('M') or float(temp[0]) < cut:
    continue
  else:       # Constant t_length checked above
    dat.append(float(temp[1]))

# Construct C(t) -- save it in case it might be useful
N = len(dat)
ave = sum(dat) / float(N)
C0 = C(0, dat, ave)
gamma = [1]
tau = 0.5 * sep
maxtau = tau
maxW = 0
for t in range(1, int(250 / sep)):
  gamma.append(C(t, dat, ave) / C0)
#  print "Gamma(%d) = %.2g" % (t, gamma[-1]) # Check
  tau += gamma[-1] * sep
  if tau > maxtau:
    maxtau = tau
    maxW = t * sep

# Print maximum for W up to 250, normalized per MDTU
print "poly: %.4g (W = %d)" % (maxtau, maxW)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Chiral condensate has one datum every meas_gap = Ntraj * t_length
oldt = -1
dat = []
pbpfile = 'data/pbp.' + tag + '.csv'
for line in open(pbpfile):
  temp = line.split(',')
  if temp[0].startswith('M') or float(temp[0]) < cut:
    continue
  else:
    if oldt > 0:
      sep = float(temp[0]) - oldt
      if not sep == meas_gap:     # Data strided by meas_gap
        print "ERROR: measurement separation seems non-constant:"
        print "      ", sep, "vs.", meas_gap
        sys.exit(1)
    oldt = float(temp[0])
    dat.append(float(temp[1]))

# Construct C(t) -- save it in case it might be useful
N = len(dat)
ave = sum(dat) / float(N)
C0 = C(0, dat, ave)
gamma = [1]
tau = 0.5
maxtau = tau
maxW = 0
for t in range(1, int(250 / sep)):
  gamma.append(C(t, dat, ave) / C0)
#  print "Gamma(%d) = %.2g" % (t, gamma[-1]) # Check
  tau += gamma[-1] * sep
  if tau > maxtau:
    maxtau = tau
    maxW = t * sep

# Print maximum for W up to 250, normalized per MDTU
print " pbp: %.4g (W = %d)" % (maxtau, maxW)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Blocked data files have blMax + 1 data per line (including unblocked)
# We already got the "+ 1" above, but include it here as well, to check
# Just look at real parts for now (not mod, arg or imaginary parts)
for obs in ['plaqB', 'poly_rB', 'xpoly_rB']:
  oldt = -1
  dat = [[] for x in range(blMax + 1)]
  obsfile = 'data/' + obs + '.' + tag + '.csv'
  for line in open(obsfile):
    temp = line.split(',')
    if temp[0].startswith('M') or float(temp[0]) < cut:
      continue
    else:       # Check measurement separation again in case mcrg_ files missing
      if oldt > 0:
        sep = float(temp[0]) - oldt
        if not sep == meas_gap:     # Data strided by meas_gap
          print "ERROR: measurement separation seems non-constant:"
          print "      ", sep, "vs.", meas_gap
          sys.exit(1)
      oldt = float(temp[0])
      for i in range(0, blMax + 1):
        dat[i].append(float(temp[i + 1]))

  # Construct C(t) but don't bother saving it
  N = len(dat[0])
  print obs + ':',
  for i in range(0, blMax + 1):
    ave = sum(dat[i]) / float(N)
    C0 = C(0, dat[i], ave)
    tau = 0.5
    maxtau = tau
    maxW = 0
    for t in range(1, int(250 / sep)):
      tau += C(t, dat[i], ave) * sep / C0
      if tau > maxtau:
        maxtau = tau
        maxW = t * sep

    # Print maximum for W up to 250, normalized per MDTU
    if i < blMax:
      print "%.4g (W = %d)," % (maxtau, maxW),
    else:
      print "%.4g (W = %d)" % (maxtau, maxW)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Finally, topological charge only has one datum every meas_gap
# But for now there are three null entries in front of it
oldt = -1
dat = []
topofile = 'data/topo.' + tag + '.csv'
for line in open(topofile):
  temp = line.split(',')
  if temp[0].startswith('M') or float(temp[0]) < cut:
    continue
  else:   # Check measurement separation again in case Wflow_ files missing
    if oldt > 0:
      sep = float(temp[0]) - oldt
      if not sep == meas_gap:     # Data strided by meas_gap
        print "ERROR: measurement separation seems non-constant:"
        print "      ", sep, "vs.", meas_gap
        sys.exit(1)
    oldt = float(temp[0])
    dat.append(float(temp[4]))

# Construct C(t) -- save it in case it might be useful
N = len(dat)
ave = sum(dat) / float(N)
C0 = C(0, dat, ave)
gamma = [1]
tau = 0.5 * sep
maxtau = tau
maxW = 0
for t in range(1, int(250 / sep)):
  gamma.append(C(t, dat, ave) / C0)
#  print "Gamma(%d) = %.2g" % (t, gamma[-1]) # Check
  tau += gamma[-1] * sep
  if tau > maxtau:
    maxtau = tau
    maxW = t * sep

# Print maximum for W up to 250, normalized per MDTU
print "topo: %.4g (W = %d)" % (maxtau, maxW)
# ------------------------------------------------------------------
