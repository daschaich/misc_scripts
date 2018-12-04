#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph data files to construct averages and standard errors
# with given thermalization cut and block size

# Assume one ensemble per directory (overwrite results files)

# Also accumulate eigenvalues from eig.* for rho(lambda)
# TODO: Account for block size in eigenvalue accumulation
MAX = 99    # Hard code number of eigenvalues to extract
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: first is thermalization cut,
# second is block size (should be larger than auto-correlation time)
# We discard any partial blocks at the end
# I also need to specify where I took over the run
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<cut> <block> <start>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
start = int(sys.argv[3])

# First extract lattice volume from path
# For now handle only one- or two-digit L and Nt
path = os.getcwd()
temp = path.split('nt')
if '/' in temp[0][-2:]:
  L = int(temp[0][-1:])    # Last digit before 'nt'
else:
  L = int(temp[0][-2:])    # Last two digits before 'nt'
if '/' in temp[0][:2]:
  Nt = int(temp[1][:1])    # First digit after 'nt'
else:
  Nt = int(temp[1][:2])    # First two digits after 'nt'
vol = L**3 * Nt

# Extract number of flavors for pbp normalization
# For quenched "valence pbp" we want the 4f normalization
if '4f' in path or '0f' in path:
  pbp_norm = 1.0
elif '8f' in path:
  pbp_norm = 2.0
else:
  print "ERROR: So far only 4f and 8f set up"
  sys.exit(1)

# Set the maximum number of RG blockings
blMax = 0
if L <= Nt:
  temp = L
else:
  temp = Nt
while temp % 2 == 0 and temp > 2:
  temp /= 2
  blMax += 1

# Check that we actually have data to average
# and convert thermalization cut from MDTU to trajectory number
# (Note unit-length trajectories let us reuse block_size for t_block)
MDTUfile = 'data/TU.csv'
sav = 0
good = -1
for line in open(MDTUfile):
  if line.startswith('t'):
    continue
  temp = line.split(',')
  if float(temp[1]) > cut:
    good = 1
    t_cut = sav
    break
  sav = float(temp[0])

final_MDTU = float(temp[1])
if good == -1:
  print "Error: no data to analyze",
  print "since cut=%d but we only have %d MDTU" % (cut, final_MDTU)
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Plaquette is special -- two data (to be averaged) per line
count = 0
ave = 0.          # Accumulate within each block
datList = []
begin = cut       # Where each block begins, to be incremented
plaqfile = 'data/plaq.csv'
for line in open(plaqfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0]) + start
  if MDTU <= cut:
    continue
  elif MDTU >= begin and float(temp[0]) < (begin + block_size):
    ave += (float(temp[1]) + float(temp[2])) / 2
    count += 1
  elif MDTU >= (begin + block_size):  # Move on to next block
    datList.append(ave / float(count))
    begin += block_size
    count = 1                     # Next block begins with this line
    ave = (float(temp[1]) + float(temp[2])) / 2

# Now print mean and standard error, assuming N>1
dat = np.array(datList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.)
outfilename = 'results/plaq.dat'
outfile = open(outfilename, 'w')
print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Chiral condensate is also special
# It needs to be halved to be normalized per continuum flavor
count = 0
ave = 0.          # Accumulate within each block
datList = []
begin = cut       # Where each block begins, to be incremented
pbpfile = 'data/pbp.csv'
for line in open(pbpfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0]) + start
  if MDTU <= cut:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
    ave += float(temp[1])
    count += 1
  elif MDTU >= (begin + block_size):  # Move on to next block
    datList.append(ave / float(pbp_norm * count))
    begin += block_size
    count = 1                     # Next block begins with this line
    ave = float(temp[1])

# Now print mean and standard error, assuming N>1
dat = np.array(datList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
outfilename = 'results/pbp.dat'
outfile = open(outfilename, 'w')
print >> outfile, "# Normalized per continuum flavor"
print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
outfile.close()
# ----------------------------------------------------------------



# ------------------------------------------------------------------
# Polyakov loop and spatial Wilson line have only one datum per line
# Just look at real parts for now (not mod, arg or imaginary parts)
for obs in ['poly_r', 'xpoly_r']:
  count = 0
  ave = 0.0         # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0]) + start
    if MDTU <= cut:
      continue
    elif MDTU >= begin and MDTU < (begin + block_size):
      ave += float(temp[1])
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      datList.append(ave / float(count))
      begin += block_size
      count = 1                     # Next block begins with this line
      ave = float(temp[1])

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList)
  N = np.size(dat)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For algorithmic/cost quantities
# we're interested in the first datum on each line
# Need to work in terms of trajectories rather than MDTU
for obs in ['wallTU', 'cg_iters', 'accP', 'exp_dS']:
  skip = -1
  count = 0
  ave = 0.0         # Accumulate within each block
  datList = []
  begin = t_cut     # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M') or line.startswith('t'):
      continue
    temp = line.split(',')
    traj = float(temp[0])
    if traj <= t_cut:
      continue
    elif traj > begin and traj < (begin + block_size):
      ave += float(temp[1])
      count += 1
    elif traj >= (begin + block_size):      # Move on to next block
      if count == 0:
        print "WARNING: no %s data to average at %d traj" % (obs, int(traj))
        skip = 1
        break
      datList.append(ave / float(count))
      begin += block_size
      count = 1                             # Next block begins here
      ave = float(temp[1])

  if len(datList) == 0:
    skip = 1
  if skip > 0:
    continue

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList)
  N = np.size(dat)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Blocked data files have blMax + 1 data per line (including unblocked)
# We already got the "+ 1" above, but include it here as well, to check
# Just look at real parts for now (not mod, arg or imaginary parts)
for obs in ['plaqB', 'poly_rB', 'xpoly_rB']:
  skip = -1
  count = 0
  ave = [0.0 for x in range(blMax + 1)] # Accumulate within each block
  datList = [[] for x in range(blMax + 1)]
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0]) + start
    if MDTU <= cut:
      continue
    elif MDTU >= begin and MDTU < (begin + block_size):
      for i in range(0, blMax + 1):
        ave[i] += float(temp[i + 1])
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      if count == 0:
        print "WARNING: no %s data to average at %d traj" % (obs, int(traj))
        skip = 1
        break
      for i in range(0, blMax + 1):
        datList[i].append(ave[i] / float(count))
        ave[i] = float(temp[i + 1])
      begin += block_size
      count = 1                     # Next block begins with this line

  if len(datList[0]) == 0:
    skip = 1
  if skip > 0:
    continue

  # Now print mean and standard error, assuming N>1
  outfilename = 'results/' + obs +  '.dat'
  outfile = open(outfilename, 'w')
  for i in range(0, blMax + 1):
    dat = np.array(datList[i])
    N = np.size(dat)
    ave = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
    print >> outfile, "%.8g %.4g" % (ave, err),
  print >> outfile, "# %d" % N
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the plaquette and link differences, the last datum on each line
# is the euclidean norm I'm interested in for now
# Also need to normalize the link differences by (half) the volume
for obs in ['plaq_diff', 'link_diff']:
  skip = -1
  count = 0
  ave = 0.          # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0]) + start
    if MDTU <= cut:
      continue
    elif MDTU >= begin and MDTU < (begin + block_size):
      if obs == 'plaq_diff':
        ave += float(temp[5])
      elif obs == 'link_diff':
        ave += float(temp[5]) * 2. / float(vol)
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      if count == 0:
        print "WARNING: no %s data to average at %d traj" % (obs, int(traj))
        skip = 1
        break
      datList.append(ave / float(count))
      begin += block_size
      count = 1                     # Next block begins with this line
      if obs == 'plaq_diff':
        ave = float(temp[5])
      elif obs == 'link_diff':
        ave = float(temp[5]) * 2.0 / float(vol)

  if len(datList) == 0:
    skip = 1
  if skip > 0:
    continue

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList)
  N = np.size(dat)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
  outfile.close()
# ----------------------------------------------------------------



# ----------------------------------------------------------------
# For the Wilson flow, the fourth datum on each line is the c=0.3
# running coupling I'm interested in for now
# This still needs the finite volume volume correction (1+delta)~0.97
count = 0
ave = 0.          # Accumulate within each block
datList = []
begin = cut       # Where each block begins, to be incremented
vol_corr = 0.97
flowfile = 'data/Wflow.csv'
for line in open(flowfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0]) + start
  if MDTU <= cut:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
    ave += float(temp[3]) / vol_corr
    count += 1
  elif MDTU >= (begin + block_size):  # Move on to next bloc
    datList.append(ave / float(count))
    begin += block_size
    count = 1                     # Next block begins with this line
    ave = float(temp[3]) / vol_corr

# Now print mean and standard error, assuming N>1
dat = np.array(datList)
N = np.size(dat)
ave = np.mean(dat, dtype = np.float64)
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.)
outfilename = 'results/Wflow_gSq.dat'
outfile = open(outfilename, 'w')
print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
outfile.close()
# ------------------------------------------------------------------



# ----------------------------------------------------------------
# For the Wilson-flowed Polyakov loop,
# let's print out all four of c=0.2, 0.3, 0.4 and 0.5
count = 0
ave2 = 0.0        # Accumulate within each block
ave3 = 0.0
ave4 = 0.0
ave5 = 0.0
datList2 = []
datList3 = []
datList4 = []
datList5 = []
begin = cut       # Where each block begins, to be incremented
flowfile = 'data/Wpoly.csv'
for line in open(flowfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0]) + start
  if MDTU <= cut:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
    ave2 += float(temp[1])
    ave3 += float(temp[2])
    ave4 += float(temp[3])
    ave5 += float(temp[4])
    count += 1
  elif MDTU >= (begin + block_size):  # Move on to next bloc
    datList2.append(ave2 / float(count))
    datList3.append(ave3 / float(count))
    datList4.append(ave4 / float(count))
    datList5.append(ave5 / float(count))
    begin += block_size
    count = 1                     # Next block begins with this line
    ave2 = float(temp[1])
    ave3 = float(temp[2])
    ave4 = float(temp[3])
    ave5 = float(temp[4])

# Now print mean and standard error, assuming N>1
dat2 = np.array(datList2)
dat3 = np.array(datList3)
dat4 = np.array(datList4)
dat5 = np.array(datList5)
N = np.size(dat2)
ave2 = np.mean(dat2, dtype = np.float64)
err2 = np.std(dat2, dtype = np.float64) / np.sqrt(N - 1.)
ave3 = np.mean(dat3, dtype = np.float64)
err3 = np.std(dat3, dtype = np.float64) / np.sqrt(N - 1.)
ave4 = np.mean(dat4, dtype = np.float64)
err4 = np.std(dat4, dtype = np.float64) / np.sqrt(N - 1.)
ave5 = np.mean(dat5, dtype = np.float64)
err5 = np.std(dat5, dtype = np.float64) / np.sqrt(N - 1.)
outfilename = 'results/Wpoly.dat'
outfile = open(outfilename, 'w')
print >> outfile, "# c=0.2 err c=0.3 err c=0.4 err c=0.5 err # Nblocks"
print >> outfile, "%.8g %.4g %.8g %.4g" % (ave2, err2, ave3, err3),
print >> outfile, "%.8g %.4g %.8g %.4g # %d" % (ave4, err4, ave5, err5, N)
outfile.close()
# ------------------------------------------------------------------



# ----------------------------------------------------------------
# For the Wilson flow anisotropy E_ss / E_st,
# again print out all four of c=0.2, 0.3, 0.4 and 0.5
count = 0
ave2 = 0.0        # Accumulate within each block
ave3 = 0.0
ave4 = 0.0
ave5 = 0.0
datList2 = []
datList3 = []
datList4 = []
datList5 = []
begin = cut       # Where each block begins, to be incremented
flowfile = 'data/Wflow_aniso.csv'
for line in open(flowfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0]) + start
  if MDTU <= cut:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
    ave2 += float(temp[1])
    ave3 += float(temp[2])
    ave4 += float(temp[3])
    ave5 += float(temp[4])
    count += 1
  elif MDTU >= (begin + block_size):  # Move on to next bloc
    datList2.append(ave2 / float(count))
    datList3.append(ave3 / float(count))
    datList4.append(ave4 / float(count))
    datList5.append(ave5 / float(count))
    begin += block_size
    count = 1                     # Next block begins with this line
    ave2 = float(temp[1])
    ave3 = float(temp[2])
    ave4 = float(temp[3])
    ave5 = float(temp[4])

# Now print mean and standard error, assuming N>1
dat2 = np.array(datList2)
dat3 = np.array(datList3)
dat4 = np.array(datList4)
dat5 = np.array(datList5)
N = np.size(dat2)
ave2 = np.mean(dat2, dtype = np.float64)
err2 = np.std(dat2, dtype = np.float64) / np.sqrt(N - 1.)
ave3 = np.mean(dat3, dtype = np.float64)
err3 = np.std(dat3, dtype = np.float64) / np.sqrt(N - 1.)
ave4 = np.mean(dat4, dtype = np.float64)
err4 = np.std(dat4, dtype = np.float64) / np.sqrt(N - 1.)
ave5 = np.mean(dat5, dtype = np.float64)
err5 = np.std(dat5, dtype = np.float64) / np.sqrt(N - 1.)
outfilename = 'results/Wflow_aniso.dat'
outfile = open(outfilename, 'w')
print >> outfile, "# c=0.2 err c=0.3 err c=0.4 err c=0.5 err # Nblocks"
print >> outfile, "%.8g %.4g %.8g %.4g" % (ave2, err2, ave3, err3),
print >> outfile, "%.8g %.4g %.8g %.4g # %d" % (ave4, err4, ave5, err5, N)
outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Finally, eigenvalues
# Here we accumulate eigenvalues instead of averaging
# For now we also ignore the block size
# TODO: It would be good to account for that in the future
# We also use the actual output files, not dygraph data files
# These files may include various numbers of (squared) eigenvalues
# We only consider files with at least 100, and ignore all past 100
# Finally, we use the minimum 100th eigenvalue as a maximum cutoff
files = 'eig.*'
allFiles = glob.glob(files)
cfgs = []
for filename in allFiles:
  cfg = int((filename.split('.'))[-1])  # Number after last .
  if cfg >= cut:
    cfgs.append(cfg)
cfgs.sort()

count = 0   # Check that we got correct number of maximum eigenvalues
eig = []
max_eig = []
for i in cfgs:
  toOpen = 'eig.' + str(i)
  for line in open(toOpen):
    if line.startswith('Nvecs'):
      Nvecs = int((line.split())[1])
      if Nvecs < MAX:    # Did not measure enough eigenvalues
        print "Warning: skipping %s, which measured only %d eigenvalues" \
              % (toOpen, Nvecs)
        break
    elif line.startswith('EIG'):
      temp = line.split()
      eig.append(np.sqrt(float(temp[2])));
      if int(temp[1]) == MAX:    # This is a maximum eigenvalue
        max_eig.append(np.sqrt(float(temp[2])));
        break
  count += 1

# If we don't have anything to print, then we're done
if len(eig) == 0:
  sys.exit(0)

# Normalize by volume and number of eigenvalue measurements
if len(max_eig) != count:
  print "ERROR: only have %d of %d maximum eigenvalues" \
        % (len(max_eig), count)
  sys.exit(1)
norm = count * vol

# Print out sorted eigenvalues and their sum, up to the min of MAX
eig.sort()
minmax = min(max_eig)
tot = 0
outfilename = 'results/eig'
outfile = open(outfilename, 'w')
for dat in eig:
  if dat > minmax:
    break
  tot += 2.0 / float(norm)      # Factor of two for complex conjugate pairs
  # Redundant but easy to include constant norm in output
  print >> outfile, "%.8g %.8g %.4g" % (dat, tot, norm)
outfile.close()
# ------------------------------------------------------------------
