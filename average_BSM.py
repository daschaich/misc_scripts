#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph data files, if they exist, and pbp.dat otherwise
# to construct averages and standard errors
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
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<cut> <block>"
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])

# First extract lattice volume from path
# For now assume L and Nt are both two-digit numbers
path = os.getcwd()
path = path.replace('mnt', '')    # Accommodate Barkla filesystem
temp = path.split('nt')
L = int(temp[0][-2:])    # Last two digits before 'nt'
Nt = int(temp[1][:2])    # First two digits after 'nt'
vol = L**3 * Nt
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now analyze based on what files/directories we have
if os.path.isdir('data'):
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
    MDTU = float(temp[0])
    if MDTU < cut:
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
  # ----------------------------------------------------------------

  # ----------------------------------------------------------------
  # Chiral condensate is also special
  # From MILC it needs to be halved to be normalized per continuum flavor
  count = 0
  ave = 0.          # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  pbpfile = 'data/pbp.csv'
  for line in open(pbpfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU < cut:
      continue
    elif MDTU >= begin and MDTU < (begin + block_size):
      ave += float(temp[1])
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      datList.append(ave / float(2. * count))
      begin += block_size
      count = 1                     # Next block begins with this line
      ave = float(temp[1])

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList)
  N = np.size(dat)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.)
  outfilename = 'results/pbp.dat'
  outfile = open(outfilename, 'w')
  print >> outfile, "# Normalized per continuum flavor"
  print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
  outfile.close()
  # ----------------------------------------------------------------

  # ----------------------------------------------------------------
  # Polyakov loop and spatial Wilson line have only one datum per line
  # Just look at real parts for now (not mod, arg or imaginary parts)
  for obs in ['poly_r', 'xpoly_r']:
    count = 0
    ave = 0.          # Accumulate within each block
    datList = []
    begin = cut       # Where each block begins, to be incremented
    obsfile = 'data/' + obs + '.csv'
    for line in open(obsfile):
      if line.startswith('M'):
        continue
      temp = line.split(',')
      MDTU = float(temp[0])
      if MDTU < cut:
        continue
      elif MDTU >= begin and MDTU < (begin + block_size):
        ave += float(temp[1])
        count += 1
      elif MDTU >= (begin + block_size):  # Move on to next block
        datList.append(ave / count)
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
  # ----------------------------------------------------------------

  # ----------------------------------------------------------------
  # For the plaquette and link differences, the last datum on each line
  # is the euclidean norm I'm interested in for now
  # Also need to normalize the link differences by (half) the volume
  for obs in ['plaq_diff', 'link_diff']:
    count = 0
    ave = 0.          # Accumulate within each block
    datList = []
    begin = cut       # Where each block begins, to be incremented
    obsfile = 'data/' + obs + '.csv'
    for line in open(obsfile):
      if line.startswith('M'):
        continue
      temp = line.split(',')
      MDTU = float(temp[0])
      if MDTU < cut:
        continue
      elif MDTU >= begin and MDTU < (begin + block_size):
        if obs == 'plaq_diff':
          ave += float(temp[5])
        elif obs == 'link_diff':
          ave += float(temp[5]) * 2. / float(vol)
        count += 1
      elif MDTU >= (begin + block_size):  # Move on to next bloc
        datList.append(ave / float(count))
        begin += block_size
        count = 1                     # Next block begins with this line
        if obs == 'plaq_diff':
          ave = float(temp[5])
        elif obs == 'link_diff':
          ave = float(temp[5]) * 2. / float(vol)

    # Now print mean and standard error, assuming N>1
    dat = np.array(datList)
    N = np.size(dat)
    ave = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.)
    outfilename = 'results/' + obs + '.dat'
    outfile = open(outfilename, 'w')
    print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
    outfile.close()
  # ----------------------------------------------------------------

  # ----------------------------------------------------------------
  # For the Wilson flow, the fourth datum on each line
  # is the c=0.3 running coupling I'm interested in for now
  # Add rough finite-volume perturbative correction (1+delta)~0.97
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
    MDTU = float(temp[0])
    if MDTU < cut:
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



# ------------------------------------------------------------------
# If we don't have a data/ directory,
# we should still have pbp.dat and Wflow.dat
# The FUEL pbp normalization is for four flavors, rather than the 2f above
elif os.path.isfile('pbp.dat'):
  count = 0
  ave = 0.          # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  pbpfile = 'pbp.dat'
  for line in open(pbpfile):
    temp = line.split()
    if temp[1] == 'null':   # Tag for missing measurements
      continue
    MDTU = float(temp[0])
    if MDTU < cut:
      continue
    elif MDTU >= begin and MDTU < (begin + block_size):
      ave += float(temp[1])
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      datList.append(ave / float(4. * count))
      begin += block_size
      count = 1                     # Next block begins with this line
      ave = float(temp[1])

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList)
  N = np.size(dat)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.)
  outfilename = 'results/pbp.dat'
  outfile = open(outfilename, 'w')
  print >> outfile, "# Normalized per continuum flavor"
  print >> outfile, "%.8g %.4g # %d" % (ave, err, N)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# As above, for the Wilson flow we're interested in the fourth column
# for the c=0.3 running coupling with (1+delta)~0.97
if os.path.isfile('Wflow.dat') and not os.path.isdir('data'):
  count = 0
  ave = 0.          # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  vol_corr = 0.97
  flowfile = 'Wflow.dat'
  for line in open(flowfile):
    temp = line.split()
    MDTU = float(temp[0])
    if MDTU < cut:
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
