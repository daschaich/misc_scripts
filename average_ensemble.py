#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph data files to construct averages and standard errors
# with given thermalization cut and block size

# We consider only one ensemble at a time,
# appending to results files, rather than overwriting them
# These files record results as functions of the mass and of the coupling

# Also accumulate eigenvalues from Out/eig_* for rho(lambda)
# TODO: Account for block size in eigenvalue accumulation
MAX = 99    # Hard code number of eigenvalues to extract
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: first is ensemble tag, second is thermalization cut,
# third is block size (should be larger than auto-correlation time)
# We discard any partial blocks at the end
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<tag> <cut> <bin>"
  sys.exit(1)
tag = str(sys.argv[1])
cut = int(sys.argv[2])
block_size = int(sys.argv[3])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)
if not 'Run_' in os.getcwd():
  print "ERROR: need to call from a Run_* directory"
  sys.exit(1)

# Extract lattice volume from first output file
firstFile = 'Out/out_' + tag + '.1'
for line in open(firstFile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  elif line.startswith('iseed '):
    break   # Done scanning through file
vol = L**3 * Nt

# Set the maximum number of RG blockings
blMax = 0
if L <= Nt:
  temp = L
else:
  temp = Nt
while temp % 2 == 0 and temp > 2:
  temp /= 2
  blMax += 1

# Use tag to set beta, mass and start
temp = tag.split('_')
start = temp[0]
#vol = temp[1]    # Would overwrite int above
beta = float(temp[2])
ratio = temp[3]
mass = float(temp[4])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Translate MDTU cut into corresponding output file, for eigenvalues
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
  print tag, cut, block_size
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Plaquette is special -- two data (to be averaged) per line
count = 0
ave = 0.          # Accumulate within each block
datList = []
begin = cut       # Where each block begins, to be incremented
plaqfile = 'data/plaq.' + tag + '.csv'
for line in open(plaqfile):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU < cut:
    continue
  elif MDTU >= begin and MDTU < (begin + block_size):
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
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1)
outfilename = 'results/plaq_m' + str(mass) + '.' + start
outfile = open(outfilename, 'a')
print >> outfile, "%.8g %.8g %.4g # %d" % (beta, ave, err, N)
outfile.close()
outfilename = 'results/plaq_b' + str(beta) + '.' + start
outfile = open(outfilename, 'a')
print >> outfile, "%.8g %.8g %.4g # %d" % (mass, ave, err, N)
outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Chiral condensate is also special
# It needs to be halved to be normalized per continuum flavor
count = 0
ave = 0.          # Accumulate within each block
datList = []
begin = cut       # Where each block begins, to be incremented
pbpfile = 'data/pbp.' + tag + '.csv'
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
err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1.0)

massfilename = 'results/pbp_m' + str(mass) + '.' + start
if not os.path.isfile(massfilename):
  massfile = open(massfilename, 'a')
  print >> massfile, "# Normalized per continuum flavor"
else:
  massfile = open(massfilename, 'a')

betafilename = 'results/pbp_b' + str(beta) + '.' + start
if not os.path.isfile(betafilename):
  betafile = open(betafilename, 'a')
  print >> betafile, "# Normalized per continuum flavor"
else:
  betafile = open(betafilename, 'a')

print >> massfile, "%.8g %.8g %.4g # %d" % (beta, ave, err, N)
print >> betafile, "%.8g %.8g %.4g # %d" % (mass, ave, err, N)
massfile.close()
betafile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Polyakov loop and spatial Wilson line have only one datum per line
# Just look at real parts for now (not mod, arg or imaginary parts)
# At the moment we don't have xpoly_r without RG-blocked measurements
mcrgFiles = 'Out/mcrg_' + tag + '.*'
mcrgCheck = glob.glob(mcrgFiles)
if len(mcrgCheck) == 0:
  polyObs = ['poly_r']
else:
  polyObs = ['poly_r']
for obs in polyObs:
  count = 0
  ave = 0.          # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.' + tag + '.csv'
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
      datList.append(ave / float(count))
      begin += block_size
      count = 1                     # Next block begins with this line
      ave = float(temp[1])

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList)
  N = np.size(dat)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1)
  outfilename = 'results/' + obs + '_m' + str(mass) + '.' + start
  outfile = open(outfilename, 'a')
  print >> outfile, "%.8g %.8g %.4g # %d" % (beta, ave, err, N)
  outfile.close()
  outfilename = 'results/' + obs + '_b' + str(beta) + '.' + start
  outfile = open(outfilename, 'a')
  print >> outfile, "%.8g %.8g %.4g # %d" % (mass, ave, err, N)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Blocked data files have blMax + 1 data per line (including unblocked)
# We already got the "+ 1" above, but include it here as well, to check
# Just look at real parts for now (not mod, arg or imaginary parts)
if len(mcrgCheck) == 0:
  mcrgObs = []
else:
  mcrgObs = ['plaqB', 'poly_rB', 'xpoly_rB']
for obs in mcrgObs:
  count = 0
  ave = [0. for x in range(blMax + 1)] # Accumulate within each block
  datList = [[] for x in range(blMax + 1)]
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.' + tag + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU < cut:
      continue
    elif MDTU >= begin and MDTU < (begin + block_size):
      for i in range(0, blMax + 1):
        ave[i] += float(temp[i + 1])
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      for i in range(0, blMax + 1):
        datList[i].append(ave[i] / float(count))
        ave[i] = float(temp[i + 1])
      begin += block_size
      count = 1                     # Next block begins with this line

  # Now print mean and standard error, assuming N>1
  massfilename = 'results/' + obs + '_m' + str(mass) + '.' + start
  massfile = open(massfilename, 'a')
  betafilename = 'results/' + obs + '_b' + str(beta) + '.' + start
  betafile = open(betafilename, 'a')
  print >> massfile, "%.8g" % beta,
  print >> betafile, "%.8g" % mass,
  for i in range(0, blMax + 1):
    dat = np.array(datList[i])
    N = np.size(dat)
    ave = np.mean(dat, dtype = np.float64)
    err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1)
    print >> massfile, "%.8g %.4g" % (ave, err),
    print >> betafile, "%.8g %.4g" % (ave, err),
  print >> massfile, "%d" % N
  print >> betafile, "%d" % N
  massfile.close()
  betafile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the plaquette and link differences, the last datum on each line
# is the euclidean norm I'm interested in for now
# Also need to normalize the link differences by (half) the volume
for obs in ['plaq_diff', 'link_diff']:
  count = 0
  ave = 0.          # Accumulate within each block
  datList = []
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.' + tag + '.csv'
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
    elif MDTU >= (begin + block_size):  # Move on to next block
      datList.append(ave / float(count))
      begin += block_size
      count = 1                     # Next block begins with this line
      if obs == 'plaq_diff':
        ave = float(temp[5])
      elif obs == 'link_diff':
        ave = float(temp[5]) * 2.0 / float(vol)

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList)
  N = np.size(dat)
  ave = np.mean(dat, dtype = np.float64)
  err = np.std(dat, dtype = np.float64) / np.sqrt(N - 1)
  outfilename = 'results/' + obs + '_m' + str(mass) + '.' + start
  outfile = open(outfilename, 'a')
  print >> outfile, "%.8g %.8g %.4g # %d" % (beta, ave, err, N)
  outfile.close()
  outfilename = 'results/' + obs + '_b' + str(beta) + '.' + start
  outfile = open(outfilename, 'a')
  print >> outfile, "%.8g %.8g %.4g # %d" % (mass, ave, err, N)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Wilson-flowed Polyakov loop, considering both c=0.2 and 0.3 for now
# We also only care about these for fixed mass, not fixed beta
WpolyFiles = 'Out/Wpoly_' + tag + '.*'
WpolyCheck = glob.glob(WpolyFiles)
if len(WpolyCheck) == 0:
  WpolyObs = []
else:
  WpolyObs = ['Wpoly']
for obs in WpolyObs:
  count = 0
  ave2 = 0.0        # Accumulate within each block
  ave3 = 0.0
  datList2 = []
  datList3 = []
  begin = cut       # Where each block begins, to be incremented
  obsfile = 'data/' + obs + '.' + tag + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU < cut:
      continue
    elif MDTU >= begin and MDTU < (begin + block_size):
      ave2 += float(temp[1])
      ave3 += float(temp[2])
      count += 1
    elif MDTU >= (begin + block_size):  # Move on to next block
      datList2.append(ave2 / float(count))
      datList3.append(ave3 / float(count))
      begin += block_size
      count = 1                     # Next block begins with this line
      ave = float(temp[1])
      ave2 = float(temp[1])
      ave3 = float(temp[2])

  # Now print mean and standard error, assuming N>1
  dat2 = np.array(datList2)
  dat3 = np.array(datList3)
  N = np.size(dat2)
  ave2 = np.mean(dat2, dtype = np.float64)
  err2 = np.std(dat2, dtype = np.float64) / np.sqrt(N - 1.)
  ave3 = np.mean(dat3, dtype = np.float64)
  err3 = np.std(dat3, dtype = np.float64) / np.sqrt(N - 1.)
  outfilename = 'results/' + obs + '_m' + str(mass) + '.' + start
  outfile = open(outfilename, 'a')
  print >> outfile, "%.8g %.8g %.4g %.8g %.4g # %d" \
                    % (beta, ave2, err2, ave3, err3, N)
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
print "Seem to be missing some eigenvalues?..."
files = 'Out/eig_*_' + tag + '.*'
allFiles = glob.glob(files)
cfgs = []
for filename in allFiles:
  cfg = int((filename.split('.'))[-1])  # Number after last .
  if cfg >= file_cut:
    cfgs.append(cfg)
cfgs.sort()

# If multiple files, take the one with the largest number of eigenvalues
# Supposing that the extra iterations may have helped refine the first MAX
count = 0   # Check that we got correct number of maximum eigenvalues
eig = []
max_eig = []
for i in cfgs:
  files = 'Out/eig_*_' + tag + '.' + str(i)
  allfiles = glob.glob(files)
  Nvecs = int((allfiles[0].split('_'))[1])
  if len(allfiles) > 1:
    for filename in allfiles:
      temp = int((filename.split('_'))[1])
      if temp > Nvecs:
        Nvecs = temp
  if Nvecs <= MAX:
    continue

  toOpen = 'Out/eig_' + str(Nvecs) + '_' + tag + '.' + str(i)
  for line in open(toOpen):
    if line.startswith('EIG'):
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
outfilename = 'results/eig.' + tag
outfile = open(outfilename, 'w')
for dat in eig:
  if dat > minmax:
    break
  tot += 2.0 / float(norm)      # Factor of two for complex conjugate pairs
  # Redundant but easy to include constant norm in output
  print >> outfile, "%.8g %.8g %.4g" % (dat, tot, norm)
outfile.close()
# ------------------------------------------------------------------
