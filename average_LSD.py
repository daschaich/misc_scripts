#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import emcee.autocorr as acor
# ------------------------------------------------------------------
# Parse dygraph data files to construct averages and standard errors
# with given thermalization cut and block size

# Assume one ensemble per directory (overwrite results files)

# Also accumulate eigenvalues from eig.* for rho(lambda)
# TODO: Account for block size in eigenvalue accumulation
MAX = 99    # Hard code number of eigenvalues to extract
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: first is thermalization cut, second is block size
# Will check block size is larger than Wpoly_mod auto-correlation time
# We discard any partial blocks at the end
# I also need to specify where I took over the run
if len(sys.argv) < 4:
  print("Usage:", str(sys.argv[0]), "<cut> <block> <start>")
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
start = int(sys.argv[3])

# First extract lattice volume from path
# For now handle only one- or two-digit L and Nt
path = os.getcwd()
path = path.replace('mnt', '')    # Accommodate Barkla filesystem
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

# Set pbp normalization depending on number of flavors
# For quenched "valence pbp" we want the 4f normalization
if '4f' in path or '0f' in path:
  pbp_norm = 1.0
elif '8f' in path:
  pbp_norm = 0.5
else:
  print("ERROR: So far only 0f, 4f and 8f set up")
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
  print("Error: no data to analyze ", end='')
  print("since cut=%d but we only have %d MDTU" % (cut, final_MDTU))
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Check that block size is larger than Wpoly_mod auto-correlation time
# (Format: MDTU,c=0.2,0.3,0.4,0.5)
dat = []
sep = 10
prev = start
for line in open('data/Wpoly_mod.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = int(temp[0])

  # Need to check separation and update prev before skipping to therm cut
  if not MDTU - prev == sep:
    print("Error: Wpoly_mod meas at %d and %d not separated by %d" \
          % (prev, MDTU, sep))
    sys.exit(1)
  prev = MDTU

  if MDTU <= cut:
    continue
  dat.append(float(temp[4]))

# Arguments explained in emcee.readthedocs.io/en/stable/user/autocorr/
#                    and emcee.readthedocs.io/en/stable/tutorials/autocorr/
# 'c' is step size for window search (default 5)
#     Larger c should reduce noise, but can add bias...
# 'tol' is min ratio between data and tau (default 50)
# 'Quiet' prints warning rather than shutting down if tol not satisfied
tau = acor.integrated_time(np.array(dat), c=5, tol=10, quiet=True)
tau *= sep
if tau > block_size:
  print("Error: Wpoly_mod autocorrelation time %d " % tau, end='')
  print("is larger than block size %d " % block_size, end='')
  print("in %s" % path)
  sys.exit(1)

# Record Wpoly_mod auto-correlation time for future reference
# Include effective number of independent measurements
eff_stat = np.floor(len(dat) * sep / tau)
outfilename = 'results/Wpoly_mod.autocorr'
outfile = open(outfilename, 'w')
print("%d # %d" % (tau, eff_stat), file=outfile)
outfile.close()

# Also keep track of the pbp autocorrelation time
# Print warnings if larger than the Wpoly_mod tau computed above
# Format: MDTU,Tr(X1)^2,...,Tr(X9)^2
dat = []
sep = 1       # Measured after each trajectory
prev = 0
for line in open('data/pbp.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = int(temp[0])

  # Need to check separation and update prev before skipping to therm cut
  if not MDTU - prev == sep:
    print("Error: pbp meas at %d and %d not separated by %d" \
          % (prev, MDTU, sep))
    sys.exit(1)
  prev = MDTU

  if MDTU <= cut:
    continue
  dat.append(float(temp[1]))

# Arguments discussed above --- make this advice rather than requirement
tau_pbp = acor.integrated_time(np.array(dat), c=5, tol=10, quiet=True)
tau_pbp *= sep
if tau_pbp > block_size:
  print("Warning: pbp autocorrelation time %d " % tau_pbp, end='')
  print("is larger than block size %d " % block_size, end='')
  print("in %s" % path)

# Record pbp auto-correlation time for future reference
# Include effective number of independent measurements
eff_stat = np.floor(len(dat) * sep / tau_pbp)
outfilename = 'results/pbp.autocorr'
outfile = open(outfilename, 'w')
print("%d # %d" % (tau_pbp, eff_stat), file=outfile)
if tau_pbp > tau:
  print("Warning: %d for pbp " % tau_pbp, end='', file=outfile)
  print("larger than %d for Wpoly_mod" % tau, file=outfile)
outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For plaquette, average two data per line
# For pbp we may need to normalize per continuum flavor
# For (x)poly(_mod), have only one datum per line
# Currently not looking at arg or imaginary parts
for obs in ['plaq', 'pbp', 'poly_r', 'poly_mod', 'xpoly_r', 'xpoly_mod']:
  ave = 0.0         # Accumulate within each block
  count = 0
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

    # Accumulate within block
    elif MDTU > begin and MDTU <= (begin + block_size):
      if obs == 'plaq':
        tr = 0.5 * (float(temp[1]) + float(temp[2]))
      elif obs == 'pbp':
        tr = pbp_norm * float(temp[1])
      elif obs == 'poly_r' or obs == 'poly_mod' \
        or obs == 'xpoly_r' or obs == 'xpoly_mod':
        tr = float(temp[1])
      ave += tr
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU == (begin + block_size):
        datList.append(ave / float(count))
        begin += block_size
        ave = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList, dtype = np.float64)
  N = np.size(dat)
  ave = np.mean(dat)
  err = np.std(dat) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  if obs == 'pbp':
    print("# Normalized per continuum flavor", file=outfile)
  print("%.8g %.4g # %d" % (ave, err, N), file=outfile)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For algorithmic/cost quantities
# we're interested in the first datum on each line
# Need to work in terms of trajectories rather than MDTU
for obs in ['wallTU', 'cg_iters', 'accP', 'exp_dS']:
  if '0f' in path and obs == 'cg_iters':
    continue        # No CG iterations for pure-gauge runs

  ave = 0.0         # Accumulate within each block
  count = 0
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

    # Accumulate within block
    elif traj > begin and traj <= (begin + block_size):
      ave += float(temp[1])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if traj == (begin + block_size):
        datList.append(ave / float(count))
        begin += block_size
        ave = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif traj > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList, dtype = np.float64)
  N = np.size(dat)
  ave = np.mean(dat)
  err = np.std(dat) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print("%.8g %.4g # %d" % (ave, err, N), file=outfile)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Blocked data files have blMax + 1 data per line (including unblocked)
# We already got the "+ 1" above, but include it here as well, to check
# Just look at real parts for now (not mod, arg or imaginary parts)
for obs in ['plaqB', 'poly_rB', 'xpoly_rB']:
  ave = [0.0 for x in range(blMax + 1)] # Accumulate within each block
  count = 0
  datList = [[] for x in range(blMax + 1)]
  begin = cut       # Where each block begins, to be incremented

  # Only run if we have blocked data
  obsfile = 'data/' + obs + '.csv'
  if not os.path.isfile(obsfile):
    continue
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0]) + start
    if MDTU <= cut:
      continue

    # Accumulate within block
    elif MDTU > begin and MDTU <= (begin + block_size):
      for i in range(blMax + 1):
        ave[i] += float(temp[i + 1])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU == (begin + block_size):
        for i in range(blMax + 1):
          datList[i].append(ave[i] / float(count))
          ave[i] = 0.0
        begin += block_size
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  outfilename = 'results/' + obs +  '.dat'
  outfile = open(outfilename, 'w')
  for i in range(blMax + 1):
    dat = np.array(datList[i], dtype = np.float64)
    N = np.size(dat)
    ave = np.mean(dat)
    err = np.std(dat) / np.sqrt(N - 1.0)
    print("%.8g %.4g " % (ave, err), end='', file=outfile)
  print("# %d" % N, file=outfile)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the plaquette and link differences, the last datum on each line
# is the euclidean norm I'm interested in for now
# Also need to normalize the link differences by (half) the volume
for obs in ['plaq_diff', 'link_diff']:
  ave = 0.0         # Accumulate within each block
  count = 0
  datList = []
  begin = cut       # Where each block begins, to be incremented

  # Only run if we have S4b data
  obsfile = 'data/' + obs + '.csv'
  if not os.path.isfile(obsfile):
    continue
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0]) + start
    if MDTU <= cut:
      continue

    # Accumulate within block
    elif MDTU > begin and MDTU <= (begin + block_size):
      if obs == 'plaq_diff':
        tr = float(temp[5])
      elif obs == 'link_diff':
        tr = float(temp[5]) * 2.0 / float(vol)
      ave += tr
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU == (begin + block_size):
        datList.append(ave / float(count))
        begin += block_size
        ave = 0.0
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  dat = np.array(datList, dtype = np.float64)
  N = np.size(dat)
  ave = np.mean(dat)
  err = np.std(dat) / np.sqrt(N - 1.0)
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print("%.8g %.4g # %d" % (ave, err, N), file=outfile)
  outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the Wilson flow, the fourth datum on each line
# is the c=0.3 running coupling I'm interested in for now
# Add rough finite-volume perturbative correction (1+delta)~0.97
ave = 0.0         # Accumulate within each block
count = 0
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

  # Accumulate within block
  elif MDTU > begin and MDTU <= (begin + block_size):
    ave += float(temp[3]) / vol_corr
    count += 1

    # If that "<=" is really "==" then we are done with this block
    # Record it and re-initialize for the next block
    if MDTU == (begin + block_size):
      datList.append(ave / float(count))
      begin += block_size
      ave = 0.0
      count = 0

  # This doesn't happen for ensembles I generate
  # May need to be revisited for more general applicability
  elif MDTU > (begin + block_size):
    print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
    sys.exit(1)

# Now print mean and standard error, assuming N>1
dat = np.array(datList, dtype = np.float64)
N = np.size(dat)
ave = np.mean(dat)
err = np.std(dat) / np.sqrt(N - 1.0)
outfilename = 'results/Wflow_gSq.dat'
outfile = open(outfilename, 'w')
print("%.8g %.4g # %d" % (ave, err, N), file=outfile)
outfile.close()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For the Wilson-flowed Polyakov loop (real part and modulus)
# and anisotropy E_ss / E_st,
# let's print out all four of c=0.2, 0.3, 0.4 and 0.5
for obs in ['Wpoly', 'Wpoly_mod', 'Wflow_aniso']:
  ave = [0.0 for x in range(4)]     # Accumulate within each block
  count = 0
  datList = [[] for x in range(4)]
  begin = cut       # Where each block begins, to be incremented
  flowfile = 'data/' + obs + '.csv'
  for line in open(flowfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0]) + start
    if MDTU <= cut:
      continue

    # Accumulate within block
    elif MDTU > begin and MDTU <= (begin + block_size):
      for i in range(4):
        ave[i] += float(temp[i + 1])
      count += 1

      # If that "<=" is really "==" then we are done with this block
      # Record it and re-initialize for the next block
      if MDTU == (begin + block_size):
        for i in range(4):
          datList[i].append(ave[i] / float(count))
          ave[i] = 0.0
        begin += block_size
        count = 0

    # This doesn't happen for ensembles I generate
    # May need to be revisited for more general applicability
    elif MDTU > (begin + block_size):
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Now print mean and standard error, assuming N>1
  outfilename = 'results/' + obs + '.dat'
  outfile = open(outfilename, 'w')
  print("# c=0.2 err c=0.3 err c=0.4 err c=0.5 err # Nblocks", file=outfile)
  for i in range(4):
    dat = np.array(datList[i], dtype = np.float64)
    N = np.size(dat)
    ave = np.mean(dat)
    err = np.std(dat) / np.sqrt(N - 1.0)
    print("%.8g %.4g " % (ave, err), end='', file=outfile)
  print("# %d" % N, file=outfile)
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
        print("Warning: skipping %s, which measured only %d eigenvalues" \
              % (toOpen, Nvecs))
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
  print("ERROR: only have %d of %d maximum eigenvalues" \
        % (len(max_eig), count))
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
  print("%.8g %.8g %.4g" % (dat, tot, norm), file=outfile)
outfile.close()
# ------------------------------------------------------------------
