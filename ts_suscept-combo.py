#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# Parse dygraph time-series data file to compute susceptibility,
# skewness and (excess) kurtosis
# for the plaquette, (Wilson-flowed) Polyakov loop and chiral condensate
# (Other targets may be added in the future)

# When there are multiple measurements in a block,
# we can't use the naive sum_i(dat_i - vev)^n expressions!
# Need to accumulate dat^n within each block instead

# Combine high- and low-start runs
# Save output in both runs' results directories

# Parse arguments: first specify ensemble by beta and mass
# Then give both thermalization cuts, and finally
# a combined block size (larger than combined time...)
# We discard any partial blocks at the end of each run
if len(sys.argv) < 6:
  print "Usage:", str(sys.argv[0]), "<beta> <mass>",
  print "<high cut> <low cut> <block>"
  sys.exit(1)
beta = str(sys.argv[1])
mass = str(sys.argv[2])
hi_cut = int(sys.argv[3])
lo_cut = int(sys.argv[4])
block_size = int(sys.argv[5])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Set up directories from path and input
hi_dir = 'b' + beta + '_high_m' + mass
lo_dir = 'b' + beta + '_low_m' + mass

# First make sure we're calling this from the right place
toCheck = [hi_dir + '/data', lo_dir + '/data']
for i in toCheck:
  if not os.path.isdir(i):
    print "ERROR:", i, "does not exist"
    sys.exit(1)

# Extract number of flavors for pbp normalization
# For quenched "valence pbp" we want the 4f normalization
path = os.getcwd()
if '4f' in path or '0f' in path:
  pbp_norm = 1.0
elif '8f' in path:
  pbp_norm = 0.5
else:
  print "ERROR: So far only 4f and 8f set up"
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# For plaquette, average two data per line
# For pbp we may need to normalize per continuum flavor
# For poly(_mod), have only one datum per line
# For Wpoly(_mod), use the last (fourth) number for c=0.5
for obs in ['plaq', 'pbp', 'Wpoly', 'Wpoly_mod', 'poly_r', 'poly_mod']:
  ave = 0.0         # Accumulate within each block
  aveSq = 0.0
  aveCu = 0.0
  aveFo = 0.0
  count = 0
  datList = []
  sqList = []
  cuList = []
  foList = []

  # First high-start run
  begin = hi_cut    # Where each block begins, to be incremented
  obsfile = hi_dir + '/data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= hi_cut:
      continue

    # Accumulate within block
    elif MDTU > begin and MDTU <= (begin + block_size):
      if obs == 'plaq':
        tr = 0.5 * (float(temp[1]) + float(temp[2]))
      elif obs == 'pbp':
        tr = pbp_norm * float(temp[1])
      elif obs == 'poly_r' or obs == 'poly_mod':
        tr = float(temp[1])
      elif obs == 'Wpoly' or obs == 'Wpoly_mod':
        tr = float(temp[-1])
      ave += tr
      aveSq += tr * tr
      aveCu += tr**3
      aveFo += tr**4
      count += 1

      # If that "<=" is really "==" then we are done
      # Record this block and re-initialize for the next block
      if MDTU == (begin + block_size):
        datList.append(ave / float(count))
        sqList.append(aveSq / float(count))
        cuList.append(aveCu / float(count))
        foList.append(aveFo / float(count))

        begin += block_size
        ave = 0.0
        aveSq = 0.0
        aveCu = 0.0
        aveFo = 0.0
        count = 0

    # This should never happen
    elif MDTU > (begin + block_size):
      print "ERROR: Unexpected behavior in %s, aborting" % obsfile
      sys.exit(1)

  # Clear any partial block from the end of the high-start run
  # and add the low-start run
  ave = 0.0         # Accumulate within each block
  aveSq = 0.0
  aveCu = 0.0
  aveFo = 0.0
  count = 0
  begin = lo_cut    # Where each block begins, to be incremented
  obsfile = lo_dir + '/data/' + obs + '.csv'
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= lo_cut:
      continue

    # Accumulate within block
    elif MDTU > begin and MDTU <= (begin + block_size):
      if obs == 'plaq':
        tr = 0.5 * (float(temp[1]) + float(temp[2]))
      elif obs == 'pbp':
        tr = pbp_norm * float(temp[1])
      elif obs == 'poly_r' or obs == 'poly_mod':
        tr = float(temp[1])
      elif obs == 'Wpoly' or obs == 'Wpoly_mod':
        tr = float(temp[-1])
      ave += tr
      aveSq += tr * tr
      aveCu += tr**3
      aveFo += tr**4
      count += 1

      # If that "<=" is really "==" then we are done
      # Record this block and re-initialize for the next block
      if MDTU == (begin + block_size):
        datList.append(ave / float(count))
        sqList.append(aveSq / float(count))
        cuList.append(aveCu / float(count))
        foList.append(aveFo / float(count))

        begin += block_size
        ave = 0.0
        aveSq = 0.0
        aveCu = 0.0
        aveFo = 0.0
        count = 0

    # This should never happen
    elif MDTU > (begin + block_size):
      print "ERROR: Unexpected behavior in %s, aborting" % obsfile
      sys.exit(1)

  # Require multiple blocks, N>1
  N = len(datList)
  if N < 2:
    print "ERROR: need multiple blocks to take average"
    sys.exit(1)

  # Now construct jackknife samples through single-block elimination
  # Do any (e.g., volume) normalization offline
  # Expanding (1/N) sum_i (dat_i - vev)^n, we have
  #   chi = <obs^2> - <obs>^2
  #   S   = [<obs^3> - 3<obs^2>*<obs> + 2<obs>^3] / chi^{3/2}
  #   ka  = [<obs^4> - 4<obs^3>*<obs> + 6<obs^2>*<obs>^2 - 3<obs>^4] / chi^2
  # Finally 'excess' kurtosis is ka - 3
  dat = np.array(datList, dtype = np.float64)
  sq = np.array(sqList, dtype = np.float64)
  cu = np.array(cuList, dtype = np.float64)
  fo = np.array(foList, dtype = np.float64)
  chi = np.zeros(N, dtype = np.float64)
  S   = np.zeros(N, dtype = np.float64)
  ka  = np.zeros(N, dtype = np.float64)
  for i in range(N):    # Jackknife samples
    vev = np.mean(np.delete(dat, i))
    sq_vev = np.mean(np.delete(sq, i))
    cu_vev = np.mean(np.delete(cu, i))
    fo_vev = np.mean(np.delete(fo, i))
    chi[i] = sq_vev - vev * vev
    num = cu_vev - 3.0 * sq_vev * vev + 2.0 * vev**3
    S[i] = num / np.power(chi[i], 1.5)
    num = fo_vev - 4.0 * cu_vev * vev + 6.0 * sq_vev * vev**2 - 3.0 * vev**4
    ka[i] = num / chi[i]**2 - 3.0

  # Sanity check -- compare against averages computed separately
#  print obs, "ave = %.8g" % np.mean(dat)

  # Now we can average over jackknife samples and print out results
  outfilename = hi_dir + '/results/' + obs + '.suscept-combo'
  outfile_hi = open(outfilename, 'w')
  outfilename = lo_dir + '/results/' + obs + '.suscept-combo'
  outfile_lo = open(outfilename, 'w')

  ave = np.mean(chi)
  var = (N - 1.0) * np.mean((chi - ave)**2)
  print >> outfile_hi, "suscept %.8g %.4g # %d" % (ave, np.sqrt(var), N)
  print >> outfile_lo, "suscept %.8g %.4g # %d" % (ave, np.sqrt(var), N)

  ave = np.mean(S)
  var = (N - 1.0) * np.mean((S - ave)**2)
  print >> outfile_hi, "skewness %.8g %.4g # %d" % (ave, np.sqrt(var), N)
  print >> outfile_lo, "skewness %.8g %.4g # %d" % (ave, np.sqrt(var), N)

  ave = np.mean(ka)
  var = (N - 1.0) * np.mean((ka - ave)**2)
  print >> outfile_hi, "kurtosis %.8g %.4g # %d" % (ave, np.sqrt(var), N)
  print >> outfile_lo, "kurtosis %.8g %.4g # %d" % (ave, np.sqrt(var), N)
  outfile_hi.close()
  outfile_lo.close()
# ------------------------------------------------------------------
