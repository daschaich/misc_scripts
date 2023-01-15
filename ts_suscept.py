#!/usr/bin/python3
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

# Parse arguments: first is thermalization cut,
# second is block size (should be larger than autocorrelation time)
# We discard any partial blocks at the end
if len(sys.argv) < 3:
  print("Usage:", str(sys.argv[0]), "<cut> <block>")
  sys.exit(1)
cut = int(sys.argv[1])
block_size = int(sys.argv[2])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print("ERROR: data/ does not exist")
  sys.exit(1)

# Check that we actually have data to average
MDTUfile = 'data/TU.csv'
good = -1
for line in open(MDTUfile):
  if line.startswith('t'):
    continue
  temp = line.split(',')
  if float(temp[1]) > cut:
    good = 1
    break

if good == -1:
  print("Error: no data to analyze ", end='')
  print("since cut=%d but we only have %d MDTU" % (cut, float(temp[1])))
  sys.exit(1)

# Extract number of flavors for pbp normalization
# For quenched "valence pbp" we want the 4f normalization
path = os.getcwd()
if '4f' in path or '0f' in path:
  pbp_norm = 1.0
elif '8f' in path:
  pbp_norm = 0.5
else:
  print("ERROR: So far only 0f, 4f and 8f set up")
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
  begin = cut       # Where each block begins, to be incremented

  # Only run if we have these data (not all saved by LargeN-YM code)
  obsfile = 'data/' + obs + '.csv'
  if not os.path.isfile(obsfile):
    continue
  for line in open(obsfile):
    if line.startswith('M'):
      continue
    temp = line.split(',')
    MDTU = float(temp[0])
    if MDTU <= cut:
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
      print("ERROR: Unexpected behavior in %s, aborting" % obsfile)
      sys.exit(1)

  # Require multiple blocks, N>1
  N = len(datList)
  if N < 2:
    print("ERROR: need multiple blocks to take average")
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
#  print(obs, "ave = %.8g" % np.mean(dat))

  # Now we can average over jackknife samples and print out results
  outfilename = 'results/' + obs + '.suscept'
  outfile = open(outfilename, 'w')

  ave = np.mean(chi)
  var = (N - 1.0) * np.mean((chi - ave)**2)
  print("suscept %.8g %.4g # %d" % (ave, np.sqrt(var), N), file=outfile)

  ave = np.mean(S)
  var = (N - 1.0) * np.mean((S - ave)**2)
  print("skewness %.8g %.4g # %d" % (ave, np.sqrt(var), N), file=outfile)

  ave = np.mean(ka)
  var = (N - 1.0) * np.mean((ka - ave)**2)
  print("kurtosis %.8g %.4g # %d" % (ave, np.sqrt(var), N), file=outfile)
  outfile.close()
# ------------------------------------------------------------------
