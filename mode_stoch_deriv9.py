#!/usr/bin/python
import glob
import os
import sys
import time
import numpy as np
# ------------------------------------------------------------------
# Run jackknife calculations of <nu'> and <nu''>,
# to find lambda <nu''> / <nu'> = alpha and gamma = 4 / (alpha + 1) - 1
# Fix nine-point first and second derivatives

# Only read in averages over stochastic sources
# (i.e., one datum per configuration per lambda)
# ignoring the stochastic errors since we're doing a jackknife analysis

# Require that all configurations being analyzed consider the same lambda
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: all files to analyze
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<files>"
  sys.exit(1)

files = []
for i in range(1, len(sys.argv)):
  files.append(str(sys.argv[i]))
runtime = -time.time()
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct list of which configurations have been analyzed
cfgs = []
for filename in files:
  cfg = int((filename.split('.'))[-1])  # Number after last .
  if cfg not in cfgs:
    cfgs.append(cfg)
cfgs.sort()

if len(cfgs) == 0:
  print "ERROR: specified files not found"
  sys.exit(1)

# Now use first configuration (including _fine, _more, etc. files)
# to construct list of la at which we have nu
first = '.' + str(cfgs[0])
all_lambda = []
for filename in files:
  if filename.endswith(first):
    for line in open(filename):
      if line.startswith('nu('):    # Want la from 'nu(la): # ( ave over 5)'
        temp = line.split('): ')
        la = ((temp[0]).split('('))[1]    # Save as string
        if la in all_lambda:
          print "ERROR: lambda=%.4g measured twice" % float(la)
          sys.exit(1)
        else:
          all_lambda.append(la)
all_lambda.sort()
N_la = len(all_lambda)

if N_la == 0:
  print "ERROR: no data found"
  sys.exit(1)

# Now we can construct arrays of all measurements of each nu(la)
dat = [[] for x in range(N_la)]
datSq = [[] for x in range(N_la)]
for filename in files:
  for line in open(filename):
    if line.startswith('nu('):    # Want # from 'nu(la): # ( ave over 5)'
      temp = line.split('): ')
      # Figure out which element of dat this goes to
      la = ((temp[0]).split('('))[1]    # Save as string
      toAdd = all_lambda.index(la)
      temp = float(((temp[1]).split())[0])
      dat[toAdd].append(temp)
      datSq[toAdd].append(temp**2)

# Check that all configurations being analyzed consider the same lambda
# Also save lambdas as floats in addition to strings
Nmeas = len(dat[0])
la_dat = np.empty(N_la, dtype = np.float)
for i in range(N_la):
  la_dat[i] = float(all_lambda[i])
  if not len(dat[i]) == Nmeas:
    print "ERROR: lambda=%.4g only measured on" % la_dat[i],
    print "%d of %d configs" % (len(dat[i]), Nmeas)
    sys.exit(1)

# !!! We assume that all lambda are evenly spaced
# Because floating-point comparisons are too annoying for the moment
epsilon = la_dat[1] - la_dat[0]
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we are ready to construct jackknife samples through single elimination
tot = np.array([sum(dat[i]) for i in range(N_la)])

# We can't construct nine-point derivatives at the first and last four points
Npts = N_la - 8
first = np.empty((Npts, Nmeas), dtype = np.float)
secon = np.empty((Npts, Nmeas), dtype = np.float)
ratio = np.empty((Npts, Nmeas), dtype = np.float)
toAve = np.empty((N_la, Nmeas), dtype = np.float)
for i in range(Nmeas):        # Jackknife samples
  nu = np.empty(N_la, dtype = np.float)
  for la in range(N_la):
    toAve[la][i] = (tot[la] - dat[la][i]) / (Nmeas - 1.)
    nu[la] = (tot[la] - dat[la][i]) / (Nmeas - 1.)

  # Set up first and second nine-point derivatives
  # Ignore common factor of 5040epsilon that will cancel out in the ratio
  # Left with only relative factor of epsilon in the second derivative
  # Errors will come from jackknife estimates below
  for la in range(4, N_la - 4):
    j = la - 4
    first[j][i] =   18. * (nu[la - 4] - nu[la + 4]) \
                -  192. * (nu[la - 3] - nu[la + 3]) \
                + 1008. * (nu[la - 2] - nu[la + 2]) \
                - 4032. * (nu[la - 1] - nu[la + 1])
    secon[j][i] = -14350. * nu[la] \
                +   8064. * (nu[la - 1] + nu[la + 1]) \
                -   1008. * (nu[la - 2] + nu[la + 2]) \
                +    128. * (nu[la - 3] + nu[la + 3]) \
                -      9. * (nu[la - 4] + nu[la + 4])
    secon[j][i] /= epsilon
    ratio[j][i] = la_dat[la] * secon[j][i] / first[j][i]

# Save averages for test below
nuAve = np.empty(N_la, dtype = np.float)
nuErr = np.empty(N_la, dtype = np.float)
for la in range(N_la):
  nuAve[la] = np.average(toAve[la])
  var = (Nmeas - 1.) * np.sum((toAve[la] - nuAve[la])**2) / float(Nmeas)
  nuErr[la] = np.sqrt(var)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples, construct the ratio
# lambda <nu''> / <nu'> = alpha
# and print out results with simple/simplistic error propagation
print "# Jackknifed results from", Nmeas, "measurements",
print "(listed at bottom)"
print "# Nine-point first and second derivatives"
print "# lambda gamma delta_gamma"

# We lose the first and last four points by taking nine-point derivatives
for la in range(Npts):
  x = la_dat[la + 4]
  print "%.4g" % x,
  # First try jackknifing the ratio
  alpha = np.average(ratio[la])
  var = (Nmeas - 1.) * np.sum((ratio[la] - alpha)**2) / float(Nmeas)
  alpha_err = np.sqrt(var)
  gamma = 4. / (alpha + 1.) - 1.
  err = 4. * alpha_err / (alpha + 1.)**2
  print "%.6g %.4g" % (gamma, err),

  # Now try taking the ratio of the jackknifed derivatives
  num = np.average(secon[la])
  den = np.average(first[la])
  var = (Nmeas - 1.) * np.sum((secon[la] - num)**2) / float(Nmeas)
  num_err = np.sqrt(var)
  var = (Nmeas - 1.) * np.sum((first[la] - den)**2) / float(Nmeas)
  den_err = np.sqrt(var)
  alpha = x * num / den
#  print "%.4g %.4g %.4g %.4g #" % (x, num, den, alpha),
  alpha_err = np.abs(alpha) * np.sqrt((num_err / num)**2 + (den_err / den)**2)
  gamma = 4. / (alpha + 1.) - 1.
  err = 4. * alpha_err / (alpha + 1.)**2
  print "%.6g %.4g" % (gamma, err),

  # Finally, re-calculate derivatives from jackknifed data
  # la runs over Npts = N_la + 4: shift all indices up by four
  num = -14350. * nuAve[la + 4] \
      +   8064. * (nuAve[la + 3] + nuAve[la + 5]) \
      -   1008. * (nuAve[la + 2] + nuAve[la + 6]) \
      +    128. * (nuAve[la + 1] + nuAve[la + 7]) \
      -      9. * (nuAve[la] + nuAve[la + 8])
  num /= epsilon
  var = (14350. * nuErr[la + 4])**2 \
      +  (8064. * nuErr[la + 3])**2 + (8064. * nuErr[la + 5])**2 \
      +  (1008. * nuErr[la + 2])**2 + (1008. * nuErr[la + 6])**2 \
      +   (128. * nuErr[la + 1])**2 + (128. * nuErr[la + 7])**2 \
      +     (9. * nuErr[la])**2 + (9. * nuErr[la + 8])**2
  num_err = np.sqrt(var) / epsilon

  den =   18. * (nuAve[la] - nuAve[la + 8]) \
      -  192. * (nuAve[la + 1] - nuAve[la + 7]) \
      + 1008. * (nuAve[la + 2] - nuAve[la + 6]) \
      - 4032. * (nuAve[la + 3] - nuAve[la + 5])
  var =   (18. * nuErr[la])**2 + (18. * nuErr[la + 8])**2 \
      +  (192. * nuErr[la + 1])**2 + (192. * nuErr[la + 7])**2 \
      + (1008. * nuErr[la + 2])**2 + (1008. * nuErr[la + 6])**2 \
      + (4032. * nuErr[la + 3])**2 + (4032. * nuErr[la + 5])**2
  den_err = np.sqrt(var)

  alpha = x * num / den
  alpha_err = np.abs(alpha) * np.sqrt((num_err / num)**2 + (den_err / den)**2)
  gamma = 4. / (alpha + 1.) - 1.
  err = 4. * alpha_err / (alpha + 1.)**2
  print "%.6g %.4g" % (gamma, err)

runtime += time.time()
print "# Runtime: %.2g seconds" % runtime

# Print list of files
print "# Files analyzed:"
for filename in files:
  print "#", filename
# ------------------------------------------------------------------
