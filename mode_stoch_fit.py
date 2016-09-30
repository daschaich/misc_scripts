#!/usr/bin/python
import os
import sys
import glob
import time
import numpy as np
from scipy import optimize
# ------------------------------------------------------------------
# Run jackknifed power-law fits of nu vs. la

# Read in number of points to use in each fit, and run a fit for each
# point until we run out of data ("sliding window" fits)
# Print out midpoint and width of each fit range for horizontal error bars

# Only read in averages over stochastic sources
# (i.e., one datum per configuration per lambda)
# ignoring the stochastic errors since we're doing a jackknife analysis

# Require that all configurations being analyzed consider the same lambda
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: start with all files to analyze,
# then give number of points per fit
if len(sys.argv) < 4:
  print "Usage:", str(sys.argv[0]), "<files> <Npts>"
  sys.exit(1)

files = []
for i in range(1, len(sys.argv) - 1):
  files.append(str(sys.argv[i]))
Npts = int(sys.argv[-1])
runtime = -time.time()

# errfunc will be minimized via least-squares optimization
fitfunc = lambda p, x: p[0] * np.power(x, p[1]) + p[2]
errfunc = lambda p, x, y, err: (fitfunc(p, x) - y) / err
p_in = [1., 1., 1.]   # Order-of-magnitude initial guesses
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
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we are ready to construct jackknife samples through single elimination
tot = np.array([sum(dat[i]) for i in range(N_la)])
totSq = np.array([sum(datSq[i]) for i in range(N_la)])

# Fit results for all jk samples and all fits each with Npts points
Nfits = N_la - Npts + 1
jkgamma = np.empty((Nfits, Nmeas), dtype = np.float)
for i in range(Nmeas):        # Jackknife samples
  nu = np.empty(N_la, dtype = np.float)
  nu_err = np.empty(N_la, dtype = np.float)
  for la in range(N_la):
    nu[la] = (tot[la] - dat[la][i]) / (Nmeas - 1.)
    nu_err[la] = (totSq[la] - datSq[la][i]) / (Nmeas - 1.)
    nu_err[la] = np.sqrt(nu_err[la] - nu[la]**2)

  # All fits for this jackknife sample
  # The least-squares fit returns alpha+1 as p_out[1]
  # The weight will be squared in the fit...
  # Ignore the errors (covariance matrix) from the fit
  for la in range(Nfits):
    X = np.array([la_dat[j] for j in range(la, la + Npts)])
    to_fit = np.array([nu[j] for j in range(la, la + Npts)])
    to_err = np.array([nu_err[j] for j in range(la, la + Npts)])
    p_out, success = optimize.leastsq(errfunc, p_in[:],
                                      args=(X, to_fit, to_err))
    jkgamma[la][i] = 4. / p_out[1] - 1.
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can average over jackknife samples and print out results
# Output format x y de_x de_y for gnuplot
print "# Jackknifed results from", Nmeas, "measurements",
print "(listed at bottom)"
print "#", Npts, "points per fit"
print "# lambda gamma delta_lambda delta_gamma"

for la in range(Nfits):
  # We lost the first and last all_lambda by taking derivatives
  # Increment both indices by 1 to compensate
  de_x = (float(all_lambda[la + Npts]) - float(all_lambda[la + 1])) / 2
  x = (float(all_lambda[la + Npts]) + float(all_lambda[la + 1])) / 2

  ga_ave = np.average(jkgamma[la])
  var = (Nmeas - 1.) * np.sum((jkgamma[la] - ga_ave)**2) / float(Nmeas)
  print "%.6g %.6g %.4g %.4g" % (x, ga_ave, de_x, np.sqrt(var))

runtime += time.time()
print "# Runtime: %.2g seconds" % runtime

# Print list of files
print "# Files analyzed:"
for filename in files:
  print "#", filename
# ------------------------------------------------------------------
