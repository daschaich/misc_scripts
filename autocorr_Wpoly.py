#!/usr/bin/python
import os
import sys
import glob
import numpy as np
import acor         # Uses "the Kubo formula" to compute autocorrelation time
# ------------------------------------------------------------------
# Compute and print autocorrelation (and resulting mean, err)
# of (c=0.5 Wilson-flowed) Polyakov loop magnitudes
# Parse dygraph data file with given thermalization cut

# Parse argument: Thermalization cut
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<cut>"
  sys.exit(1)
cut = int(sys.argv[1])

# First make sure we're calling this from the right place
if not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

# For now consider only modulus of Wilson-flowed Polyakov loop with c=0.5
# (Format: MDTU,c=0.2,0.3,0.4,0.5)
dat = []
for line in open('data/Wpoly_mod.csv'):
  if line.startswith('M'):
    continue
  temp = line.split(',')
  MDTU = float(temp[0])
  if MDTU <= cut:
    continue
  dat.append(float(temp[4]))

# !!! Assume each measurement separated by 10 MDTU
# Print out autocorrelation time in MDTU,
# along with effective number of independent measurements
tau, mean, sigma = acor.acor(np.array(dat))
eff_stat = np.floor(len(dat) / tau)

print "Analyzing %d thermalized measurements" % len(dat)
print "tau = %d MDTU" % (tau * 10.0)
print "%.8g %.4g # %d" % (mean, sigma, eff_stat)
# ------------------------------------------------------------------
