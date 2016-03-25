#!/usr/bin/env python
import os
import sys
import numpy as np
# ------------------------------------------------------------------
# Estimate mass and errors from cosh effective mass

# Parse argument file to analyze
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
filename = str(sys.argv[1])
if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)

for line in open(filename):
  # Looking for "# meff = dat +/- err"
  if line.startswith('# meff = '):
    temp = line.split()
    in_dat = float(temp[3])
    in_err = float(temp[5])
    out_dat = np.arccosh(in_dat)
    out_err = (abs(out_dat - np.arccosh(in_dat - in_err))
            +  abs(out_dat - np.arccosh(in_dat + in_err))) / 2
    print "M = %.4g +/- %.2g" % (out_dat, out_err)
    sys.exit(0)   # Done
# ------------------------------------------------------------------
