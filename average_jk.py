#!/usr/bin/python
import os
import sys
import numpy as np
# ------------------------------------------------------------------
# Extract final central value and uncertainty
# from recorded jackknife sample results

# Parse argument: the file to analyze (FORMAT: block# dat)
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<file>"
  sys.exit(1)
filename = str(sys.argv[1])

if not os.path.isfile(filename):
  print "ERROR:", filename, "does not exist"
  sys.exit(1)

# Read, parse and process data
# Assumed format: block# dat
datList = []
for line in open(filename):
  if len(line) == 1 or line.startswith('#') or line.startswith('!'):
    continue
  temp = line.split()
  datList.append(float(temp[1]))
dat = np.array(datList)
nblocks = len(datList)

# Average over jackknife samples and print out results
ave = np.average(dat)
var = (nblocks - 1.) * np.sum((dat - ave)**2) / float(nblocks)
print "%.6g %.4g" % (ave, np.sqrt(var))
# ------------------------------------------------------------------
