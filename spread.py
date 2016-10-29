#!/usr/bin/python
import sys
# ------------------------------------------------------------------
# Extract the maximum and minimum from a list of numbers
# Print their average and spread around average

# Parse arguments: at least two numbers
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<numbers>"
  sys.exit(1)
dat = []
for i in range(1, len(sys.argv)):
  dat.append(float(sys.argv[i]))

big = max(dat)
small = min(dat)
ave = 0.5 * (big + small)
spread = 0.5 * (big - small)
print "%.4g %.4g" % (ave, spread)
# ------------------------------------------------------------------
