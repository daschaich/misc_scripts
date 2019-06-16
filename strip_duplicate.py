#!/usr/bin/python
import os
import sys
import glob
from subprocess import call
# ------------------------------------------------------------------
# Strip duplicate output from files where HMC accidentally run twice
# Should have no effect on files that are already okay

# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Cycle over measurement files
for filename in glob.glob('Out/out.*'):
  outfile = open('TEMP', 'w')
  check = -1
  for line in open(filename):
    if line.startswith('=== MPI application finished '):
      check = 1

    # Only start copying file after first 'application finished' tag
    if check == 1
      print >> outfile, line.rstrip()

  if check == -1:
    print filename, "did not complete"
    continue
  outfile.close()
  os.rename('TEMP', filename)
# ------------------------------------------------------------------

