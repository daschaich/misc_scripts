#!/usr/bin/python
import os
import sys
import glob
import math
# ------------------------------------------------------------------
# This script extracts the real part of the Wilson-flowed Polyakov loop
# after flow times c=sqrt(8t)/L=0.2, 0.3, 0.4 and 0.5
# Needs to be called from the directory with the output files

# First make sure target doesn't already exist
if os.path.isfile('Wpoly.dat'):
  print "ERROR: Wpoly.dat already exists"
  sys.exit(1)

# Construct and sort list of output files
cfgs = []
temp = 'Wflow.*'
for filename in glob.glob(temp):
  cfgs.append(int((filename.split('.'))[-1]))  # Number after last .
cfgs.sort()

# Check to make sure the arguments are appropriate
if len(cfgs) == 0:
  print "ERROR: no files named", temp
  sys.exit(1)

# Extract lattice volume from first output file
firstFile = 'Wflow.' + str(cfgs[0])
for line in open(firstFile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  elif line.startswith('iseed '):
    break   # Done scanning through file
if L > Nt:
  L = Nt    # Take minimum
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now just grab and print for c=sqrt(8t)/L=0.2, 0.3, 0.4 and 0.5
# Encountered some round-off issues for Nt=24...
outfilename = 'Wpoly.dat'
outfile = open(outfilename, 'w')
for i in cfgs:
  print >> outfile, i,
  toOpen = 'Wflow.' + str(i)
  check = -1
  for line in open(toOpen):
    if line.startswith('POLYA ORIG '):
      t = float((line.split())[2])
      c = math.sqrt(8 * t) / L
      if c > 0.19 and c < 0.21:
        print >> outfile, float((line.split())[3]),
      elif c > 0.29 and c < 0.31:
        print >> outfile, float((line.split())[3]),
      elif c > 0.39 and c < 0.41:
        print >> outfile, float((line.split())[3]),
      elif c > 0.49 and c < 0.51:
        print >> outfile, float((line.split())[3])
    elif line.startswith('RUNNING COMPLETED'):
      check = 1
  if check == -1:
    print toOpen[0], "did not complete"
    sys.exit(1)
outfile.close()
# ------------------------------------------------------------------
