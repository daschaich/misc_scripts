#!/usr/bin/python
import glob
import math
import os
import sys
# ------------------------------------------------------------------
# This script extracts S4b order parameters from Sout files
# in the current directory

# First make sure target doesn't already exist
if os.path.isfile('S4b.dat'):
  print "ERROR: S4b.dat already exists"
  sys.exit(1)

# Construct and sort list of output files
cfgs = []
temp = filetag + '.*'
for filename in glob.glob(temp):
  cfgs.append(int((filename.split('.'))[-1]))  # Number after last .
cfgs.sort()

# Check to make sure the arguments are appropriate
if len(cfgs) == 0:
  print "ERROR: no files named", temp
  sys.exit(1)

# Extract lattice volume from first output file
firstFile = filetag + '.' + str(cfgs[0])
for line in open(firstFile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  elif line.startswith('iseed '):
    break   # Done scanning through file
volume = float(L**3 * Nt)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Construct and print euclidean means of plaquette and link differences
outfilename = 'S4b.dat'
outfile = open(outfilename, 'w')
for i in cfgs:
  pd = 0.0
  ld = 0.0
  print >> outfile, i,
  toOpen = filetag + '.' + str(i)
  for line in open(toOpen):
    if line.startswith('StaggPlaq t '):
      pd += (float((line.split())[2]) - float((line.split())[3]))**2
    elif line.startswith('StaggPlaq x '):
      pd += (float((line.split())[2]) - float((line.split())[3]))**2
    elif line.startswith('StaggPlaq y '):
      pd += (float((line.split())[2]) - float((line.split())[3]))**2
    elif line.startswith('StaggPlaq z '):
      pd += (float((line.split())[2]) - float((line.split())[3]))**2
    elif line.startswith('pbpt: '):
      ld += (float((line.split())[3]) - float((line.split())[4]))**2
    elif line.startswith('pbpx: '):
      ld += (float((line.split())[3]) - float((line.split())[4]))**2
    elif line.startswith('pbpy: '):
      ld += (float((line.split())[3]) - float((line.split())[4]))**2
    elif line.startswith('pbpz: '):
      ld += (float((line.split())[3]) - float((line.split())[4]))**2

  print >> outfile, math.sqrt(pd), math.sqrt(ld) / volume
outfile.close()
# ------------------------------------------------------------------
