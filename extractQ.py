#!/usr/bin/python
import glob
import os
import sys
# ------------------------------------------------------------------
# This script extracts the topological charge from Wilson flow output,
# flowing out to c=sqrt(8t)/L=0.5 --> t=(L*0.5)^2/8
# and optionally correcting the normalization

# Parse arguments: first is file name, up to # at the end
# (require this to be called from the output directory)
# Second is a flag that indicates we should fix the normalization
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<tag [out.Wflow]>",
  print "[fix normalization (optional)]"
  sys.exit(1)
filetag = str(sys.argv[1])

norm = 1
if len(sys.argv) > 2:
  norm = 1 / (16.0 * 3.14159 * 3.14159 * 3.14159 * 3.14159)   # 1 / (4pi^2)^2
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure target doesn't already exist
if os.path.isfile('Q.dat'):
  print "ERROR: Q.dat already exists"
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
vol = L**3 * Nt
if len(sys.argv) > 2:
  norm *= vol   # volume / (4pi^2)^2

# Set t at which c=sqrt(8t)/L=0.5
# Set the maximum number of RG blockings
if L <= Nt:
  tmax = L * 0.5 * L * 0.5 / 8.
else:
  tmax = Nt * 0.5 * Nt * 0.5 / 8.
target = 'WFLOW ' + str(tmax)
if '.0' in target:
  target = target[:-2]   # Strip '.0' from end of target
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now just grab and print properly normalized topological charge
outfilename = 'Q.dat'
outfile = open(outfilename, 'w')
for i in cfgs:
  toOpen = filetag + '.' + str(i)
  for line in open(toOpen):
    if line.startswith(target):
      toprint = float((line.split())[7]) * norm
      print >> outfile, "%d %.4g" % (i, toprint)
outfile.close()
# ------------------------------------------------------------------
