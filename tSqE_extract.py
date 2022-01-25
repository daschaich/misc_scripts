#!/usr/bin/python
import glob
import math
import os
import sys
# ------------------------------------------------------------------
# This script extracts t^2E from Wilson flow output,
# flowing out to c=sqrt(8t)/L=0.5 --> t=(L*0.5)^2/8

# Some differences compared to Wflow_extract.py
# * Allow different Wflow file names
# * Assume correct normalization (for meas at Argonne, Livermore, Liverpool)
# * Fix c=0.2, 0.3, 0.4 and 0.5 (no 0.25 or 0.35)

# Parse arguments: file name, up to # at the end
# (require this to be called from the output directory)
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<tag [out.Wflow]>",
  sys.exit(1)
filetag = str(sys.argv[1])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure target doesn't already exist
if os.path.isfile('tSqE.dat'):
  print "ERROR: tSqE.dat already exists"
  sys.exit(1)
if os.path.isfile('topo.dat'):
  print "ERROR: topo.dat already exists"
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
if L > Nt:
  L = Nt    # Take minimum
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now just grab and print for c=sqrt(8t)/L=0.2, 0.3, 0.4 and 0.5
outfilename = 'tSqE.dat'
outfile = open(outfilename, 'w')
topofilename = 'topo.dat'
topofile = open(topofilename, 'w')
for i in cfgs:
  cOld = -1.0
  print >> outfile, i,
  toOpen = filetag + '.' + str(i)
  for line in open(toOpen):
    if line.startswith('WFLOW '):
      t = float((line.split())[1])
      c = math.sqrt(8.0 * t) / L
      if cOld < 0.2 and c >= 0.2:
        print >> outfile, float((line.split())[4]),
        print >> topofile, float((line.split())[7]),
      elif cOld < 0.3 and c >= 0.3:
        print >> outfile, float((line.split())[4]),
        print >> topofile, float((line.split())[7]),
      elif cOld < 0.4 and c >= 0.4:
        print >> outfile, float((line.split())[4]),
        print >> topofile, float((line.split())[7]),
      elif cOld < 0.5 and c >= 0.5:
        print >> outfile, float((line.split())[4])
        print >> topofile, float((line.split())[7])
      cOld = c
outfile.close()
topofile.close()
# ------------------------------------------------------------------
