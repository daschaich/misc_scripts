#!/usr/bin/python
import glob
import math
import os
import sys
# ------------------------------------------------------------------
# This script extracts pbp from QHMC output files,
# assuming some naming conventions that may be Argonne-specific

# Parse argument: destination directory
# third is autocorrelation time
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<destination dir>"
  sys.exit(1)
dest = str(sys.argv[1])
# ------------------------------------------------------------------




# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('done'):
  print "ERROR: done/ does not exist"
  sys.exit(1)
if not 'runs' in os.getcwd():
  print "ERROR: need to call from a 'runs' directory"
  sys.exit(1)

# Check target
if not os.path.isdir(dest):
  print "ERROR:", dest, "does not exist"
  sys.exit(1)
outfilename = dest + 'pbp.dat'
if os.path.isfile(outfilename):
  print "ERROR:", outfilename, "already exists"
  sys.exit(1)

# Construct and sort list of output files
# Avoid duplicate output files, one in the done/ directory
cfgs = []
for filename in glob.glob('done/job*.output'):
  temp = (filename.split('.'))[0]             # Everything before .
  cfg = int((temp.split('job'))[1])           # Number after 'job'
  if cfg not in cfgs:
    cfgs.append(cfg)
for filename in glob.glob('job*.output'):
  temp = (filename.split('.'))[0]             # Everything before .
  cfg = int((temp.split('job'))[1])           # Number after 'job'
  if cfg not in cfgs:
    cfgs.append(cfg)
cfgs.sort()

# Check to make sure the arguments are appropriate
if len(cfgs) == 0:
  print "ERROR: no files named", temp
  sys.exit(1)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now just grab and print
print "Warning: Assuming unit-length trajectories",
print "with pbp measured after each trajectory"
toprint = []
for i in cfgs:
  # Starting lattice may not be labelled 0 (it is often 6)
  MDTU = cfgs[0] + len(toprint)

  # This is annoying:
  # The next job may load the wrong lattice
  # There may be duplicate output files, one in the done/ directory
  # Preferentially load the output file in the current directory
  toOpen = 'job' + str(i) + '.output'
  if not os.path.isfile(toOpen):
    toOpen = 'done/job' + str(i) + '.output'
    if not os.path.isfile(toOpen):
      print "ERROR: can't find", toOpen
      sys.exit(1)

  # If job loaded the wrong lattice (MDTU > i)
  # pop off the last (MDTU - i) points
  if MDTU > i:
    print "Warning: Opening %s after counting %d MDTU" % (toOpen, MDTU)
    count = 0
    while MDTU > i:
      toprint.pop()
      MDTU = cfgs[0] + len(toprint)
      count += 1
    print "         Discarding previous %d measurements" % count
  # If we are missing an output file, just skip to next
  elif MDTU < i:
    print "Warning: Opening %s after counting %d MDTU" % (toOpen, MDTU)
    while MDTU < i:
      MDTU += 1
      toprint.append("null")
    print "         Skipping to %d MDTU" % MDTU

  temp = []
  for line in open(toOpen):
    if line.startswith('pbp mass '):
      pbp = float((line.split())[4])
      temp.append(pbp)

    # This is annoying:
    # There may be extra pbp measurements after the last lattice is saved
    # which should be ignored
    if line.startswith('saving lattice file '):
      for X in temp:
        MDTU += 1
        toprint.append(X)
      temp = []  # Reset

# Print to file, overwriting current version
outfile = open(outfilename, 'w')
MDTU = cfgs[0]
for pbp in toprint:
  MDTU += 1   # pbp printed after accept/reject step
  print >> outfile, MDTU, pbp
outfile.close()
# ------------------------------------------------------------------
