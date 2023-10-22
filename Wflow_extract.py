#!/usr/bin/python
import os
import sys
import glob
import numpy as np
# ------------------------------------------------------------------
# This script extracts the Wilson flow running coupling,
# mimicking parse_BSM.pl without relying on out.* files
# No arguments; for now fix c=0.2, 0.25, 0.3 and 0.35
# Also extract topological charge for c=0.2, 0.3, 0.4 and 0.5

# See if this is necessary
if os.path.isfile('data/Wflow.csv'):
  print "Just use data/Wflow.csv"
  sys.exit(0)
if os.path.isfile('data/topo.csv'):
  print "Just use data/topo.dat"
  sys.exit(0)

# Construct and sort list of output files
cfgs = []
for filename in glob.glob('Out/Wflow.[0-9]*'):
  cfg = int((filename.split('.'))[1])             # Number after '.'
  if cfg not in cfgs:
    cfgs.append(cfg)
cfgs.sort()

# Check to make sure the arguments are appropriate
if len(cfgs) == 0:
  print "ERROR: no files named Out/Wflow.[0-9]*"
  sys.exit(1)

# Extract lattice volume from path
# For now assume L and Nt are both two-digit numbers
path = os.getcwd()
path = path.replace('mnt', '')    # Accommodate Barkla filesystem
temp = path.split('nt')
L = int(temp[0][-2:])    # Last two digits before 'nt'
Nt = int(temp[1][:2])    # First two digits after 'nt'
vol = L**3 * Nt

# Proportionality factor, ignoring finite-volume correction delta_c
if 'SU4' in path:
  gprop = 128. * np.pi**2 / (3. * 15.);
else:
  gprop = 128. * np.pi**2 / (3. * 8.);
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now just grab and print, overwriting existing file
flowfilename = 'data/Wflow.dat'
flowfile = open(flowfilename, 'w')
topofilename = 'data/topo.dat'
topofile = open(topofilename, 'w')
gSq = np.zeros(4, dtype = np.float)
topo = np.zeros(4, dtype = np.float)
for i in cfgs:
  # Reset
  topprop = 1
  gSq[:] = 0;     topo[:] = 0;      c = -1.;      cOld = 1.
  toOpen = 'Out/Wflow.' + str(i)
  for line in open(toOpen):
    # This is annoying: measurements at Fermilab have the wrong normalization
    # Measurements at Argonne or Livermore have the right normalization
    if 'fnal.gov' in line:
      topprop = float(vol) / (4. * np.pi**2)**2
      continue
    elif not line.startswith('WFLOW '):
      continue

    # Format: WFLOW  t  plaq  E  t^2*E  t*d(t^2*E)/dt  12t^2*(3-plaq)  topo
    temp = line.split()
    c = np.sqrt(8. * float(temp[1])) / float(L)
    if cOld < 0.2 and c >= 0.2:
      gSq[0] = gprop * float(temp[4])
      topo[0] = topprop * float(temp[7])
    elif cOld < 0.25 and c >= 0.25:
      gSq[1] = gprop * float(temp[4])
    elif cOld < 0.3 and c >= 0.3:
      gSq[2] = gprop * float(temp[4])
      topo[1] = topprop * float(temp[7])
    elif cOld < 0.35 and c >= 0.35:
      gSq[3] = gprop * float(temp[4])
    elif cOld < 0.4 and c >= 0.4:
      topo[2] = topprop * float(temp[7])
    elif cOld < 0.5 and c >= 0.5:
      topo[3] = topprop * float(temp[7])
    cOld = c

  print >> flowfile, i,
  print >> topofile, i,
  for dat in gSq:
    print >> flowfile, "%.6g" % dat,
  for dat in topo:
    print >> topofile, "%.6g" % dat,
  print >> flowfile, "\n",
  print >> topofile, "\n",

flowfile.close()
topofile.close()
# ------------------------------------------------------------------
