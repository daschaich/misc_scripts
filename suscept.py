#!/usr/bin/python
import glob
import math
import os
import sys
# ------------------------------------------------------------------
# Jackknife analysis of chiral susceptibilities
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Standard error from list (since beowulf doesn't have numerical python)
def std(x):
  N = len(x)
  mean = sum(x) / float(N)
  err = 0
  for i in x:
    err += (i - mean)**2
  err = math.sqrt(err / float(N * (N - 1)))
  return err
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Parse arguments: first is ensemble tag, second is thermalization cut,
# third is autocorrelation time
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<tag> <cut> <bin>" 
  sys.exit(1)
tag = str(sys.argv[1])
cut = int(sys.argv[2])
bin = int(sys.argv[3])
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First make sure we're calling this from the right place
if not os.path.isdir('Out'):
  print "ERROR: Out/ does not exist"
  sys.exit(1)

# Construct and sort list of output files
files = 'Out/susc_' + tag + '.'
allFiles = files + '*'
cfgs = []
for filename in glob.glob(allFiles):
  cfg = int((filename.split('.'))[-1])  # Number after last .
  if cfg >= cut:
    cfgs.append(cfg)
cfgs.sort()

# Extract lattice volume from first output file
firstFile = files + str(cfgs[0])
for line in open(firstFile):
  if line.startswith('nx '):
    L = int((line.split())[1])
  elif line.startswith('nt '):
    Nt = int((line.split())[1])
  elif line.startswith('iseed '):
    break   # Done scanning through file
vol = L**3 * Nt
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Go through output files
# Record two susceptibility terms, and chiral condensate
# Note we discard any partial bin at the end
iter = 0
first_dat = 0     # Accumulated within each bin
second_dat = 0
pbp_dat = 0
first_all = []    # Record data from all bins
second_all = []   # For later construction of jackknife samples
pbp_all = []
for i in cfgs:
  toOpen = files + str(i)
  for line in open(toOpen):
    if line.startswith('susc: '):
      temp = line.split()   # Convert line into list
      first_dat += float(temp[1])
      second_dat += float(temp[2])
    elif line.startswith('pbp: '):
      temp = line.split()
      pbp_dat += float(temp[3]) + float(temp[4])

  # Accumulate results every bin
  iter += 1
  if iter == bin:
    first_all.append(first_dat / iter)
    second_all.append(second_dat / iter)
    pbp_all.append(pbp_dat / iter)
    first_dat = 0
    second_dat = 0
    pbp_dat = 0
    iter = 0
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Now we can print the mean and construct the jackknife error
# Cf. www.physics.utah.edu/~detar/phycs6730/handouts/jackknife/jackknife/
print "# Thermalization cut at configuration", cut, "of", cfgs[-1]
print "#", len(cfgs), "files -->", len(first_all), "measurements"
N = len(first_all)
first_tot = sum(first_all)  # Let's only compute these once
second_tot = sum(second_all)
pbp_tot = sum(pbp_all)
first = first_tot / float(N)
second = second_tot / float(N)
temp = pbp_tot / float(N)
pbp = temp**2 / float(vol)   # pbp contribution to chi
full = first + second - pbp
print "# Individual terms:",
print "%.2f(%.2f)," % (first, std(first_all)),
print "%.2f(%.2f)," % (second, std(second_all)),
print "%.2f(%.2f)" % (pbp, 2 * pbp_tot * std(pbp_all) / float(N * vol))
print "# Overall result:",
print "%.2f + %.2f - %.2f =" % (first, second, pbp),

err = 0
for i in range(0, N):
  first = (first_tot - first_all[i]) / float(N - 1)
  second = (second_tot - second_all[i]) / float(N - 1)
  temp = (pbp_tot - pbp_all[i]) / float(N - 1)
  pbp = temp**2 / float(vol)
  sample = first + second - pbp
  err += (sample - full)**2
print "%.2f(%.2f)" % (full, math.sqrt(float(N - 1) * err / float(N)))
# ------------------------------------------------------------------
