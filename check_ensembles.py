#!/usr/bin/python
import math
import os
import sys
# ------------------------------------------------------------------
# Print out summary of all ensembles in current directory

# First make sure we're calling this from the right place
if not os.path.isfile('ensembles'):
  print "ERROR: ensembles list does not exist"
  sys.exit(1)
elif not os.path.isdir('data'):
  print "ERROR: data/ does not exist"
  sys.exit(1)

tags = []
for ensemble in open('ensembles'):
  tags.append(ensemble.rstrip())

# Print header
print "|---------------------------------------------------------------------|"
print "| Ensemble (Total: %3d)       | errors    | Wflow     | MDTU  | therm | Notes" % len(tags)
print "|---------------------------------------------------------------------|"

# Now cycle through ensembles, determining errors, Wflow, MCRG and MDTU
errors = 0
Wflow = 0
therm = 0
note = 0
for tag in tags:
  fulltag = tag + '.'   # Avoid false positives
  errorfile = 'ERRORS.' + tag
  if not os.path.isfile(errorfile):   # Check presence of error file
    print "ERROR: missing file ERRORS.%s" % tag
    continue
#    sys.exit(1)
  for line in open(errorfile):
    if fulltag in line:
      errors += 1

  missing = 'MISSING.' + tag
  if not os.path.isfile(missing):
    print "ERROR: missing file MISSING.%s" % tag
    continue
#    sys.exit(1)
  # Now check for missing Wilson flow measurements
  for line in open(missing):
    if 'Wflow' in line:
      Wflow += 1

  # Grab MDTU from the appropriate file in data/
  TUfile = 'data/TU.' + tag + '.csv'
  temp = ((os.popen("tail -n 1 %s" % TUfile).read()).split(','))[1]
  if 'MDTU' in temp:
    MDTU = 0
  else:
    MDTU = math.floor(float(temp))

  # Check to see if we have set a thermalization cut for this ensemble
  if os.path.isfile('therm'):
    for line in open('therm'):
      if line.startswith('python') or line.startswith('  python'):
        temp = line.split()
        if tag in temp[2] and temp[2] in tag and int(temp[4]) < 500:
          therm = int(temp[3])
          if line.rstrip().endswith('Placeholder'):
            note = 1
            toprint = "Placeholder thermalization cut"

  # Check to see if we have any notes for this ensemble
  # This determines where to end the line
  # Note specific name of file and tag:notes format of its lines
  if os.path.isfile('NOTES'):
    for line in open('NOTES'):
      if line.startswith(tag):
        note = 1
        toprint = ((line.split(':'))[1]).rstrip()
        break

  # Print ensemble name, errors, Wflow, MCRG and MDTU
  # ljust determined by longest tag, in 8f 16nt32
  print "| %s" % tag.ljust(27),
  if errors > 0:
    print "| %3d TODO " % errors,
  else:
    print "|     DONE ",
  if Wflow > 0:
    print "| %3d TODO " % Wflow,
  else:
    print "|     DONE ",
  print "| %5d" % MDTU,
  if therm == 0:
    print "| unset",
  else:
    print "| %5d" % therm,
  if note == 0:
    print "|"
  else:
    print "| %s" % toprint

  # Reset for next ensemble
  errors = 0
  Wflow = 0
  therm = 0
  note = 0
print "|---------------------------------------------------------------------|"
errorfile = 'ERRORS.CFGs'
cfg_errors = 0
if not os.path.isfile(errorfile):   # Check presence of error file
  print "ERROR: missing file ERRORS.CFGs"
for line in open(errorfile):
  cfg_errors += 1
print "| Configuration time stamp mismatches: %5d                          |" % cfg_errors
print "|---------------------------------------------------------------------|"
# ------------------------------------------------------------------
