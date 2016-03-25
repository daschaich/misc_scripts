#!/usr/bin/python
import glob
import os
import sys
# ------------------------------------------------------------------
# Label eigenvalue output files by the number of eigenvalues calculated
# Target format: eig_Nvecs_start_vol_beta_ratio_mass.#
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# First deal with files that don't yet have Nvecs in their name
# Remove any symlinks that I set up in the past to run these files
# through my file-parsing scripts
for start in ['low', 'high', 'mix']:
  allFiles = 'eig_' + start + '_*'
  for filename in glob.glob(allFiles):
    if os.path.islink(filename):
      print "Removing link", filename
      os.remove(filename)
      continue    # Done with this file, which no longer exists

    # Extract number of eigenvalues
    Nvecs = -1
    for line in open(filename):
      if line.startswith('Nvecs ') \
      or line.startswith('Number_of_eigenvals '):
        Nvecs = (line.split())[1]
        break   # Done scanning through file

    if Nvecs == -1:
      print "ERROR: Nvecs not defined in", filename
#      sys.exit(1)
      continue

    # Split filename, reorganize and move
    temp = filename.split('_')
    # start is temp[1]
    vol = temp[2]
    beta = temp[3]
    ratio = temp[4]
    massN = temp[5]
    newname = 'eig_' + Nvecs + '_' + start + '_' + vol + '_' + beta \
              + '_' + ratio + '_' + massN

    if os.path.isfile(newname):
      print "ERROR:", newname, "already exists"
#      sys.exit(1)
      continue
    print filename, "-->", newname
    os.rename(filename, newname)
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Now fix eig_Nvecs_vol_beta_ratio_mass_start format
  allFiles = 'eig_*' + start + '.[0-9]*'
  for filename in glob.glob(allFiles):
    if os.path.islink(filename):
      print "Removing link", filename
      os.remove(filename)
      continue    # Done with this file, which no longer exists

    # Split filename, reorganize and move
    temp = filename.split('_')
    Nvecs = temp[1]
    vol = temp[2]
    beta = temp[3]
    ratio = temp[4]
    mass = temp[5]
    startN = temp[6]
    N = startN.split('.')[1]

    # Make sure we don't mess up eig_*block* files
    # Those should all have at least one more label,
    # pushing the start tag to later in the filename
    if start == startN.split('.')[0]:
      newname = 'eig_' + Nvecs + '_' + start + '_' + vol + '_' + beta \
                + '_' + ratio + '_' + mass + '.' + N

      if os.path.isfile(newname):
        print "ERROR:", newname, "already exists"
#        sys.exit(1)
        continue
      print filename, "-->", newname
      os.rename(filename, newname)
  # ----------------------------------------------------------------



  # ----------------------------------------------------------------
  # Since we're going through this trouble, let's deal with the
  # MCRG-blocked output files as well
  allFiles = 'mcrg_' + start + '_*'
  thisdir = os.getcwd()
  if 'Run_a-0.4' in thisdir:
    ratio = '-0.4'
  elif 'Run_un' in thisdir or 'Run_plaq' in thisdir:
    ratio = '0.0'
  else:
    ratio = '-0.25'
  for filename in glob.glob(allFiles):
    # Split filename, reorganize and move
    temp = filename.split('_')
    # start is temp[1]
    vol = temp[2]
    beta = temp[3]
    massN = temp[4]
    if massN != ratio:
      newname = 'mcrg_' + start + '_' + vol + '_' + beta \
                + '_' + ratio + '_' + massN

      if os.path.isfile(newname):
        print "ERROR:", newname, "already exists"
#        sys.exit(1)
        continue
      print filename, "-->", newname
      os.rename(filename, newname)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Finally, let's check that we really are done
for start in ['low', 'high', 'mix']:
  allFiles = 'eig_' + start + '_*'
  print "ls", allFiles
  for filename in glob.glob(allFiles):
    print filename
  allFiles = 'eig_*' + start + '.[0-9]*'
  print "ls", allFiles
  for filename in glob.glob(allFiles):
    if start == (((filename.split('_'))[6]).split('.'))[0]:
      print filename
  allFiles = 'mcrg_' + start + '_* | grep -v ' + ratio
  print "ls", allFiles
  for filename in glob.glob(allFiles):
    if ratio != (filename.split('_'))[4]:
      print filename
# ------------------------------------------------------------------
