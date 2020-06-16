#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import h5py
# ------------------------------------------------------------------
# Package SU(4) data and results/attributes into publishable HDF5 file

# Cycle over all streams and write to ~/LSD/SU4/SU4_data.hdf5
# Attributes for each stream:
#   Nf, L, Nt, beta_F, (valence) mass,
#   thermalization cut, block size
# Data for each stream:
#   plaquette, chiral condensate, ...
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Start working in desired place, if it exists
if not os.path.isdir("~/LSD/SU4"):
  print("ERROR: ~/LSD/SU4 not found")
  sys.exit(1)
os.chdir("~/LSD/SU4")
f = h5py.File("SU4_data.h5", 'w')

# Top-level groups for each Nf
for Nf in glob.glob('*f'):
  f.create_group(Nf)

  # Second- and third-level groups for each (Nt, L) volume
  os.chdir("~/LSD/SU4/" + Nf)
  for vol in glob.glob('*nt*'):
    if vol == '24nt48':       # Do zero-temperature runs separately
      continue

    temp = vol.split('nt')
    L = 'L' + temp[0]
    Nt = 'Nt' + temp[1]
    this_vol = Nf + '/' + Nt + '/' + L
    f.create_group(this_vol)

    # Fourth- and fifth-level groups for each (beta_F, mass) ensemble
    # For 0f ensembles this is the valence mass for pbp
    # Then sixth-level groups for each stream ('high', 'low', 'combo')
    os.chdir("~/LSD/SU4/" + Nf + '/' + vol)
    for stream in glob.glob('b[1-9]*'):   # Avoid 'backup' dirs
      temp = stream.split('_')
      beta = temp[0]
      start = temp[1]
      mass = temp[2]
      this_ens = this_vol + '/' + mass + '/' + beta
      f.create_group(this_ens + '/' + start)
      if start == 'low':
        f.create_group(this_ens + '/combo')

      # Now set up data and attributes/results for this stream
      # First get numbers of trajectories and Wilson flow measurements
      os.chdir("~/LSD/SU4/" + Nf + '/' + vol + '/' + stream)
      for line in open('data/plaq.csv'):
        Ntraj = int(line.split(',')[0])
      for line in open('data/Wpoly.csv'):
        NWflow = int(line.split(',')[0])

      # First observables measured after every trajectory=MDTU
      # Try to figure out how many data to include from each line
      for obs in ['plaq']:
        name = this_ens + '/' + obs
        obsfile = 'data/' + obs + '.csv'
        i = 0
        for line in open(obsfile):
          temp = line.split(',')
          if line.startswith('M'):
            Ndat = len(temp)
            dset = f.create_dataset(name, (Ntraj,Ndat,), dtype='f')
            continue
          for j in range(Ndat):
            dset[i][j] = float(temp[j + 1])
          i += 1
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Zero-temperature runs are a bit simpler
os.chdir("~/LSD/SU4/4f/24nt48")
for ens in glob.glob('b[1-9]*'):
  temp = stream.split('_')
  beta = temp[0]
  mass = temp[1]
  f.create_group('4f/24nt48/' + ens)
# ------------------------------------------------------------------
