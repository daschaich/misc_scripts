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
# Annoyingly, need to use absolute Barkla path
path = "/users/schaich/LSD/SU4/"
if not os.path.isdir(path):
  print("ERROR: " + path + " not found")
  sys.exit(1)
os.chdir(path)
f = h5py.File("SU4_data.h5", 'w')

# Top-level groups for each Nf
for Nf in glob.glob('*f'):
  f.create_group(Nf)

  # Second- and third-level groups for each (Nt, L) volume
  os.chdir(path + Nf)
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
    os.chdir(path + Nf + '/' + vol)
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
      # First pass through plaq.dat to set Ntraj
      os.chdir(path + Nf + '/' + vol + '/' + stream)
      traj_arr = []
      plaq_arr = [[], []]
      for line in open('data/plaq.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        traj_arr.append(int(temp[0]))
        plaq_arr[0].append(float(temp[1]))    # ss
        plaq_arr[1].append(float(temp[2]))    # st
      Ntraj = len(traj_arr)
      traj = np.array(traj_arr)
      plaq = np.array(plaq_arr)
      f.create_dataset(this_ens + '/traj', data=traj)
      f.create_dataset(this_ens + '/plaq', data=plaq)

      # TODO: Attributes and results... set Nblocks

      # Now all other observables measured every trajectory=MDTU
      # Dynamically figure out how many data on each line
      # TODO: Add others... pbp header problematic...
      for obs in ['exp_dS']:
        name = this_ens + '/' + obs
        obsfile = 'data/' + obs + '.csv'
        traj = 0
        for line in open(obsfile):
          temp = line.split(',')
          if line.startswith('M') or line.startswith('t'):
            Ndat = len(temp)
            dset = f.create_dataset(name, (Ntraj,Ndat,), dtype='f')
            continue
          for j in range(Ndat - 1):       # Skipping MDTU label
            dset[traj][j] = float(temp[j + 1])
          traj += 1

        # TODO: Attributes and results... check Nblocks

      # First pass through topo.dat to set NWflow
      # Only recording c=0.5 topological charge for clean time-series plot
      Wmeas_arr = []
      topo_arr = []
      for line in open('data/topo.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        Wmeas_arr.append(int(temp[0]))
        topo_arr.append(float(temp[4]))
      NWflow = len(Wmeas_arr)
      Wmeas = np.array(Wmeas_arr)
      topo = np.array(topo_arr)
      f.create_dataset(this_ens + '/Wmeas', data=Wmeas)
      f.create_dataset(this_ens + '/topo', data=topo)

      # Now all other observables measured every trajectory=MDTU
      # Try to figure out how many data to include from each line
      # TODO: Add others... (different c for Wflow vs. Wpoly...)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Zero-temperature runs are a bit simpler
os.chdir(path + "4f/24nt48")
for ens in glob.glob('b[1-9]*'):
  temp = stream.split('_')
  beta = temp[0]
  mass = temp[1]
  f.create_group('4f/24nt48/' + mass + '/' + beta)
# ------------------------------------------------------------------
