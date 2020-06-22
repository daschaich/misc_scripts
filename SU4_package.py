#!/usr/bin/python3
import os
import sys
import glob
import numpy as np
import h5py
# ------------------------------------------------------------------
# Package SU(4) data and results/attributes into publishable HDF5 file

# File will hopefully be smaller than ~5.5GB total disk usage
# $ du -h --total *f/*nt*/b[1-9]*/data *f/*nt*/b[1-9]*/results | tail -n 1
#   5.6G total

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
#  for vol in glob.glob('*nt*'): # !!!TODO: Accelerate for testing
  for vol in glob.glob('*nt12'):
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
      toCheck = stream + '/results/Wpoly_mod.autocorr'
      if not os.path.isfile(toCheck):     # Skip unfinished streams
        continue

      temp = stream.split('_')
      beta = temp[0]
      start = temp[1]
      mass = temp[2]
      this_ens = this_vol + '/' + mass + '/' + beta
      this_str = this_ens + '/' + start
      this_grp = f.create_group(this_str)
      if start == 'low':
        combine = True
        f.create_group(this_ens + '/combo')

      # ------------------------------------------------------------
      # Record acceptance and Nblocks for stream
      os.chdir(path + Nf + '/' + vol + '/' + stream)
      for line in open('results/accP.dat'):
      # Now set up data and attributes/results for this stream
      # First do plaq, setting Ntraj & Nblocks
      plaq_arr = [[], []]
      for line in open('data/plaq.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        plaq_arr[0].append(float(temp[1]))    # ss
        plaq_arr[1].append(float(temp[2]))    # st
      Ntraj = len(plaq_arr[0])
      this_grp.attrs['Ntraj'] = Ntraj

      plaq = np.array(plaq_arr)
      this_grp.create_dataset('traj', data=traj)
      dset = this_grp.create_dataset('plaq', data=plaq)
      dset.attrs['columns'] = ['ss', 'st']

      # Set Nblocks and record results as attributes
      for line in open('results/plaq.dat'):
        temp = line.split()
        Nblocks = int(temp[-1])
        this_grp.attrs['Nblocks'] = Nblocks
        dset.attrs['ave'] = float(temp[0])
        dset.attrs['err'] = float(temp[1])
      for line in open('results/plaq.suscept'):
        temp = line.split()
        dset.attrs[temp[0]] = float(temp[1])
        dset.attrs[temp[0] + '_err'] = float(temp[2])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
          print("%s vs %d in plaq.suscept" % (temp[-1], Nblocks))
          sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Next do simple observables measured every trajectory=MDTU
      #   pbp, exp_dS
      #   Skip xpoly, wall_time
      # Dynamically figure out how many data on each line
      # TODO: Add others... pbp header problematic...
      for obs in ['pbp', 'exp_dS']:
        obsfile = 'data/' + obs + '.csv'
        traj = 0
        for line in open(obsfile):
          temp = line.split(',')
          if line.startswith('M') or line.startswith('t'):
            Ndat = len(temp) - 1       # Skipping MDTU label
            if Ndat == 1:
              dset = this_grp.create_dataset(obs, (Ntraj,), dtype='f')
            elif Ndat > 1:
              dset = this_grp.create_dataset(obs, (Ntraj,Ndat,), dtype='f')
            else:
              print("ERROR: Ndat=%d for %s" % (Ndat, obs))
            continue
          for j in range(Ndat):
            dset[traj][j] = float(temp[j + 1])
          traj += 1

        resfile = 'results/' + obs + '.dat'
        for line in open(resfile):
          if line.startswith('#'):
            continue
          temp = line.split()
          dset.attrs['ave'] = float(temp[0])
          dset.attrs['err'] = float(temp[1])
          if not int(temp[-1]) == Nblocks:
            print("ERROR: Nblocks mismatch in %s, " % this_str, end='')
            print("%s vs %d in %s" % (temp[-1], Nblocks, resfile))
            sys.exit(1)


        # TODO: Attributes and results... check Nblocks
      # ------------------------------------------------------------

      # ------------------------------------------------------------
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
      this_str.create_dataset('Wmeas', data=Wmeas)
      this_str.create_dataset('topo', data=topo)

      # Now all other observables measured every trajectory=MDTU
      # Try to figure out how many data to include from each line
      # TODO: Add others... (different c for Wflow vs. Wpoly...)

    # TODO: Record thermalization cuts for each volume
    for line in open('therm.sh'):
      temp = line.split(',')
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
