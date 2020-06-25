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
# Group paths will specify Nf, Nt, L, beta_F, (valence) mass and start
# Attributes for each stream:
#   Number of trajectories, acceptance rate, (Wpoly) autocorrelation time,
#   additional pbp autocorrelation time, thermalization cut, block size,
#   number of blocks
# Datasets for each stream will each with ave, err (suscept) as attributes
#   plaquette, chiral condensate, exp(-Delta S), TODO...
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
      os.chdir(path + Nf + '/' + vol + '/' + stream)
      toCheck = 'results/Wpoly_mod.autocorr'
      if not os.path.isfile(toCheck):     # Skip unfinished streams
        print("Skipping %s" % Nf + '/' + vol + '/' + stream)
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
      # Record acceptance, auto-correlation and Nblocks for stream
      for line in open('results/accP.dat'):
        temp = line.split()
        Nblocks = int(temp[-1])
        this_grp.attrs['Nblocks'] = Nblocks
        this_grp.attrs['acceptance'] = float(temp[0])
      for line in open('results/Wpoly_mod.autocorr'):
        temp = line.split()
        this_grp.attrs['autocorrelation time'] = float(temp[0])
      for line in open('results/pbp.autocorr'):
        if not line.startswith('W'):
          temp = line.split()
          this_grp.attrs['pbp autocorrelation time'] = float(temp[0])

      # Now set up data and attributes/results for this stream
      # First do plaq, setting Ntraj
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
      dset = this_grp.create_dataset('plaq', data=plaq)
      dset.attrs['columns'] = ['ss', 'st']

      # Record results as attributes
      for line in open('results/plaq.dat'):
        temp = line.split()
        dset.attrs['ave'] = float(temp[0])
        dset.attrs['err'] = float(temp[1])
        if not int(temp[-1]) == Nblocks:
          print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
          print("%s vs %d in plaq.suscept" % (temp[-1], Nblocks))
          sys.exit(1)
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
      # TODO: Add others...
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
          if Ndat == 1:
            dset[traj] = float(temp[1])
          else:
            for j in range(Ndat):
              dset[traj][j] = float(temp[j + 1])
          traj += 1

        # TODO: Attributes and results... check Nblocks
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

        # Susceptibility, skewness, kurtosis if present
        suscfile = 'results/' + obs + '.suscept'
        if not os.path.isfile(suscfile):     # No suscept for some obs
          continue
        for line in open(suscfile):
          temp = line.split()
          dset.attrs[temp[0]] = float(temp[1])
          dset.attrs[temp[0] + '_err'] = float(temp[2])
          if not int(temp[-1]) == Nblocks:
            print("ERROR: Nblocks mismatch in %s: " % this_str, end='')
            print("%s vs %d in %s.suscept" % (temp[-1], Nblocks, obs))
            sys.exit(1)
      # ------------------------------------------------------------

      # ------------------------------------------------------------
      # Now for Wilson-flowed observables
      # First set Nflow by going through topo.dat to set NWflow
      # Only have c=0.5 topological charge (for clean time-series plot)
      # Also save trajectories at which Wilson flow measurements were run
      topo_arr = []
      Wmeas_arr = []
      for line in open('data/topo.csv'):
        if line.startswith('M'):
          continue
        temp = line.split(',')
        Wmeas_arr.append(int(temp[0]))
        topo_arr.append(float(temp[4]))
      NWflow = len(Wmeas_arr)
      Wmeas = np.array(Wmeas_arr)
      topo = np.array(topo_arr)
      this_grp.create_dataset('Wmeas', data=Wmeas)
      this_grp.create_dataset('topo', data=topo)

      # TODO: Now all other Wilson-flowed observables
      # Try to figure out how many data to include from each line

    # TODO: Record thermalization cuts for each volume
    therm = path + Nf + '/' + vol + '/therm.sh'
    for line in open(therm):
      temp = line.split(',')
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Zero-temperature runs are a bit simpler
os.chdir(path + '4f/24nt48')
for ens in glob.glob('b[1-9]*'):
  temp = ens.split('_')
  beta = temp[0]
  mass = temp[1]
  this_grp = f.create_group('4f/24nt48/' + mass + '/' + beta)

  # Record acceptance and Nblocks for stream
  os.chdir(path + '4f/24nt48/' + ens)
  for line in open('results/accP.dat'):
    temp = line.split()
    Nblocks = int(temp[-1])
    this_grp.attrs['Nblocks'] = Nblocks
    this_grp.attrs['acceptance'] = float(temp[0])
# ------------------------------------------------------------------
