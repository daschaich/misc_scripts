#!/usr/bin/python
import sys
import numpy as np
# ------------------------------------------------------------------
# Compute tree-level finite-volume perturbative correction to t^2 E
# from Eq. 6.2 of arXiv:1406.0827, given L, Nt and range of t
# Use plaquette discretization of E(t) (Eq. 3.5 with c=0):
#   Tr[S^e] = sum_mu [phat^2 - phat_mu**2]

# Parse arguments: first is L, second is Nt,
# third and fourth define range of t to cover
if len(sys.argv) < 5:
  print "Usage:", str(sys.argv[0]), "<L> <Nt> <t-start> <t-end>"
  sys.exit(1)
L = int(sys.argv[1])
Nt = int(sys.argv[2])
start = float(sys.argv[3])
end = float(sys.argv[4])

# For later convenience
one_ov_L = np.array([1.0 / float(L), 1.0 / float(L),
                     1.0 / float(L), 1.0 / float(Nt)], dtype = np.float)
twopiOvL = 2.0 * np.pi * one_ov_L
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# May need to manually check whether or not the endpoint is included
for t in np.arange(start, end, 0.01):
  if t <= 0:
    continue

  # Integrate over (n1, n2, n3, n4), each from 0 to L-1 or Nt-1,
  # except for (0, 0, 0, 0), which is treated separately
  # Be lazy and re-compute (almost) everything within the lowest-level loop
  tSqE = 0.0
  for n1 in range(L):
    for n2 in range(L):
      for n3 in range(L):
        for n4 in range(Nt):
          if n1 + n2 + n3 + n4 == 0:    # Zero-mode contribution is separate
            continue

          # Define p and phat
          n_mu = np.array([n1, n2, n3, n4], dtype = np.float)
          p_mu = twopiOvL * n_mu
          phat_mu = 2.0 * np.sin(p_mu / 2.0)
          phatSq = (phat_mu**2).sum()

          # Tr[S^e] = sum_mu [phat^2 - phat_mu**2]
          TrS = (phatSq - phat_mu**2).sum()

          # Sum up exp(-2t phatSq) * TrS / phatSq
          tSqE += np.exp(-2.0 * t * phatSq) * TrS / phatSq

  # Add zero-mode contribution
  tSqE += 2.0

  # Done after overall factor of 64pi^2 t^2 / (3 L^3 Nt)
  tSqE *= (64.0 * np.pi**2 * t**2) / (3.0 * float(L**3 * Nt))
  print "%g %.8g" % (t, tSqE)
# ------------------------------------------------------------------
