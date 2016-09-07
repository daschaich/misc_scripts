#!/usr/bin/python
import sys
import numpy as np
# ------------------------------------------------------------------
# Compute tree-level finite-volume perturbative correction to t^2 E
# from Eq. 6.2 of arXiv:1406.0827, given L and c=sqrt(8t)/L

# Parse arguments: first is L, second in c
if len(sys.argv) < 3:
  print "Usage:", str(sys.argv[0]), "<L> <c>"
  sys.exit(1)
L = int(sys.argv[1])
c = float(sys.argv[2])
t = L**2 * c**2 / 8.0

# For later convenience
twopiOvL = 2.0 * np.pi / float(L)
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Integrate over (n1, n2, n3, n4), each from 0 to L-1,
# except for (0, 0, 0, 0), which is treated separately
# Be lazy and re-compute (almost) everything within the lowest-level loop
tSqE = 0.0
for n1 in range(L):
  for n2 in range(L):
    for n3 in range(L):
      for n4 in range(L):
        if n1 + n2 + n3 + n4 == 0:    # Zero-mode contribution is separate
          continue

        # Define p, phat and ptwiddle
        n_mu = np.array([n1, n2, n3, n4], dtype = np.float)
        p_mu = twopiOvL * n_mu
        phat_mu = 2.0 * np.sin(p_mu / 2.0)
        ptw_mu = np.sin(p_mu)

        phatSq = (phat_mu**2).sum()
        ptwSq = (ptw_mu**2).sum()

        # Tr[S^e] = sum_mu [ptw^2 - ptw_mu**2] [cos(p_mu / 2)**2]
        TrS = ((ptwSq - ptw_mu**2) * (np.cos(p_mu / 2))**2).sum()

        # Sum up exp(-2t phatSq) * TrS / phatSq
        tSqE += np.exp(-2.0 * t * phatSq) * TrS / phatSq

# Add zero-mode contribution
tSqE += 2.0

# Done after overall factor of 64pi^2 t^2 / (3 L^4)
tSqE *= (64.0 * np.pi**2 * t**2) / (3.0 * float(L**4))
print "%d %.4g %.6g %.8g" % (L, c, t, tSqE)

# Match format in Anna's Ct_pert_c0.38
# TODO: How is pert related to tSqE and zeromode?
#pert = ??? * (64.0 * np.pi**2 * t**2) / (3.0 * float(L**4))
#zeromode = 128 * np.pi**2 * t**2 / (3.0 * float(L**4))
#print "Lat %d c %.4g tc %.6g" % (L, c, t),
#print "C(t) %.8g" % tSqE,
#print "Pert %.8g" % pert,
#print "Zeromode %.8g" % zeromode
# ------------------------------------------------------------------
