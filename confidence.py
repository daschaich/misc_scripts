#!/usr/bin/python
import sys
from scipy import special
# ------------------------------------------------------------------
# Confidence level corresponding to given chiSq and dof
# This is just the (complement of the) incomplete Gamma function
#   Q(a, x) = 1 - P(a, x) = [1 / Gamma(a)] int_x^{\infty} e^{-t} t^{a - 1} dt
# where a = dof / 2 > 0 and x = chiSq / 2 >= 0
# SciPy's special.gammainc is
#                 P(a, x) = [1 / Gamma(a)] int_0^x e^{-t} t^{a - 1} dt
# Reproduces Numerical Recipes-based C code from Tom DeGrand

# Parse arguments: first is chiSq, second is dof
if len(sys.argv) < 2:
  print "Usage:", str(sys.argv[0]), "<chiSq> <dof>"
  sys.exit(1)
chiSq = float(sys.argv[1])
dof = int(sys.argv[2])

if dof <= 0:
  print "ERROR: dof must be positive!"
  sys.exit(1)
if chiSq < 0:
  print "ERROR: chiSq must be non-negative!"
  sys.exit(1)

print "%.4g" % (1.0 - special.gammainc(0.5 * dof, 0.5 * chiSq))
# ------------------------------------------------------------------
