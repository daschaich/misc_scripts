#!/usr/bin/python
import numpy as np
from scipy import special
# ------------------------------------------------------------------
# Test when confidence level underflows to zero
for chiSq in range(100):
  CL = 1.0 - special.gammainc(0.5, 0.5 * chiSq)
#  CL = np.longdouble(1.0 - special.gammainc(0.5, 0.5 * chiSq))
  print "%d %.4g" % (chiSq, CL)
# ------------------------------------------------------------------
