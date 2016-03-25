import numpy as np
# ------------------------------------------------------------------
# Houdayer--Hartmann curve collapse analysis
# Requires that data has already been parsed into all_L, mf, MH and err
def FSS_HH(data, params):
  all_L, mf, MH, err = data
  ymin, ymax, ystep, y0, cut, c0 = params

  Npts = int((ymax - ymin) / ystep) + 1
  all_ym = np.linspace(ymin, ymax, Npts)
  S = np.zeros(Npts, dtype = np.float)
  Nover = np.zeros(Npts, dtype = np.int)

  # Maximum possible number of overlapping points
  maxPts = 2 * (len(all_L)  - 1)
  xL = np.zeros(maxPts, dtype = np.float)
  yL = np.zeros(maxPts, dtype = np.float)
  wL = np.zeros(maxPts, dtype = np.float)   # Weights from errors on yL

  # Scan over gamma_m
  for index in range(Npts):
    ym = all_ym[index]
    omega = -y0 / ym    # Need to recalculate scale and scalep
    # Cycle over each point j in each data set i
    for i_str in all_L:
      i = all_L.index(i_str)
      L = float(i_str)
      for j in range(len(mf[i])):
#        print "L[%d]=%d, m[%d][%d]=%.2g" % (i, L, i, j, mf[i][j])
        xij = L * np.power(mf[i][j], 1 / ym)
        scale = 1 + c0 * np.power(mf[i][j], omega)
        if xij < cut:
          continue            # Skip small-mass region on small volumes
        yij = L * MH[i][j] / scale
        yijErr = L * err[i][j] / scale

        # Find and save bracketing points in all other data sets
        # First clear previous round
        count = 0;          cntP1 = 1
        xL[:] = 0.;         yL[:] = 0.;           wL[:] = 0.
        for ip_str in all_L:
          ip = all_L.index(ip_str)
          Lp = float(ip_str)
          if ip == i:   # Looking only at different L
            continue
          for jp in range(len(mf[ip]) - 1):
            thisX = Lp * np.power(mf[ip][jp], 1 / ym)
            nextX = Lp * np.power(mf[ip][jp + 1], 1 / ym)
            if thisX <= xij and xij <= nextX:
#              print "  L[%d]=%d, m[%d][%d]=%.2g" % (ip, Lp, ip, jp, mf[ip][jp])
#              print "  %.2g <= %.2g <= %.2g" % (thisX, xij, nextX)
              scale0 = 1 + c0 * np.power(mf[ip][jp], omega)
              scale1 = 1 + c0 * np.power(mf[ip][jp + 1], omega)
              xL[count] = thisX
              xL[cntP1] = nextX
              yL[count] = Lp * MH[ip][jp] / scale0
              yL[cntP1] = Lp * MH[ip][jp + 1] / scale1
              wL[count] = np.power(Lp * err[ip][jp] / scale0, -2)
              wL[cntP1] = np.power(Lp * err[ip][jp + 1] / scale1, -2)
              count += 2
              cntP1 += 2

              # Avoid rare case that nextX = xij and would be counted twice
              break

        # Use bracketing points to construct fiducial curve F
        # Convert lists to numpy arrays for easy multiplication
        if count > 0:
          Nover[index] += 1     # We have overlap
          Ktot = np.sum(wL)     # Total weights
          Kx = np.sum(wL * xL)
          Ky = np.sum(wL * yL)
          Kxx = np.sum(wL * xL * xL)
          Kxy = np.sum(wL * xL * yL)
          Delta = Ktot * Kxx - Kx * Kx

          # Fiducial curve and uncertainty
          F = Kxx * Ky - Kx * Kxy + xij * (Ktot * Kxy - Kx * Ky)
          F /= Delta
          dFSq = (Kxx - 2 * xij * Kx + xij * xij * Ktot) / Delta

          # Deviation weighted by uncertainties
          S[index] += np.power(yij - F, 2) / (dFSq + np.power(yijErr, 2))

    # Scale by number of overlapping points
    # Set large dummy deviation if there are no overlapping points
    # Hopefully that won't be a common problem
    if Nover[index] == 0:
      print "Warning: no overlap for y_m=%g" % ym
      S[index] = 1e6
    else:
      S[index] /= Nover[index]
#    print "S=%.4g at y_m=%g" % (S[index], ym)

  return S, Nover
# ------------------------------------------------------------------



# ------------------------------------------------------------------
# Find range of ym for which S - S_min < 1
# This is interpreted as a 'one-sigma' confidence interval
def find_range(S, ymin, ymax, ystep):
  target = np.amin(S) + 1
  start = np.argmin(S)
  Npts = S.shape[0]

  lo = ymin   # Lowest possible
  hi = ymax   # Highest possible
  for i in range(start):
    j = start - i
    if S[j] > target:
      lo = ymin + j * ystep
      break
  for i in range(start, Npts):
    if S[i] > target:
      hi = ymin + i * ystep
      break

  return lo, hi
# ------------------------------------------------------------------
