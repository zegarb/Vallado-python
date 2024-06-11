#     -----------------------------------------------------------------
#
#                              Ex10_5.py
#
#  this file demonstrates example 10-5.
#
#                          companion code for
#             fundamentals of astrodyanmics and applications
#                                 2007
#                            by david vallado
#
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#            12 oct 11  david vallado
#                         original
#  changes :
#             7 oct 11  david vallado
#                         original baseline
#
#     *****************************************************************

# problem 5
import numpy as np
from space_constants import *
from space_conversions import razel2rv, rv2coe
from spacetime_utils import jday, convtime
from orbit_utils import site, hgibbs, pkepler, findatwaatwb
from spacemath_utils import mag
import os

print('Example 10-5 -------------------------\n')
np.set_printoptions(precision=8, floatmode='fixed')

# read in GEOS data
fn = os.path.join(os.path.dirname(__file__), "data", "geos6a.dat")
filedat = np.loadtxt(fn)
numobs = filedat.shape[0]

#load data into x y z arrays
yeararr = filedat[:, 3]
monarr = filedat[:, 4]
dayarr = filedat[:, 5]
hrarr = filedat[:, 6]
minarr = filedat[:, 7]
secarr = filedat[:, 8]
rngarr = filedat[:, 9]

azarr = filedat[:, 10] * deg2rad

elarr = filedat[:, 11] * deg2rad

latgd = 21.572056 * deg2rad

lon = -158.266578 * deg2rad
alt = 0.3002

dat = 29

dut1 = 0.3261068

lod = 0.0
xp = -0.11554 * arcsec2rad

yp = 0.48187 * arcsec2rad
rs, vs = site(latgd, lon, alt)
obsrecarr = []
for j in range(numobs):
    obsrec = {}
    obsrec['time'], obsrec['timef'] = jday(yeararr[j], monarr[j], dayarr[j],
                                           hrarr[j], minarr[j], secarr[j])
    obsrec['latgd'] = latgd
    obsrec['lon'] = lon
    obsrec['alt'] = alt
    ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb,\
    jdtdb, jdtdbfrac = convtime(yeararr[j], monarr[j], dayarr[j],
                                hrarr[j], minarr[j], secarr[j], 0, dut1, dat)
    obsrec['ttt'] = ttt
    obsrec['jdut1'] = jdut1 + jdut1frac
    obsrec['xp'] = xp
    obsrec['yp'] = yp
    obsrec['noiserng'] = 0.0925
    obsrec['noiseaz'] = 0.0224 * deg2rad
    obsrec['noiseel'] = 0.0139 * deg2rad
    obsrec['obstype'] = 2
    obsrec['rsecef'] = rs
    obsrec['rng'] = rngarr[j] #- 0.08
    obsrec['az'] = azarr[j] #- 0.0081 * deg2rad
    obsrec['el'] = elarr[j] #- 0.0045 * deg2rad
    obsrecarr.append(obsrec)


firstobs = 0
lastobs = 9
# -------  form nominal vector
reci = np.array([0, 0, 0])
veci = np.array([0, 0, 0])
for obsktr in range(firstobs + 1, lastobs - 1):
    currobsrec = obsrecarr[obsktr - 1]
    re1, ve1 = razel2rv(currobsrec['rng'], currobsrec['az'], currobsrec['el'], 0.0, 0.0, 0.0,
                       currobsrec['latgd'], currobsrec['lon'], currobsrec['alt'],
                       currobsrec['ttt'], currobsrec['jdut1'], 0.0, currobsrec['xp'],
                       currobsrec['yp'], 2, 0.0, 0.0)
    currobsrec = obsrecarr[obsktr]
    re2, ve2 = razel2rv(currobsrec['rng'], currobsrec['az'], currobsrec['el'], 0.0, 0.0, 0.0,
                       currobsrec['latgd'], currobsrec['lon'], currobsrec['alt'],
                       currobsrec['ttt'], currobsrec['jdut1'], 0.0, currobsrec['xp'],
                       currobsrec['yp'], 2, 0.0, 0.0)
    currobsrec = obsrecarr[obsktr + 1]
    re3, ve3 = razel2rv(currobsrec['rng'], currobsrec['az'], currobsrec['el'], 0.0, 0.0, 0.0,
                       currobsrec['latgd'], currobsrec['lon'], currobsrec['alt'],
                       currobsrec['ttt'], currobsrec['jdut1'], 0.0, currobsrec['xp'],
                       currobsrec['yp'], 2, 0.0, 0.0)
    #        [ve2, theta, theta1, copa, error] = gibbs( re1, re2, re3);
    ve2, theta, theta1, copa, error = \
        hgibbs(re1, re2, re3, obsrecarr[obsktr - 1]['jdut1'],
               obsrecarr[obsktr]['jdut1'], obsrecarr[obsktr + 1]['jdut1'])
    # move back to 1st time, not central time
    dtsec = (obsrecarr[obsktr]['time'] + obsrecarr[obsktr]['timef']
             - obsrecarr[1]['time'] - obsrecarr[1]['timef']) * 86400.0
    #        [reci1, veci1] =  kepler ( re2, ve2, -dtsec );
    reci1, veci1 = pkepler(re2, ve2, - dtsec, 0.0, 0.0)
    print(f're {reci1} km ')
    print(f've {veci1} km {dtsec:11.7f} dtsec\n')
    reci = reci + reci1
    veci = veci + veci1

reci = reci / (lastobs - firstobs - 1)
veci = veci / (lastobs - firstobs - 1)
print(f'rnom {reci} km ')
print(f'vnom {veci} km \n')
jdepoch = obsrecarr[0]['time'] + obsrecarr[0]['timef']
# use a vector that's farther off to see iteratins take effect
reci = np.array([5975.2904, 2568.64, 3120.5845])
veci = np.array([3.983846, -2.071159, -5.917095])
print('input: \n')
print('rnom %16.8f %16.8f %16.8f km ' % (reci[0], reci[1], reci[2]))
print('vnom %16.8f %16.8f %16.8f km \n' % (veci[0], veci[1], veci[2]))
xnom = np.zeros((6, 1))
xnom[0, 0] = reci[0]
xnom[1, 0] = reci[1]
xnom[2, 0] = reci[2]
xnom[3, 0] = veci[0]
xnom[4, 0] = veci[1]
xnom[5, 0] = veci[2]
r1 = np.zeros(3)
v1 = np.zeros(3)
for j in range(5):
    # ---- accumulate obs and assemble matrices
    percentchg = 0.01
    deltaamtchg = 0.01
    atwa, atwb, atw, b, drng2, daz2, del2 = \
        findatwaatwb(firstobs, lastobs, obsrecarr, 6, percentchg, deltaamtchg, xnom)
    if j == 0:
        print('atwa = \n')
        print(f'{atwa}')
        print('atwb = \n')
        print(f'{atwb}')
    deltax = np.linalg.inv(atwa) @ atwb
    print(f'dx{deltax} \n')
    xnom = xnom + deltax
    r1[0] = xnom[0, 0]
    r1[1] = xnom[1, 0]
    r1[2] = xnom[2, 0]
    v1[0] = xnom[3, 0]
    v1[1] = xnom[4, 0]
    v1[2] = xnom[5, 0]
    # answer in km and km/s
    print('output: \n')
    print(f'r1 {r1} km ')
    print(f'v1 {v1} km \n')
    if 2 == (obsrecarr[1]['obstype']):
        sigmanew = np.sqrt((drng2 + daz2 + del2) / (lastobs - firstobs))
        #            3 : sigmanew= SQRT( (DRng2 + DAz2 + DEl2 + DDRng2) / NumWork );
#            4 : sigmanew= SQRT( (DRng2 + DAz2 + DEl2 + DDRng2 + DDAz2 + DDEl2)
#                                 / NumWork );
#            5 : sigmanew= SQRT( (DTRtAsc2 + DTDecl2) / NumWork );
#            6 : sigmanew= SQRT( DRng2 / NumWork );
    print('rms  %16.8f  \n' % (sigmanew))

p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(r1, v1)
print('          p km       a km      ecc      incl deg     raan deg     '
      'argp deg      nu deg      m deg      arglat   truelon    lonper\n')
print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n'
      % (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg,
         nu * rad2deg, m * rad2deg))
if arglat:
    print(arglat * rad2deg)
elif truelon:
    print(truelon * rad2deg)
elif lonper:
    print(lonper * rad2deg)

rans = np.array([5748.607, 2679.9324, 3442.7902])

vans = np.array([4.33046, - 1.922874, - 5.726562])
r2shans = np.array([5746.4693, 2679.944362, 3442.41616])

v2shans = np.array([4.33993195, - 1.91866466, - 5.72186789])
r3ans = np.array([5748.398438, 2680.135601, 3442.296315])

v3ans = np.array([4.32743417, - 1.9201418, - 5.72605888])

dr = r2shans - r1
print(f'diff odtk short bad init {dr} {mag(dr):16.8f}  km \n')
dr = rans - r1
print(f'diff gtds w bad init      {dr} {mag(dr):16.8f}   km \n')
dr = r3ans - r1
print(f'diff odtk long w bad init {dr} {mag(dr):16.8f}  km \n')
p = np.linalg.inv(atwa)
print('covariance p = \n')
print(f'{p} \n')
print('cov r, 1 sigma  %16.8f  %16.8f  %16.8f  km \n'
      % (np.sqrt(p[0, 0]) * 1000, np.sqrt(p[1, 1]) * 1000, np.sqrt(p[2, 2]) * 1000))
d, eigenaxes = np.linalg.eig(p)
eigenvalues = np.sqrt(d) * 1000


# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(1, 6+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(7, 12+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(13, 18+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(19, 24+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(25, 30+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(31, 36+1))))

print(f'{eigenaxes.T}')
print('eigenvalues   %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  m \n'
      % (eigenvalues[0], eigenvalues[1], eigenvalues[2], eigenvalues[3], eigenvalues[4], eigenvalues[5]))
print('mag eigenvalues  %16.8f m  \n'
      % (np.sqrt(eigenvalues[0] ** 2
                 + eigenvalues[1] ** 2 + eigenvalues[2] ** 2)))
print('mag eigenvalues  %16.8f m/s \n'
      % (np.sqrt(eigenvalues[3] ** 2
                 + eigenvalues[4] ** 2 + eigenvalues[5] ** 2)))
rms = np.sqrt(b * np.transpose(b) / 4)
rmsold = rms
#     fprintf(1, '0 dx #11.7f  #11.7f ans #11.7f  #11.7f rms #11.7f \n', ans(1), ans(2), alpha, beta, rms);

# --------------------------------------- now do sequential part  -----------------------------------
firstobs = 10
lastobs = 17
atwaold = atwa
atwbold = atwb
numobs = 8
jdepoch = obsrecarr[0]['time'] + obsrecarr[0]['timef']
# use vector from first 10 points
#   reci = [5975.2904  2568.6400  3120.5845];
#   veci = [3.983846  -2.071159  -5.917095];

xnom[0, 0] = r1[0]
xnom[1, 0] = r1[1]
xnom[2, 0] = r1[2]
xnom[3, 0] = v1[0]
xnom[4, 0] = v1[1]
xnom[5, 0] = v1[2]
print('input: \n')
print(f'rnom {r1} km ')
print(f'vnom {v1} km \n')
for j in range(3):
    # ---- accumulate obs and assemble matrices
    percentchg = 0.01
    deltaamtchg = 0.01
    atwa, atwb, atw, b, drng2, daz2, del2 = findatwaatwb(firstobs, lastobs, obsrecarr, 6, percentchg, deltaamtchg, xnom)
    if j == 1:
        print('atwa = \n')
        print(f'{atwa}\n')
        #            fprintf(1, '#10.1f #10.1f #10.1f #10.1f #10.1f #10.1f\n', atwa+atwaold );
        print('atwb = \n')
        print(f'{atwb}\n')
    deltax = np.linalg.inv(atwa + atwaold) @ (atwb + atwbold)
    print(f'dx {deltax} \n')
    xnom = xnom + deltax
    r1[0] = xnom[0, 0]
    r1[1] = xnom[1, 0]
    r1[2] = xnom[2, 0]
    v1[0] = xnom[3, 0]
    v1[1] = xnom[4, 0]
    v1[2] = xnom[5, 0]
    # answer in km and km/s
    print('output: \n')
    print(f'r1 {r1} km ')
    print(f'v1 {v1} km \n')
    if 2 == (obsrecarr[1]['obstype']):
        sigmanew = np.sqrt((drng2 + daz2 + del2) / (lastobs - firstobs))
        #            3 : sigmanew= SQRT( (DRng2 + DAz2 + DEl2 + DDRng2) / NumWork );
#            4 : sigmanew= SQRT( (DRng2 + DAz2 + DEl2 + DDRng2 + DDAz2 + DDEl2)
#                                 / NumWork );
#            5 : sigmanew= SQRT( (DTRtAsc2 + DTDecl2) / NumWork );
#            6 : sigmanew= SQRT( DRng2 / NumWork );
    print('rms  %16.8f  \n' % (sigmanew))

p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(r1, v1)
print('          p km       a km      ecc      incl deg     raan deg     '
      'argp deg      nu deg      m deg      arglat   truelon    lonper\n')
print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n'
      % (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg,
         nu * rad2deg, m * rad2deg))
if arglat:
    print(arglat * rad2deg)
elif truelon:
    print(truelon * rad2deg)
elif lonper:
    print(lonper * rad2deg)

rans = np.array([5748.607, 2679.9324, 3442.7902])

vans = np.array([4.33046, - 1.922874, - 5.726562])
r2shans = np.array([5746.4693, 2679.944362, 3442.41616])

v2shans = np.array([4.33993195, - 1.91866466, - 5.72186789])
r3ans = np.array([5748.398438, 2680.135601, 3442.296315])

v3ans = np.array([4.32743417, - 1.9201418, - 5.72605888])

dr = r2shans - r1
print(f'diff odtk short bad init {dr} {mag(dr):16.8f}  km \n')
dr = rans - r1
print(f'diff gtds w bad init      {dr} {mag(dr):16.8f}   km \n')
dr = r3ans - r1
print(f'diff odtk long w bad init {dr} {mag(dr):16.8f}  km \n')
p = np.linalg.inv(atwa)
print('covariance p = \n')
print(f'{p} \n')
print('cov r, 1 sigma  %16.8f  %16.8f  %16.8f  km \n'
      % (np.sqrt(p[0, 0]) * 1000, np.sqrt(p[1, 1]) * 1000, np.sqrt(p[2, 2]) * 1000))
d, eigenaxes = np.linalg.eig(p)
eigenvalues = np.sqrt(d) * 1000

# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(1, 6+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(7, 12+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(13, 18+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(19, 24+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(25, 30+1))))
# print('eigenaxes  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  km^2 \n'
#       % (eigenaxes(np.arange(31, 36+1))))

print(f'{eigenaxes.T}')
print('eigenvalues   %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  %16.8f  m \n'
      % (eigenvalues[0], eigenvalues[1], eigenvalues[2], eigenvalues[3],
          eigenvalues[4], eigenvalues[5]))
print('mag eigenvalues  %16.8f m  \n'
      % (np.sqrt(eigenvalues[0] ** 2
                 + eigenvalues[1] ** 2 + eigenvalues[2] ** 2)))
print('mag eigenvalues  %16.8f m/s \n'
      % (np.sqrt(eigenvalues[3] ** 2
                 + eigenvalues[4] ** 2 + eigenvalues[5] ** 2)))
rms = np.sqrt(b * np.transpose(b) / 4)
rmsold = rms
#     fprintf(1, '0 dx #11.7f  #11.7f ans #11.7f  #11.7f rms #11.7f \n', ans(1), ans(2), alpha, beta, rms);