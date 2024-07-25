#     -----------------------------------------------------------------
#
#                              Ex7_1
#
#  this file demonstrates example 7-1.
#
#                          companion code for
#             fundamentals of astrodynamics and applications
#                                 2007
#                            by david vallado
#
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#             7 jun 07  david vallado
#                         original
#  changes :
#            13 feb 07  david vallado
#                         original baseline
#
#     *****************************************************************


import numpy as np
import spacetime_utils as stu
import orbit_utils as obu
from space_constants import *
import space_conversions as sc

# ---- site
print('\n-------- site test \n')
latgd = 39.007 * deg2rad
lon = -104.883 * deg2rad
alt = 2.187 #km
rs, vs = obu.site(latgd, lon, alt)
print('site %14.7f%14.7f%14.7f%14.7f%14.7f%14.7f\n'
      % (rs[0], rs[1], rs[2], vs[0], vs[1], vs[2]))
print('--------------- razel tests ----------------------------\n')
rho = 604.68 #km
az = 205.6 * deg2rad
el = 30.7 * deg2rad
drho = 2.08 #km/s
daz = 0.15 * deg2rad
del_ = 0.17 * deg2rad
year = 1995
mon = 5
day = 20
hr = 3
min = 17
sec = 2.0
dut1 = 0.0
dat = 29
xp = 0.0
yp = 0.0
lod = 0.0
timezone = 0
terms = 0
ddpsi = 0.0
ddeps = 0.0
utc = sec
ut1 = utc + dut1
tai = utc + dat
tt = tai + 32.184
jdut1, jdut1frac = stu.jday(year, mon, day, hr, min, ut1)
jdtt, jdttfrac = stu.jday(year, mon, day, hr, min, tt)
ttt = (jdtt - 2451545.0) / 36525.0
print('year %5i ' % (year))
print('mon %4i ' % (mon))
print('day %3i ' % (day))
print('hr %3i:%2i:%8.6f\n' % (hr, min, sec))
print('dut1 %8.6f s' % (dut1))
print(' dat %3i s' % (dat))
print(' xp %8.6f "' % (xp))
print(' yp %8.6f "' % (yp))
print(' lod %8.6f s\n' % (lod))
print('           range km        az deg      el    deg     rngrt km/s      '
      'azrate deg/s  elrate deg/s')
print('rvraz %14.7f%14.7f%14.7f%14.7f%14.7f%14.7f\n'
      % (rho, az * rad2deg, el * rad2deg, drho, daz * rad2deg, del_ * rad2deg))
reci, veci = sc.razel2rv(rho, az, el, drho, daz, del_, latgd, lon, alt, ttt,
                        jdut1 + jdut1frac, lod, xp, yp, terms, ddpsi, ddeps)
print('r  \n', (reci))
print('v  \n', (veci))
rho, az, el, drho, daz, del_ = sc.rv2razel(reci, veci, latgd, lon, alt, ttt,
                                      jdut1 + jdut1frac, lod, xp, yp, terms, ddpsi,
                                      ddeps)
print('\n         rho          az                el            drho          daz         del_\n')
print('rvraz %14.7f%14.7f%14.7f%14.7f%14.7f%14.7f\n'
      % (rho, az * rad2deg, el * rad2deg, drho, daz * rad2deg, del_ * rad2deg))

p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = sc.rv2coe(reci, veci)

print('          p km       a km        ecc      incl deg     raan deg     '
      'argp deg      nu deg      m deg      arglat\n')
print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f\n'
      % (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg,
         nu * rad2deg, m * rad2deg, arglat * rad2deg))

























