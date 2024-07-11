import math
import numpy as np

from space_constants import *
import spacetime_utils as stu
from space_conversions import ecef2eciiau06, ecef2eci

#     -----------------------------------------------------------------
#
#                              ex3_14
#
#  this file demonstrates example 3-14.
#
#                          companion code for
#             fundamentals of astrodynamics and applications
#                                 2013
#                            by david vallado
#
#     (h)               email davallado@gmail.com
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#            16 feb 19  david vallado
#                         update for new constants
#  changes :
#            13 feb 07  david vallado
#                         original baseline
#
#     *****************************************************************

# LEO test
recef = np.array([[-1033.4793830],  [7901.2952754],  [6380.3565958]])
vecef = np.array([[-3.225636520],   [-2.872451450],  [5.531924446]])
aecef = np.array([[0.001],          [0.002],         [0.003]])

year=2004
mon = 4
day = 6
hr =  7
min= 51
sec= 28.386009

dut1 = -0.4399619  # sec
dat  = 32         # sec
xp   = -0.140682 * arcsec2rad  # " to rad
yp   =  0.333309 * arcsec2rad
lod  =  0.0015563
ddpsi = -0.052195 * arcsec2rad  # " to rad
ddeps = -0.003875 * arcsec2rad
ddx = -0.000205 * arcsec2rad  # " to rad
ddy = -0.000136 * arcsec2rad
order = 106
eqeterms = 2
timezone = 0
opt = 'c' # specify the iau00 cio approach


print('input data \n')
print(' year %5i ' % year)
print(' mon %4i ' % mon)
print(' day %3i ' % day)
print(' %3i:%2i:%8.6f\n' % (hr, min, sec))
print(' dut1 %8.6f s' % dut1)
print(' dat %3i s' % dat)
print(' xp %8.6f "' % (xp * rad2arcsec))
print(' yp %8.6f "' % (yp * rad2arcsec))
print(' lod %8.6f s\n' % lod)
print(' ddpsi %8.6f " ddeps  %8.6f\n' % (ddpsi* rad2arcsec, ddeps* rad2arcsec))
print(' ddx   %8.6f " ddy    %8.6f\n' % (ddx* rad2arcsec, ddy* rad2arcsec))
print(' order %3i  eqeterms %3i  opt %3s \n' % (order, eqeterms, opt))
print('units are km and km/s and km/s2\n')

# -------- convtime    - convert time from utc to all the others
print('convtime results\n')
ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, \
    jdtdb, jdtdbfrac = stu.convtime(year, mon, day, hr, min, sec, timezone,
                                    dut1, dat)
print('ut1 %8.6f tut1 %16.12f jdut1 %18.11f ' % (ut1, tut1, jdut1+jdut1frac))
h, m, s = stu.sec2hms(ut1)
print('hms %3i %3i %8.6f \n' % (h, m, s))
print('utc %8.6f ' % utc)
h, m, s = stu.sec2hms(utc)
print('hms %3i %3i %8.6f \n' % (h, m, s))
print('tai %8.6f' % tai)
h, m, s = stu.sec2hms(tai)
print('hms %3i %3i %8.6f \n' % (h, m, s))
print('tt  %8.6f ttt  %16.12f jdtt  %18.11f '% (tt, ttt, jdtt+jdttfrac))
h, m, s = stu.sec2hms(tt)
print('hms %3i %3i %8.6f \n' % (h, m, s))
print('tdb %8.6f ttdb %16.12f jdtdb %18.11f\n' % (tdb, ttdb, jdtdb+jdtdbfrac))


recigg, vecigg, aecig = ecef2eciiau06(recef, vecef, aecef, ttt, jdut1+jdut1frac,
                                      lod, xp, yp, 'c', ddx, ddy )
print('GCRF          IAU-2006 CIO:')
print(' r: \n' , (recigg))
print(' v: \n' , (vecigg))
print(' a: \n' , (aecig))

