#     -----------------------------------------------------------------
#
#                              Ex5_3
#
#  this file demonstrates example 5-3.
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

import math
import numpy as np
import spacetime_utils as stu
from orbit_utils import moon
from space_constants import *

ddpsi = -0.052195 * arcsec2rad  # " to rad
ddeps = -0.003875 * arcsec2rad

# --------  moon -------------- analytical moon ephemeris
year = 1994

mon = 4
day = 27
hr = 23
minute = 58
second = 59.816
dat = 28
dut1 = -0.0889721
timezone = 0
ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, \
jdtdb, jdtdbfrac = stu.convtime(year, mon, day, hr, minute, second,
                                timezone, dut1, dat)
print('input data \n')
print(' year %5i ' % (year))
print(' mon %4i ' % (mon))
print(' day %3i ' % (day))
print(' %3i:%2i:%8.6f\n ' % (hr, minute, second))
print(' dut1 %8.6f s' % (dut1))
print(' dat %3i s' % (dat))
print(' ddpsi %8.6f " ddeps  %8.6f\n' % (ddpsi * rad2arcsec, ddeps * rad2arcsec))
print('tt  %8.6f ttt  %16.12f jdtt  %18.11f \n' % (tt, ttt, jdtt + jdttfrac))
print('tut1  %8.6f tut1  %16.12f jdut1  %18.11f \n'
      % (tut1, tut1, jdut1 + jdut1frac))
rmoon, rtasc, decl = moon(jdtt + jdttfrac)
print('rmoon  rtasc %14.6f deg decl %14.6f deg\n' % (rtasc * rad2deg, decl * rad2deg))
print('rmoon new      %11.9f%11.9f%11.9f er\n' % (rmoon[0], rmoon[1], rmoon[2]))
print('rmoon new   %14.4f%14.4f%14.4f km\n'
      % (rmoon[0] * re, rmoon[1] * re, rmoon[2] * re))
#  1994  4 28  119377864.5535      84237799.9341      36522905.6699     150602138.5030  0.9867078
#  1994  4 28    -134038.3192       -311589.2121       -126061.1912
rmoonaa = np.array([- 134038.3192, - 311589.2121, - 126061.1912])
print('moon jplde430 -134038.3192  -311589.2121  -126061.1912 km\n')
print(rmoonaa)



