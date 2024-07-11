import math
import numpy as np

from space_constants import *
import spacetime_utils as stu


#     -----------------------------------------------------------------
#
#                              Ex3_6
#
#  this file demonstrates example 3-6.
#
#                          companion code for
#             fundamentals of astrodyanmics and applications
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


# -------- gstime       - greenwich sidereal time
year = 1992
mon = 8
day = 20
hr = 12
min = 14
sec = 0.0

jdut1, jdut1frac = stu.jday(year, mon, day, hr, min, sec)

print('jd %18.10f  %18.10f   %18.10f \n' % (jdut1, jdut1frac, jdut1+jdut1frac))

print('\n--------gstime0 test \n' )
print('input year %4i \n' % year)

gst0 = stu.gstime0(year)

print('gst0 = %14.8f deg\n' % (gst0 * rad2deg))


# -------- lstime       - local sidereal time
lon   = -104.000 * deg2rad
print('input lon = %11.7f ' % (lon * rad2deg))

lst, gst = stu.lstime(lon, jdut1+jdut1frac)

print('gst = %14.8f\n' % (gst * rad2deg))
print('lst %11.7f gst %11.7f deg\n' % (lst * rad2deg, gst * rad2deg))

