import math
import numpy as np

from space_constants import *
import spacetime_utils as stu


#     -----------------------------------------------------------------
#
#                              Ex3_12
#
#  this file demonstrates example 3-12.
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


year = 2001
mon = 3
day = 18
hr  = 12
min = 14
sec = 0.0
jd, jdfrac = stu.jday(year, mon, day, hr, min, sec)
print('input jd %14.7f %11.8f  \n\n' % (jd, jdfrac))

# check jdfrac for multiple days
if abs(jdfrac) >= 1.0:
    jd = jd + np.floor(jdfrac)
    jdfrac = jdfrac - np.floor(jdfrac)

# check for fraction of a day included in the jd
dt = jd - np.floor(jd) - 0.5
if (abs(dt) > 0.00000001):
    jd = jd - dt
    jdfrac = jdfrac + dt

# ----------------- find year and days of the year ---------------
temp   = jd - 2415019.5
tu     = temp / 365.25
year   = 1900 + np.floor(tu)
leapyrs= np.floor((year - 1901) * 0.25)
days   = np.floor(temp - ((year - 1900) * 365.0 + leapyrs))

# ------------ check for case of beginning of a year -------------
if days + jdfrac < 1.0:
    year   = year - 1
    leapyrs= np.floor((year - 1901) * 0.25)
    days   = np.floor(temp - ((year - 1900) *365.0 + leapyrs))

print('year %6i  \n' % year)
print('mon  %3i  \n' % mon)
print('days  %3i  \n\n' % days)

daysf = days + hr/24.0 + min/1440.0 + sec/86400.0
print('days  %11.7f  \n' % daysf)

utsec = (daysf - days) * 86400.0

temp   = utsec / 3600.0
hr  = np.fix(temp )
min = np.fix((temp - hr) * 60.0)
sec = (temp - hr - min/60.0) * 3600.0

print(' hr min sec %4i  %4i  %8.6f \n' % (hr, min, sec))


