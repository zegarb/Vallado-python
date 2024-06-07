import math
import numpy as np

from space_constants import *
import spacetime_utils as stu


#     -----------------------------------------------------------------
#
#                              Ex3_4.m
#
#  this file demonstrates example 3-4.
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


# -------- jday         - find julian date
year = 1996
mon = 10
day = 26
hr  = 14
minute = 20
secs = 0.00
print('\n--------jday test \n' )
print('year %4i '%year)
print('mon %4i '%mon)
print('day %3i '%day)
print('hr %3i:%2i:%8.6f\n '%(hr, minute, secs ))

jd, jdfrac = stu.jday(year, mon, day, hr, minute, secs)
print('jd %18.10f  %18.10f   %18.10f \n\n\n'%(jd, jdfrac, jd + jdfrac))


# alt tests
year, mon, day, hr, minute, secs = stu.invjday ( 2450382.5, jdfrac )
print('year %5i   mon %4i day %3i %3i:%2i:%8.6f\n'%(year, mon, day, hr, minute, secs))

year, mon, day, hr, minute, secs = stu.invjday ( 2450382.5, jdfrac - 0.2 )
print('year %5i   mon %4i day %3i %3i:%2i:%8.6f\n'%(year, mon, day, hr, minute, secs))

year, mon, day, hr, minute, secs = stu.invjday ( 2450382.5 + 1, jdfrac + 1.5 )
print('year %5i   mon %4i day %3i %3i:%2i:%8.6f\n'%(year, mon, day, hr, minute, secs))

year, mon, day, hr, minute, secs = stu.invjday ( 2450382.5, -0.5 )
print('year %5i   mon %4i day %3i %3i:%2i:%8.6f\n'%(year, mon, day, hr, minute, secs))

year, mon, day, hr, minute, secs = stu.invjday ( 2450382.5, +0.5 )
print('year %5i   mon %4i day %3i %3i:%2i:%8.6f\n'%(year, mon, day, hr, minute, secs))
