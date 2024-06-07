#     -----------------------------------------------------------------
#
#                              Ex5_6.m
#
#  this file demonstrates example 5-6.
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

import numpy as np
import orbit_utils as obu
import spacetime_utils as stu

jd, jdfrac = stu.jday(1995, 2, 15, 12, 0, 0.0)
r1 = np.array([0.0, -4464.696, -5102.509])

r2 = np.array([0.0, 5740.323, 3189.068])
los = obu.sight(r1, r2, 's')
print(los)
jd, jdfrac = stu.jday(1995, 2, 15, 0, 0, 0.0)
rsun, rtasc, decl = obu.sun(jd + jdfrac)
print('sun MOD %11.9f%11.9f%11.9f au\n' % (rsun[0], rsun[1], rsun[2]))
print('sun MOD %14.4f%14.4f%14.4f km\n' % (rsun[0] * 149597870.0,
                                           rsun[1] * 149597870.0,
                                           rsun[2] * 149597870.0))
los = obu.sight(r1, rsun * 149597870.0, 's')
print(los)
