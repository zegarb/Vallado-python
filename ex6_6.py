#     -----------------------------------------------------------------
#
#                              ex6_6.m
#
#  this file demonstrates example 6-6.
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
from space_constants import au, rad2deg, deg2rad, re, velkmps, tumin, mu
import space_conversions as sc
import spacemath_utils as smu

print('-------------------- problem ex 6-6 \n')
iinit = 55.0 * deg2rad
ifinal = 40.0 * deg2rad
ecc = 0.0
deltaomega = 45.0 * deg2rad
vinit = 5.892311 / 7.905365719
fpa = 0.0 * deg2rad
deltav = obu.iandnode(iinit, deltaomega, ifinal, vinit, fpa)
print('inclination and node changes \n')
print(' deltav %11.7f %11.7f km/s \n' % (deltav, deltav * 7.905365719))
# End of Book Example


r1, v1 = sc.coe2rv(11480.649, 0.0, 55.0 * deg2rad, 45.0 * deg2rad,
                  0.0 * deg2rad, 330.0 * deg2rad, 128.9041397 * deg2rad,
                  0.0, 0.0)
r2, v2 = sc.coe2rv(11480.649, 0.0, 40.0 * deg2rad, 90.0 * deg2rad,
                  0.0 * deg2rad, 330.0 * deg2rad, 97.3803453 * deg2rad,
                  0.0, 0.0)
vx= v2-v1
print(' form Dv vectors and look at magnitude difference \n comp %11.7f %11.7f %11.7f km/s in icrf %11.7f \n' % (vx[0], vx[1], vx[2], smu.mag(vx)))
