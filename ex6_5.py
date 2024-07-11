#     -----------------------------------------------------------------
#
#                              ex6_5
#
#  this file demonstrates example 6-5.
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

print('-------------------- problem ex 6-5 \n')

ecc = 0.0
deltaomega = 45.0 * deg2rad
vinit = 5.892311 / 7.905365719
fpa = 0.0 * deg2rad
incl = 55.0 * deg2rad
ifinal, deltav = obu.nodeonly(ecc, deltaomega, vinit, incl, fpa)
print('node only changes \n')
print(' ifinal %11.7f \n' % (ifinal * rad2deg))
print(' deltav %11.7f %11.7f \n' % (deltav, deltav * 7.905365719))
# End of book example


#                   node       argp      nu            ci           ce   ee
r1, v1 = sc.coe2rv(11480.649, 0.0, 55.0 * deg2rad, 45.0 * deg2rad,
                  30.0 * deg2rad, 330.0 * deg2rad, 103.3647275 * deg2rad, 0.0,
                  0.0)
r2, v2 = sc.coe2rv(11480.649, 0.0, 55.0 * deg2rad, 90.0 * deg2rad,
                  30.0 * deg2rad, 330.0 * deg2rad, 76.6352725 * deg2rad, 0.0,
                  0.0)
vx = v2 - v1
print(' comp %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vx[0], vx[1], vx[2], smu.mag(vx)))

ecc = 0.0
deltaomega = 0.98564736 * deg2rad

vinit = np.sqrt(398600.4418 / 42164) / 7.905365719

fpa = 0.0 * deg2rad
incl = 90.0 * deg2rad
ifinal, deltav = obu.nodeonly(ecc, deltaomega, vinit, incl, fpa)
print('node only changes \n')
print(' ifinal %11.7f \n' % (ifinal * rad2deg))
print(' deltav %11.7f %11.7f km/s \n' % (deltav, deltav * 7.905365719))
#                  node                argp      nu          ci    ce   ee
r1, v1 = sc.coe2rv(42164, 0.0, 90.0 * deg2rad, 0.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 90.0 * deg2rad, 0.0, 0.0)
r2, v2 = sc.coe2rv(42164, 0.0, 90.0 * deg2rad, 0.0 * deg2rad + deltaomega,
                  0.0 * deg2rad, 0.0 * deg2rad, 90.0 * deg2rad, 0.0, 0.0)
vx = v2 - v1
print(' manv %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vx[0], vx[1], vx[2], smu.mag(vx)))
