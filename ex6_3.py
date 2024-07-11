#     -----------------------------------------------------------------
#
#                              ex6_3
#
#  this file demonstrates example 6-3.
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
from space_constants import au, rad2deg, deg2rad, re, velkmps, tumin
import space_conversions as sc
import spacemath_utils as smu

print('-------------------- problem ex 6-3 \n' % ())
rinit = (re + 191.3411) / re
rfinal = (re + 35781.34857) / re
einit = 0.0
efinal = 0.0
nuinit = 0.0 * deg2rad
nutran = 160.0 * deg2rad
print('initial position \n')
print(' rinit  %11.7f %11.7f km \n' % (rinit, rinit * re))
print(' rfinal %11.7f %11.7f km \n' % (rfinal, rfinal * re))
print(' einit   %11.7f \n' % (einit))
print(' efinal  %11.7f \n' % (efinal))
print(' nuinit  %11.7f deg \n' % (nuinit * rad2deg))
print(' nutran %11.7f deg \n' % (nutran * rad2deg))
deltava, deltavb, dttu, etran, atran, vtrana, vtranb = \
    obu.onetang(rinit, rfinal, einit, efinal, nuinit, nutran)
print('one tangent answers \n')
print(' deltava  %11.7f  %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f  %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltav  %11.7f %11.7f   km/s \n' % (deltavb + deltava,
                                            (deltava + deltavb) * velkmps))
print(' dttu  %11.7f tu %11.7f min \n' % (dttu, dttu * tumin))
print(' tran a %11.7f %11.7f  e %11.7f \n' % (atran, atran * re, etran))

# ellip equatorial
p = re * atran * (1.0 - etran * etran)
r1, v1 = sc.coe2rv(p, etran, 0.0 * deg2rad, 0.0 * deg2rad, 0.0 * deg2rad, nutran,
                   0.0 * deg2rad, 0.0, 0.0 * deg2rad)
# p = a for circular, cir equatorial
p = re + 35781.34857
r2, v2 = sc.coe2rv(p, efinal, 0.0 * deg2rad, 0.0 * deg2rad, 0.0 * deg2rad,
                   0.0 * deg2rad, 0.0 * deg2rad, nutran, 0.0)
vx = v2 - v1
print(' comp %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vx[0], vx[1], vx[2], smu.mag(vx)))
print(' r1  %11.7f %11.7f %11.7f  km v %11.7f \n'
      % (r1[0], r1[1], r1[2], smu.mag(v1)))
print(' r2  %11.7f %11.7f %11.7f  km v %11.7f \n'
      % (r2[0], r2[1], r2[2], smu.mag(v2)))
