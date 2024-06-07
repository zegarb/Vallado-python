#     -----------------------------------------------------------------
#
#                              ex6_4.m
#
#  this file demonstrates example 6-4.
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

print('-------------------- problem ex 6-4 \n')
deltai = 15.0 * deg2rad
vinit = 5.892311 / velkmps

incl = 28.5 * deg2rad
fpa = 0.0 * deg2rad
deltavionly = obu.ionlychg(deltai, vinit, fpa)
print('inclination only changes \n')
print(' vinit   %11.7f  %11.7f \n' % (vinit, vinit * velkmps))
print(' deltavionly  %11.7f  %11.7f \n' % (deltavionly, deltavionly * velkmps))
print(' comp  %11.7f  %11.7f \n' % (-deltavionly * velkmps * np.cos(deltai),
                                    deltavionly * velkmps * np.sin(deltai)))
print('-----part 2 \n')
deltai = 15.0 * deg2rad
ecc = 0.3
p = 17858.7836 / re
nu = 330.0 * deg2rad
mu = 1
a = p / (1.0 - ecc * ecc)
r = p / (1.0 + ecc * np.cos(nu))
vinit = np.sqrt((2.0 * mu) / r - (mu / a))
fpa = np.arctan((ecc * np.sin(nu)) / (1.0 + ecc * np.cos(nu)))
deltavionly = obu.ionlychg(deltai, vinit, fpa)
print(' a  %11.7f  %11.7f km r %11.7f  %11.7f km fpa %11.7f \n\n'
      % (a, a * re, r, r * re, fpa * rad2deg))
print('inclination only changes \n')
print(' vinit   %11.7f  %11.7f \n' % (vinit, vinit * velkmps))
print(' deltavionly   %11.7f  %11.7f \n' % (deltavionly, deltavionly * velkmps))
print(' comp  %11.7f  %11.7f \n' % (-deltavionly * velkmps * np.cos(deltai),
                                    deltavionly * velkmps * np.sin(deltai)))
print('-----part 3 \n')
deltai = 15.0 * deg2rad
ecc = 0.3
p = 17858.7836 / re
nu = 150.0 * deg2rad
mu = 1
a = p / (1.0 - ecc * ecc)
r = p / (1.0 + ecc * np.cos(nu))
vinit = np.sqrt((2.0 * mu) / r - (mu / a))
fpa = np.arctan((ecc * np.sin(nu)) / (1.0 + ecc * np.cos(nu)))
deltavionly = obu.ionlychg(deltai, vinit, fpa)
print(' a  %11.7f  %11.7f km r %11.7f  %11.7f km fpa %11.7f \n\n'
      % (a, a * re, r, r * re, fpa * rad2deg))
print('inclination only changes \n')
print(' deltavionly   %11.7f  %11.7f \n' % (deltavionly, deltavionly * velkmps))
print(' vinit   %11.7f  %11.7f \n' % (vinit, vinit * velkmps))
print(' comp  %11.7f  %11.7f \n' % (-deltavionly * velkmps * np.cos(deltai),
                                    deltavionly * velkmps * np.sin(deltai)))
# End of book example


print('-----other \n')
deltai = (90.0 - 28.5) * deg2rad
vinit = np.sqrt(398600.418 / 6600) / velkmps

incl = 28.5 * deg2rad
fpa = 0.0 * deg2rad
deltavionly = obu.ionlychg(deltai, vinit, fpa)
print('inclination only changes \n')
print(' vinit   %11.7f  %11.7f \n' % (vinit, vinit * velkmps))
print(' deltavionly  %11.7f  %11.7f \n' % (deltavionly, deltavionly * velkmps))
print(' comp  %11.7f  %11.7f \n' % (-deltavionly * velkmps * np.cos(deltai),
                                    deltavionly * velkmps * np.sin(deltai)))
