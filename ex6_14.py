#     -----------------------------------------------------------------
#
#                              Ex6_14
#
#  this file demonstrates example 6-14.
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
from space_constants import au, re, velkmps, tumin, mu
import space_conversions as sc
import spacemath_utils as smu

print('-------------------- problem ex 6-14 \n')
ro = np.array([0.0, 0.0, 0.0])
vo = np.array([-0.1, -0.04, -0.02])

ralt = 590.0

dtsec = 0.0
rint, vint = obu.hillsr(ro, vo, ralt, dtsec)
print('initial interceptor position \n')
print(' r  %11.7f  %11.7f  %11.7f m \n' % (rint[0], rint[1], rint[2]))
print(' v  %11.7f  %11.7f  %11.7f m/s \n\n' % (vint[0], vint[1], vint[2]))
dtsec = 300.0
rint, vint = obu.hillsr(ro, vo, ralt, dtsec)
print(' r  %11.7f  %11.7f  %11.7f \n' % (rint[0], rint[1], rint[2]))
print(' v  %11.7f  %11.7f  %11.7f \n\n' % (vint[0], vint[1], vint[2]))
rint, vint = obu.hillsr(ro * 0.001, vo * 0.001, ralt, dtsec)
print(' r  %11.7f  %11.7f  %11.7f km \n' % (rint[0], rint[1], rint[2]))
print(' v  %11.7f  %11.7f  %11.7f km/s \n\n' % (vint[0], vint[1], vint[2]))
dtsec = 1200.0
rint, vint = obu.hillsr(ro, vo, ralt, dtsec)
print(' r  %11.7f  %11.7f  %11.7f \n' % (rint[0], rint[1], rint[2]))
print(' v  %11.7f  %11.7f  %11.7f \n\n' % (vint[0], vint[1], vint[2]))
dtsec = 600.0
rint, vint = obu.hillsr(ro, vo, ralt, dtsec)
print(' r  %11.7f  %11.7f  %11.7f \n' % (rint[0], rint[1], rint[2]))
print(' v  %11.7f  %11.7f  %11.7f \n\n' % (vint[0], vint[1], vint[2]))
