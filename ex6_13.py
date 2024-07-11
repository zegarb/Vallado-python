#     -----------------------------------------------------------------
#
#                              ex6_13
#
#  this file demonstrates example 6-13.
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
from space_constants import *
import space_conversions as sc
import spacemath_utils as smu

print('-------------------- problem ex 6-13 \n')
rinit = (re + 200.0)

rfinal = (re + 35781.343)
einit = 0.0
efinal = 0.0
iinit = 28.5 * deg2rad
ifinal = 0.0 * deg2rad
print(' rinit  %11.7f m %11.7f er \n' % (rinit, rinit / re))
print(' rfinal %11.7f m %11.7f er \n' % (rfinal, rfinal / re))
print(' einit   %11.7f \n' % (einit))
print(' efinal  %11.7f \n' % (efinal))
print(' iinit   %11.7f deg \n' % (iinit * rad2deg))
print(' ifinal  %11.7f deg \n' % (ifinal * rad2deg))
# book example
du = 6578.136
du1 = du / re
tu = np.sqrt(du1 ** 3 / 1.0) * tusec
ratio = rfinal / rinit
vacc = 0.75 * du / tu
tof = -1 / (-2e-07) * (1 - np.exp((-2e-07 * 5.8382) / 1e-05))
print('tof: %f\n' % tof)
#[deltav, tof] = lowthrust(rinit, rfinal, iinit, ifinal, 3800.0, 0.456, 0.234876);
# deltav, tof = lowthrust(rinit, rfinal, iinit, ifinal, 5000.0, -2e-07,
                            # 1e-05, -0.54)
#low_thrust_transfer_prop(rinit, rfinal, iinit, ifinal, 3800.0, 0.456, 0.234876);

# print('low thrust answers \n')
# print(' deltav  %11.7f  km/s \n' % (deltav))
# print(' dttu  %11.7f day %11.7f min \n' % (tof,tof * 1440))
