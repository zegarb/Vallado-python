#     -----------------------------------------------------------------
#
#                              ex6_8.m
#
#  this file demonstrates example 6-8.
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

print('-------------------- problem ex 6-8 \n')
rcs1 = 12756.274 / re
rcs3 = 12756.274 / re
phasei = - 20.0 * deg2rad
einit = 0.0
efinal = 0.0
nuinit = 0.0 * deg2rad
nufinal = 0.0 * deg2rad
ktgt = 1
kint = 1
print('\n rendezvous ktgt %i kint %i \n' % (ktgt, kint))
phasef, waittime, deltav = obu.rendz(rcs1, rcs3, phasei, einit, efinal, nuinit,
                                     nufinal, ktgt, kint)
print(' phasef %11.7f  \n' % (phasef * rad2deg))
print(' waittime %11.7f  %11.7f min  \n' % (waittime, waittime * tumin))
print(' deltav %11.7f  %11.7f km/s \n' % (deltav, deltav * velkmps))
ktgt = 2
kint = 2
print('\n rendezvous ktgt %i kint %i \n' % (ktgt, kint))
phasef, waittime, deltav = obu.rendz(rcs1, rcs3, phasei, einit, efinal, nuinit,
                                     nufinal, ktgt, kint)
print(' phasef %11.7f  \n' % (phasef * rad2deg))
print(' waittime %11.7f  %11.7f min  \n' % (waittime, waittime * tumin))
print(' deltav %11.7f  %11.7f km/s \n' % (deltav, deltav * velkmps))
