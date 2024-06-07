#     -----------------------------------------------------------------
#
#                              ex6_2.m
#
#  this file demonstrates example 6-2.
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

print('-------------------- problem ex 6-2 \n' % ())
rinit = (re + 191.3411) / re
rb = (re + 503873.0) / re
rfinal = (re + 376310.0) / re
einit = 0.0
efinal = 0.0
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
print('initial position \n' % ())
print(' rinit  %11.7f %11.7f km \n' % (rinit, rinit * re))
print(' rfinal %11.7f %11.7f km \n' % (rfinal, rfinal * re))
print(' einit   %11.7f \n' % (einit))
print(' efinal  %11.7f \n' % (efinal))
print(' nuinit  %11.7f deg \n' % (nuinit * rad2deg))
print(' nufinal %11.7f deg \n' % (nufinal * rad2deg))
deltava, deltavb, deltavc, dttu = \
    obu.biellip(rinit, rb, rfinal, einit, efinal, nuinit, nufinal)
print('bi-elliptic answers \n')
print(' deltava  %11.7f  %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f  %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltavc  %11.7f  %11.7f km/s \n' % (deltavc, deltavc * velkmps))
deltav = deltava + deltavb + deltavc
print(' deltav   %11.7f  %11.7f km/s \n' % (deltav, deltav * velkmps))
print(' dttu  %11.6f tu %11.6f hr %11.6f min \n'
      % (dttu, dttu * tumin / 60.0, dttu * tumin))
