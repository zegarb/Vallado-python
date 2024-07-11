#     -----------------------------------------------------------------
#
#                              ex6_1
#
#  this file demonstrates example 6-1.
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


print('-------------------- problem ex 6-1 \n' % ())
rinit = (re + 191.3411) / re
rfinal = (re + 35781.34857) / re
einit = 0.0
efinal = 0.0
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
print(' rinit  %11.6f %11.6f km \n' % (rinit, rinit * re))
print(' rfinal %11.6f %11.6f km \n' % (rfinal, rfinal * re))
print(' einit   %11.6f \n' % (einit))
print(' efinal  %11.6f \n' % (efinal))
print(' nuinit  %11.6f deg \n' % (nuinit * rad2deg))
print(' nufinal %11.6f deg \n' % (nufinal * rad2deg))
deltava, deltavb, dttu = obu.hohmann(rinit, rfinal, einit, efinal, nuinit,
                                   nufinal)
print('hohmann answers \n')
print(' deltava  %11.6f  %11.6f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.6f  %11.6f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltav   %11.6f  %11.6f km/s \n' % (deltava + deltavb, (deltava + deltavb)*velkmps))
print(' dttu  %11.6f tu %11.6f hr %11.6f min \n' % (dttu, dttu * tumin / 60.0,
                                                    dttu * tumin))
