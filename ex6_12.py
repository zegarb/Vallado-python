#     -----------------------------------------------------------------
#
#                              ex6_12
#
#  this file demonstrates example 6-12.
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

print('-------------------- problem ex 6-12 \n')
rinit = (re + 318.9)
rfinal = (re + 35781.343)
einit = 0.0
efinal = 0.0
iinit = 0.0 * deg2rad
ifinal = 0.0 * deg2rad

#eq 6-50
dvacc = 1 - np.sqrt(1 / (rfinal / rinit))

#eq 6-45
athrust = 4.501e-4
pfm = 0.25
mspec = (athrust / dvacc) * np.log(1 - pfm)
tof = pfm / np.abs(mspec)

#convert from TU to days
tof = tof * 868.0722 * 1.15741e-5

print(' rinit  %11.7f m %11.7f er \n' % (rinit, rinit / re))
print(' rfinal %11.7f m %11.7f er \n' % (rfinal, rfinal / re))
print(' einit   %11.7f \n' % (einit))
print(' efinal  %11.7f \n' % (efinal))
print(' iinit   %11.7f deg \n' % (iinit * rad2deg))
print(' ifinal  %11.7f deg \n' % (ifinal * rad2deg))
#deltav,tof = lowthrust(rinit,rfinal,iinit,ifinal,5000.0,- 2.48e-07,4e-06,0.001)
#print('low thrust answers \n' % ())
#print(' deltav  %11.7f km/s  \n' % (deltav))
print(' dttu  %11.7f day %11.7f min \n' % (tof, tof * 1440))
