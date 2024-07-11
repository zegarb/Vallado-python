#     -----------------------------------------------------------------
#
#                              Ex5_7
#
#  this file demonstrates example 5-7.
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
#            21 feb 19  david vallado
#                         update for new constants
#  changes :
#            13 feb 07  david vallado
#                         original baseline
#
#     *****************************************************************

import numpy as np
from space_constants import *
import spacemath_utils as smu

# Helio for Neptune from the Sun
paral1 = np.arcsin(1.0 / 29.664361)
print('paralax  %11.9f %11.9f arcsec \n' % (paral1 * rad2deg, paral1 * rad2arcsec))
# try exact method (neptune to earth)

rnep = np.array([-1757460712.2509, 3757470099.1416, 1576777174.1537])

#  neptune to sun vector from JPL epehmerides 1994, 5/14
rearth = np.array([-1666604612.0985, 3868340828.5807, 1624846829.1305])

paral2 = np.arccos(np.dot(rearth, rnep) / (smu.mag(rearth) * smu.mag(rnep)))
print('paralax  %11.9f %11.9f arcsec \n' % (paral2 * rad2deg, paral2 * rad2arcsec))
paral1 / paral2 * 100
# for Neptune to Earth
paral1 = np.arcsin(6378.0 / (149597870 * 29.664361))
print('paralax  %11.9f %11.9f arcsec \n' % (paral1 * rad2deg, paral1 * rad2arcsec))
# try exact method
# neptune to Earth
rnep = np.array([- 1757460712.2509, 3757470099.1416, 1576777174.1537])

# neptune to site
rsite = np.array([6378.137, 0.0, 0.0])
rearth = rnep + rsite
paral2 = np.arccos(np.dot(rnep, rearth) / (smu.mag(rnep) * smu.mag(rearth)))
print('paralax  %11.9f %11.9f arcsec \n' % (paral2 * rad2deg, paral2 * rad2arcsec))
paral1 / paral2 * 100
