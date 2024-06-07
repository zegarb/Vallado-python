#     -----------------------------------------------------------------
#
#                              Ex5_2.m
#
#  this file demonstrates example 5-2 and 5-4 for sun and moon rise/set.
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
import spacetime_utils as stu
import orbit_utils as obu
from space_constants import au, rad2deg, deg2rad


# --------  sun         - sun rise set
jd, jdfrac = stu.jday(1996, 3, 23, 0, 0, 0.0)
latgd = 40.0 * deg2rad
lon = 0.0 * deg2rad
whichkind = 's'
utsunrise, utsunset, error = obu.sunriset(jd + jdfrac, latgd, lon, whichkind)
print('sun sunrise %14.4f  %14.4f  sunset %14.4f %14.4f \n' %
      (utsunrise, (utsunrise - int(np.floor(utsunrise))) * 60, utsunset,
       (utsunset - int(np.floor(utsunset))) * 60))
jd, jdfrac = stu.jday(2011, 6, 25, 0, 0, 0.0)
latgd = 40.9 * deg2rad
lon = - 74.3 * deg2rad
whichkind = 's'
utsunrise, utsunset, error = obu.sunriset(jd + jdfrac, latgd, lon, whichkind)
print('sun sunrise %14.4f  %14.4f  sunset %14.4f %14.4f \n' %
      (utsunrise, (utsunrise - int(np.floor(utsunrise))) * 60, utsunset,
       (utsunset - int(np.floor(utsunset))) * 60))
