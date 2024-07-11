#     -----------------------------------------------------------------
#
#                              Ex3_1
#
#  this file demonstrates example 3-1.
#
#                          companion code for
#             fundamentals of astrodyanmics and applications
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
from space_constants import *
import numpy as np

# -------------------------  implementation   -----------------
latgd = 39.586667 * deg2rad
lon = - 105.64 * deg2rad
alt = 4.3476672

sinlat = np.sin(latgd)
# ------  find rdel and rk components of site vector  ---------
cearth = re / np.sqrt(1.0 - (eccearthsqrd * sinlat * sinlat))
rdel = (cearth + alt) * np.cos(latgd)
rk = ((1.0 - eccearthsqrd) * cearth + alt) * sinlat
print('c %16.8f s %11.8f km \n' % (cearth, cearth * (1.0 - eccearthsqrd)))
print('rdelta %16.8f rk %11.8f km \n' % (rdel, rk))
