import math
import numpy as np

from space_constants import *
from space_conversions import dms2rad, gd2gc
import spacemath_utils as smu


#     -----------------------------------------------------------------
#
#                              Ex3_2.m
#
#  this file demonstrates example 3-2.
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



deg =  -7
min =  -54
sec = -23.886
print('deg %2i min %2i sec %8.6f\n' % (deg, min, sec))
latgd = dms2rad(deg, min, sec)
print('dms = %11.7f rad\n' % (latgd * rad2deg))

deg = 345
min =  35
sec =  51.000
print('deg %2i min %2i sec %8.6f\n' % (deg, min, sec))
lon = dms2rad(deg, min, sec)
print('dms = %11.7f rad\n' % (lon * rad2deg))

alt = 0.056

# -------------------------  implementation   -----------------
sinlat      = math.sin(latgd)
#earthrate[0]= 0.0
#earthrate[1]= 0.0
#earthrate[2]= omegaearth

# ------  find rdel and rk components of site vector  ---------
cearth= re / math.sqrt(1.0 - (eccearthsqrd * sinlat * sinlat))
rdel  = (cearth + alt) * math.cos(latgd)
rk    = ((1.0 - eccearthsqrd) * cearth + alt) * sinlat

#(1.0-eccearthsqrd)*cearth

# ---------------  find site position vector  -----------------
rs = np.zeros((3))
rs[0] = rdel * math.cos(lon)
rs[1] = rdel * math.sin(lon)
rs[2] = rk
rs = rs.T

rsmag = smu.mag(rs)
print('site gd %16.9f %16.9f %16.9f   %16.9f \n'
      % (rs[0], rs[1], rs[2], rsmag))

latgc = gd2gc(latgd)

r = np.zeros((3))
r[0]= rsmag * math.cos(latgc) * math.cos(lon)
r[1]= rsmag * math.cos(latgc) * math.sin(lon)
r[2]= rsmag * math.sin(latgc)
r = r.T

print('site gc %16.9f %16.9f %16.9f  %16.9f  \n'
      % (r[0], r[1], r[2], latgc * rad2deg))

