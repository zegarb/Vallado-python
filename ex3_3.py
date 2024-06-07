import math
import numpy as np

from space_constants import *
from space_conversions import ecef2ll, ecef2llb, dms2rad
import spacemath_utils as smu
import orbit_utils as obu


#     -----------------------------------------------------------------
#
#                              Ex3_3.m
#
#  this file demonstrates example 3-3.
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


# --------  ecef2ll       - position to lat lon alt almanac (fastest)
r1=np.array([6524.834, 6862.875, 6448.296])

latgc, latgd, lon, hellp = ecef2ll(r1)

print('ecef2ll  gc %14.7f gd %14.7f %14.7f%14.7f\n'
      % (latgc * rad2deg, latgd * rad2deg, lon * rad2deg, hellp))

# --------  ecef2llb      - position to lat lon alt borkowski
latgc, latgd, lon, hellp = ecef2llb(r1)

print('ecef2llb gc %14.7f gd %14.7f %14.7f%14.7f\n'
      % (latgc * rad2deg, latgd * rad2deg, lon * rad2deg, hellp))

# try another case for testing the borkowski method
# Latitude
deg =  -7
min =  -54
sec = -23.886
latgd = dms2rad(deg, min, sec)

# Longitude
deg = 345
min =  35
sec =  51.000
lon = dms2rad(deg, min, sec)

# Altitude
hellp = 0.056
rs, vs = obu.site(latgd, lon, hellp)
latgc, latgd, lon, hellp = ecef2ll(rs)

print('ecef2ll  gc %14.7f gd %14.7f %14.7f%14.7f\n'
      % (latgc * rad2deg, latgd * rad2deg, lon * rad2deg, hellp))

latgcb, latgdb, lonb, hellpb = ecef2llb(rs)

print('ecef2llb gc %14.7f gd %14.7f %14.7f%14.7f\n'
      % (latgcb * rad2deg, latgdb * rad2deg, lonb * rad2deg, hellpb))

# check for any imaginary numbers
print(latgcb)
print(latgdb)
print(lonb)
print(hellpb)
