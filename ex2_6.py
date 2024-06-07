from space_conversions import coe2rv
from space_constants import *


#     -----------------------------------------------------------------
#
#                              Ex2_6.m
#
#  this file demonstrates example 2-6.
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



print('coe test ----------------------------\n' )
p     = 11067.790 # km
ecc   = 0.83285
incl  = 87.87 * deg2rad
omega = 227.89 * deg2rad
argp  = 53.38 * deg2rad
nu    = 92.335 * deg2rad
arglat = 0.0
truelon= 0.0
lonper = 0.0

a = p / (1 - ecc * ecc)

print('          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper')
print('coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f\n' % \
    (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg,
    arglat * rad2deg, truelon * rad2deg, lonper * rad2deg ))

# --------  coe2rv       - classical elements to posisiotn and velocity

### coe2rvS.m is incomplete and only just replaces || with |. -jmb
#r, v = coe2rvS(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper)
#print('r    %15.9f %15.9f %15.9f' % (r[0], r[1], r[2]))
#print(' v %15.10f %15.10f %15.10f\n' % (v[0], v[1], v[2]))

r, v = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper)
print('r    %15.9f %15.9f %15.9f' % (r[0], r[1], r[2]))
print('v %15.10f %15.10f %15.10f\n' % (v[0], v[1], v[2]))
