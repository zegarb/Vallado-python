#     -----------------------------------------------------------------
#
#                              Ex11_5.m
#
#  this file demonstrates example 11-5.
#
#                          companion code for
#             fundamentals of astrodynamics and applications
#                                 2007
#                            by david vallado
#
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#            30 aug 11  david vallado
#                         original
#  changes :
#            22 aug 11 david vallado
#                         original baseline
#
#     *****************************************************************

import numpy as np
from space_constants import *

j2 = 0.00108263
# --------  repeat gt calculations
a = 6570.3358

ecc = 0.006301
incl = 45.0 * deg2rad
p = a * (1.0 - ecc * ecc)
nanom = np.sqrt(mu / (a * a * a))
n = nanom
print('n %11.7f  %11.7f \n' % (n, n * deg2rad))
argprate = 0.75 * j2 * re ** 2 * n * (4.0 - 5.0 * (np.sin(incl)) ** 2) / (p * p)
print('argprate %11.7e  %11.7f deg/day \n' % (argprate, argprate * rad2deg * 86400))

dargpmax = np.sqrt(2 * 30 / (a * ecc) * ((1.0 + ecc) / (1.0 - ecc)))
taumax = np.abs(dargpmax / argprate)
dv = n * a * ecc / 2 * dargpmax
print('dargpmax %11.7f %11.7f taumax %11.7f %11.7f  dv %11.7f \n' % (dargpmax, dargpmax * rad2deg, taumax / 60, taumax / 86400, dv))

print('frozen eccentricity design')
dargpmax = np.sqrt(2 * 10 / (a * ecc) * ((1.0 + ecc) / (1.0 - ecc)))
taumax = np.abs(dargpmax / argprate)
dv = n * a * ecc / 2 * dargpmax
print('dargpmax %11.7f %11.7f taumax %11.7f %11.7f  dv %11.7f \n' % (dargpmax, dargpmax * rad2deg, taumax / 60, taumax / 86400, dv))
