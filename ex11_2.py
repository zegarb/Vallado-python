#     -----------------------------------------------------------------
#
#                              Ex11_2.py
#
#  this file demonstrates example 11-2.
#
#                          companion code for
#             fundamentals of astrodynamics and applications
#                                 2013
#                            by david vallado
#
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#             4 nov 15  david vallado
#                         original
#  changes :
#             4 nov 15  david vallado
#                         original baseline
#
#     *****************************************************************

import numpy as np
from space_constants import *


print("---------------------Example 11-2--------------------")

j2 = 0.00108263
rate = (360.0 / 365.2421897) * deg2rad / 86164.0989036973
print(f'{rate=}')
# intermediate conversions
print(360.0 / 365.2421897)
print((360.0 / 365.2421897) * deg2rad)
rate * tusec
reca = 1.0 / 3.5
# ---- ex 11-2 ---- }
# 800 km altitude
a = 7178.1363
ecc = 0.0
incl = np.arccos(- a ** 3.5 * 2 * rate * (1.0 - ecc ** 2) ** 2 / (3.0 * re ** 2 * j2 * np.sqrt(mu)))
print('i %11.5f a %11.5f  e %11.7f %11.7f  %11.7f \n' % (incl * rad2deg, a, ecc, a * (1.0 + ecc) - re, a * (1.0 - ecc) - re))
# change ecc
incl = 98.627 * deg2rad
ecc = 0.02
temp = - 1.5 * j2 * np.cos(incl) * re ** 2 * np.sqrt(mu) / (rate * (1.0 - ecc * ecc) ** 2)
a = temp ** reca
print('i %11.5f a %11.5f  e %11.7f %11.7f  %11.7f \n' % (incl * rad2deg, a, ecc, a * (1.0 + ecc) - re, a * (1.0 - ecc) - re))
# now verify alternate equations
a = 7179.82095
ecc = 0.02
incl = np.arccos(-a ** 3.5 * 2 * rate * (1.0 - ecc ** 2) ** 2 / (3.0 * re ** 2 * j2 * np.sqrt(mu)))
incl * rad2deg
a = 7179.82095
incl = 98.627 * deg2rad
ecc = np.sqrt(1.0 - np.sqrt(- 3.0 * re ** 2 * j2 * np.sqrt(mu) * np.cos(incl) / (2.0 * a ** 3.5 * rate)))
print(f'{ecc=}')
