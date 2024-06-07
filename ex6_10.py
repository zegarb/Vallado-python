#     -----------------------------------------------------------------
#
#                              ex6_10.m
#
#  this file demonstrates example 6-10.
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

print('-------------------- problem ex 6-10 \n')
aint = 7143.51 / re
atgt = 42159.4855 / re
iint = 28.5 * deg2rad
itgt = 0.0 * deg2rad
deltai = itgt - iint
nodeint = 45.0 * deg2rad
arglatint = 15.0 * deg2rad
phasenew = np.pi - arglatint
truelon = 200.0 * deg2rad
ktgt = 0
kint = 1
print('combined maneuver \n' % ())
ttrans, tphase, dvphase, dvtrans1, dvtrans2, aphase = \
    obu.noncoplr(phasenew, aint, atgt, ktgt, kint, arglatint, nodeint, truelon,
                 deltai)
print(' ttrans  %11.7f \n' % (ttrans))
print(' tphase  %11.7f \n' % (tphase))
print(' dvphase  %11.7f  %11.7f \n' % (dvphase, dvphase * velkmps))
print(' dvtrans1  %11.7f  %11.7f \n' % (dvtrans1, dvtrans1 * velkmps))
print(' dvtrans2  %11.7f  %11.7f \n' % (dvtrans2, dvtrans2 * velkmps))
print(' aphase  %11.7f \n' % (aphase))
