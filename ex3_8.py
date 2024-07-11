import math
import numpy as np

from space_constants import *
import space_conversions as sc


#     -----------------------------------------------------------------
#
#                              Ex3_8
#
#  this file demonstrates example 3-8.
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


# -------- dms test
deg = -35
min = -15
sec = -53.63
print('deg %4i ' % deg)
print('min %4i ' % min)
print('sec %8.6f \n' % sec)

dms = sc.dms2rad(deg, min, sec)

print('dms %11.7f \n' % dms)


deg, min, sec = sc.rad2dms(dms)

print(' deg min sec %4i  %4i  %8.6f \n' % (deg, min, sec))


