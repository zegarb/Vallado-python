import math
import numpy as np

from space_constants import *
import space_conversions as sc

#     -----------------------------------------------------------------
#
#                              Ex3_9.m
#
#  this file demonstrates example 3-9.
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


# -------- hms test
hr = 15
min = 15
sec = 53.63
print('hr %4i ' % hr)
print('min %4i ' % min)
print('sec %8.6f \n' % sec)

hms = sc.hms2rad(hr, min, sec)

print('hms %11.7f %11.7f \n' % (hms, hms * 180.0/math.pi))

hr, min, sec = sc.rad2hms(hms)

print(' hr min sec %4i  %4i  %8.6f \n' % (hr, min, sec))



