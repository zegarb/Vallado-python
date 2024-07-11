import math
import numpy as np

from space_constants import *
import spacetime_utils as stu


#     -----------------------------------------------------------------
#
#                              Ex3_10
#
#  this file demonstrates example 3-10.
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
hr = 13
min = 22
sec = 45.98
print('hr %4i ' % hr)
print('min %4i ' % min)
print('sec %8.6f \n' % sec)

utsec = stu.hms2sec(hr, min, sec)

print('hms %11.7f \n' % utsec)

temp = utsec / 3600.0
hr  = np.fix(temp )
min = np.fix((temp - hr) * 60.0)
sec = (temp - hr - min/60.0) * 3600.0

print(' hr min sec %4i  %4i  %8.6f \n' % (hr, min, sec))


