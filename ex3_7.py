import math
import numpy as np

from space_constants import *
import spacetime_utils as stu


#     -----------------------------------------------------------------
#
#                              Ex3_7.m
#
#  this file demonstrates example 3-7.
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


year = 2004
mon  =   5
day  =  14
hr   =  10
min  =  43
sec  =   0.0
dut1 = -0.463326
dat  = 32
xp   =  0.0
yp   =  0.0
lod  =  0.0
timezone= 6

# -------- convtime    - convert time from utc to all the others
#, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, \
    jdtdb, jdtdbfrac = stu.convtime(year, mon, day, hr, min,
                                    sec, timezone, dut1, dat)

print('ut1 %8.6f tut1 %16.12f jdut1 %18.11f\n' % (ut1, tut1,
                                                  jdut1 + jdut1frac))
print('utc %8.6f\n' % (utc))
print('tai %8.6f\n' % (tai))
print('tt  %8.6f ttt  %16.12f jdtt  %18.11f\n' % (tt, ttt, jdtt + jdttfrac))
print('tdb %8.6f ttdb %16.12f jdtdb %18.11f\n' % (tdb, ttdb,
                                                  jdtdb + jdtdbfrac))

