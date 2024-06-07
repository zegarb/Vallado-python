import math
import numpy as np

from space_constants import *
import spacetime_utils as stu



#     -----------------------------------------------------------------
#
#                              Ex3_13.m
#
#  this file demonstrates example 3-13.
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



# -------- 1995 June 8, 20:18:3.703691
jd = 2449877.3458762

print('jd %8.7f  \n\n' % jd)

jdfrac = jd-np.floor(jd)
jd = np.floor(jd)

print('integer %6i  \n' % jd)
print('fractional  %3.7f  \n\n' % jdfrac)

year, mon, day, hr, min, sec = stu.invjday(jd, jdfrac)

print('year %6i  \n' % year)
print('mon  %3i  \n' % mon)
print('day  %3i  \n' % day)
print('hr   %3i  \n' % hr)
print('min  %3i  \n' % min)
print('sec  %3.6f  \n' % sec)

# some stressing cases
for i in range(11):
    if (i == 1):
        jd = 2457884.5
        jdF = 0.160911640046296
    if (i == 2):
        jd = 2457884.5
        jdF = -0.160911640046296
    if (i == 3):
        jd = 2457884.5
        jdF = 0.660911640046296
    if (i == 4):
        jd = 2457884.5
        jdF = -0.660911640046296
    if (i == 5):
        jd = 2457884.5
        jdF = 2.160911640046296
    if (i == 6):
        jd = 2457884.660911640046296
        jdF = 0.0
    if (i == 7):
        jd = 2457884.0
        jdF = 2.160911640046296
    if (i == 8):
        jd = 2457884.5
        jdF = 0.0
    if (i == 9):
        jd = 2457884.5
        jdF = 0.5
    if (i == 10):
        jd = 2457884.5
        jdF = 1.0
    if (i == 0):
        jd = 2457884.3
        jdF = 1.0

    jdb = np.floor(jd + jdF) + 0.5
    mfme = (jd + jdF - jdb) * 1440.0
    if (mfme < 0.0):
        mfme = 1440.0 + mfme
    year, mon, day, hr, min, sec = stu.invjday(jd, jdF)
    if (abs(jdF) >= 1.0):
        jd = jd + np.floor(jdF)
        jdF = jdF - np.floor(jdF)
    dt = jd - np.floor(jd) - 0.5
    if (abs(dt) > 0.00000001):
        jd = jd - dt
        jdF = jdF + dt
    # this gets it even to the day
    if (jdF < 0.0):
        jd = jd - 1.0
        jdF = 1.0 + jdF
    print('%2i %4i %2i %2i %2i:%2i:%7.4f   %9.4f %9.4f %12.4f %8.5f %5.4f \n'
           % (i, year, mon, day, hr, min, sec, mfme, hr*60.0+min+sec/60.0,
              jd, jdF, dt))

# through stressing cases


jdo, jdfraco = stu.jday(2017, 8, 23, 12, 15, 16)

jdo, jdfraco = stu.jday(2017, 12, 31, 12, 15, 16)
print('%11.6f  %11.6f \n' % (jdo, jdfraco))
year, mon, day, hr, min, sec = stu.invjday(jdo + jdfraco, 0.0)
print('%4i  %3i  %3i  %2i:%2i:%6.4f  \n' % (year, mon, day, hr, min, sec))
year, mon, day, hr, min, sec = stu.invjday(jdo, jdfraco)
print('%4i  %3i  %3i  %2i:%2i:%6.4f  \n' % (year, mon, day, hr, min, sec))

for i in range(-50, 51):
    jd = jdo + i/10.0
    jdfrac = jdfraco
    year, mon, day, hr, min, sec = stu.invjday ( jd + jdfrac, 0.0 )
    print('%4i  %3i  %3i  %2i:%2i:%6.4f  \n' %(year, mon, day, hr, min, sec))
    year, mon, day, hr, min, sec = stu.invjday ( jd, jdfrac )
    dt = jd - np.floor(jd)
    print('%4i  %3i  %3i  %2i:%2i:%6.4f  \n\n' %(year, mon, day, hr, min, sec))

print('end first half \n')

for i in range(-50, 51):
    jd = jdo
    jdfrac = jdfraco + i/10.0
    year, mon, day, hr, min, sec = stu.invjday ( jd + jdfrac, 0.0 )
    print('%4i  %3i  %3i  %2i:%2i:%6.4f  \n' % (year, mon, day, hr, min, sec))
    year, mon, day, hr, min, sec = stu.invjday ( jd, jdfrac )
    print('%4i  %3i  %3i  %2i:%2i:%6.4f  \n\n' % (year, mon, day, hr, min, sec))

