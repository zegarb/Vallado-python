#     -----------------------------------------------------------------
#
#                              Ex5_4.m
#
#  this file demonstrates example 5-4.
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
import spacetime_utils as stu
import space_conversions as sc
import orbit_utils as obu
from space_constants import au, rad2deg, deg2rad

# --------  moon - moon rise/set
jd, jdfrac = stu.jday(1998, 8, 21, 0, 0, 0.0)
latgd = 40.0 * deg2rad
lon = 0.0 * deg2rad
utmoonrise, utmoonset, moonphaseang, error = obu.moonriset(jd + jdfrac, latgd,
                                                          lon)
print('moon moonrise %14.6f  %14.4f   moonset %14.6f %14.4f  \n'
      % (utmoonrise, (utmoonrise - np.floor(utmoonrise)) * 60,
         utmoonset, (utmoonset - np.floor(utmoonset)) * 60))
print('moon phase angle %14.4f   \n' % (moonphaseang))

jd, jdfrac = stu.jday(1990, 3, 5, 0, 0, 0.0)
latgd = 40.94 * deg2rad
lon = - 73.97 * deg2rad
utmoonrise, utmoonset, moonphaseang, error = obu.moonriset(jd + jdfrac, latgd,
                                                          lon)
print('moon moonrise %14.4f  %14.4f   moonset %14.4f %14.4f  \n'
      % (utmoonrise, (utmoonrise - np.floor(utmoonrise)) * 60,
         utmoonset, (utmoonset - np.floor(utmoonset)) * 60))
print('moon phase angle %14.4f   \n' % (moonphaseang))

jd, jdfrac = stu.jday(2006, 6, 28, 0, 0, 0.0)
latgd = 40.0 * deg2rad
lon = 0.0 * deg2rad
utmoonrise, utmoonset, moonphaseang, error = obu.moonriset(jd + jdfrac,
                                                          latgd, lon)
print('moon moonrise %14.4f    moonset %14.4f hrs \n'
      % (utmoonrise, utmoonset))
print('moon phase angle %14.4f   \n' % (moonphaseang))

jd = 2458291.0
jdfrac = 0.0
latgd = - 108.2802963256836 * deg2rad
lon = 32.77009963989258 * deg2rad
utmoonrise, utmoonset, moonphaseang, error = obu.moonriset(jd + jdfrac,
                                                          latgd, lon)
#fprintf(1, 'moon moonrise #14.4f    moonset #14.4f hrs \n', utmoonrise, utmoonset );
print('moon moonrise %14.4f  %14.4f   moonset %14.4f %14.4f  \n' %
      (utmoonrise, (utmoonrise - np.floor(utmoonrise)) * 60,
       utmoonset, (utmoonset - np.floor(utmoonset)) * 60))
print('moon phase angle %14.4f   \n' % (moonphaseang))

print('     40    42    44    46    48    50    52    54    56    58    60    62    64    66  \n' % ())
for i in range(8, 31):
    jd, jdfrac = stu.jday(2006, 6, i, 0, 0, 0.0)
    for j in range(0, 13):
        latgd = (40.0 + j * 2.0) * deg2rad
        lon = 0.0 * deg2rad
        utmoonrise, utmoonset, moonphaseang, error = \
            obu.moonriset(jd + jdfrac, latgd, lon)
            # if strcmp(error, 'ok') == 0 # 1 if true, 0 if false
                # fprintf(1, 'error');
        hr, min, sec = sc.rad2hms(utmoonrise * 15.0 * np.pi / 180.0)
        hr1, min1, sec1 = sc.rad2hms(utmoonset * 15.0 * np.pi / 180.0)
        #  print out header date for each section of results
        if j == 0:
            print('%2i ' % (i))
                # if utmoonrise > 9998.0 utmoonset > 9998.0
            #  fprintf(1, ' none ');
        # else
        # if utmoonrise > 9998.0 && utmoonset < 24.0
            #  fprintf(1, ' nors ');
        # else
        # if utmoonset - utmoonrise > 14.0
            # fprintf(1, ' none ');
        # else
            # fprintf(1, '#2i:#2i ', hr, min);
        # if hr >= 24
        #jd, jdfrac = stu.jday(2006, 6, i, 0, 0, 0.0)
        #el1 = obu.moonel(jd + jdfrac, latgd, lon)
        #jd, jdfrac = stu.jday(2006, 6, i + 1, 0, 0, 0.0)
        #el2 = obu.moonel(jd + jdfrac, latgd, lon)
        if hr >= 24:
            print('| nost  ')
        else:
            print('| %2i:%2i ' % (hr, min))
    print('\n')

print('     40    42    44    46    48    50    52    54    56    58    60    62    64    66  \n')
for i in range(8, 31):
    jd, jdfrac = stu.jday(2006, 6, i, 0, 0, 0.0)
    for j in range(0, 13):
        latgd = (40.0 + j * 2.0) * deg2rad
        lon = 0.0 * deg2rad
        utmoonrise, utmoonset, moonphaseang, error = \
            obu.moonriset(jd + jdfrac, latgd, lon)
            # if strcmp(error, 'ok') == 0 # 1 if true, 0 if false
                # fprintf(1, 'error');
            # end;
        hr, min, sec = sc.rad2hms(utmoonrise * 15.0 * np.pi / 180.0)
        hr1, min1, sec1 = sc.rad2hms(utmoonset * 15.0 * np.pi / 180.0)
        #  print out header date for each section of results
        if j == 0:
            print('%2i ' % (i))
        #jd, jdfrac = stu.jday(2006, 6, i, 0, 0, 0.0)
        #el1 = obu.moonel(jd + jdfrac, latgd, lon)
        #jd, jdfrac = stu.jday(2006, 6, i + 1, 0, 0, 0.0)
        #el2 = obu.moonel(jd + jdfrac, latgd, lon)
        if hr1 >= 24:
            print('| nost  ')
        else:
            print('| %2i:%2i ' % (hr1, min1))
    print('\n')


