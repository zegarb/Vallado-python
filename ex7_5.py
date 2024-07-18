#     -----------------------------------------------------------------
#
#                              Ex7_5
#
#  this file implements example 7-5.note that the intialization has changed
#  from what is shown in the book, and the example there. The empiracle
#  formula seems to give better results for the intiial evaluation, so it's
#  used here. the book values of psi = 0.0 is not used here.
#
#                          companion code for
#             fundamentals of astrodynamics and applications
#                                 2020
#                            by david vallado
#
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#            26 may 20  david vallado
#                         separate from temp codes.
#  changes :
#            13 feb 07  david vallado
#                         original baseline
#     *****************************************************************

import math
import numpy as np
import orbit_utils as obu
from space_constants import au, rad2deg, re, velkmps, mu, twopi
import space_conversions as sc
import spacemath_utils as smu
import spacetime_utils as stu

fid = 1
#directory = 'd:\codes\library\matlab\'
#outfile = open(strcat(directory, 'tlambfig.out'), 'wt')

# ---------------------------- book tests -----------------------------
r1 = np.array([2.5, 0.0, 0.0]) * re
r2 = np.array([1.9151111, 1.606969, 0.0]) * re
# original orbit, assume circular
v1 = np.array([0, 0, 0])
print('\n-------- lambert test book pg 497 short way \n')
v1[1] = math.sqrt(mu / r1[0])
ang = np.arctan(r2[1] / r2[0])
v2 = np.array([[-np.sqrt(mu / r2[1]) * np.cos(ang)], [np.sqrt(mu / r2[0]) * np.sin(ang)], [0.0]])
print('\n v1 \n', (v1))
print('\n v2 \n', (v2))
dtsec = 76.0 * 60.0
# now show all the cases
magr1 = smu.mag(r1)
magr2 = smu.mag(r2)
# this value stays constant in all calcs, vara changes with df
cosdeltanu = np.dot(r1, r2) / (magr1 * magr2)
dm = 'S'
de = 'L'
nrev = 0
tbidu, tbiru = obu.lambgettbiu(r1, r2, 5)
print(' r1 \n',  (r1))
print(' r2 \n',  (r2))
print('From universal variables \n%11.7f %11.7f s \n' % (tbidu[0, 0], tbidu[0, 1]))
print('%11.7f %11.7f s \n' % (tbidu[1, 0], tbidu[1, 1]))
print('%11.7f %11.7f s \n' % (tbidu[2, 0], tbidu[2, 1]))
print('%11.7f %11.7f s \n' % (tbidu[3, 0], tbidu[3, 1]))
print('%11.7f %11.7f s \n\n' % (tbidu[4, 0], tbidu[4, 1]))
print('%11.7f %11.7f s \n' % (tbiru[0, 0], tbiru[0, 1]))
print('%11.7f %11.7f s \n' % (tbiru[1, 0], tbiru[1, 1]))
print('%11.7f %11.7f s \n' % (tbiru[2, 0], tbiru[2, 1]))
print('%11.7f %11.7f s \n' % (tbiru[3, 0], tbiru[3, 1]))
print('%11.7f %11.7f s \n' % (tbiru[4, 0], tbiru[4, 1]))
minenergyv, aminenergy, tminenergy, tminabs = obu.lambertmin(r1, r2, 'L', 0)
print(' minenergyv %16.8f %16.8f %16.8f a %11.7f  dt %11.7f  %11.7f \n' % (minenergyv[0], minenergyv[1], minenergyv[2], aminenergy, tminenergy, tminabs))
minenergyv, aminenergy, tminenergy, tminabs = obu.lambertmin(r1, r2, 'S', 0)
print(' minenergyv %16.8f %16.8f %16.8f a %11.7f  dt %11.7f  %11.7f \n' % (minenergyv[0], minenergyv[1], minenergyv[2], aminenergy, tminenergy, tminabs))
dtwait = 0.0

print('\n-------- lambertu test \n' % ())
v1t, v2t, errorl = obu.lambertu(r1, v1, r2, dm, de, nrev, dtwait, dtsec, tbidu, fid)
print(' v1t \n',  (v1t))
print(' v2t \n',  (v2t))
# run the 6 cases
print(' ------------- new time to accomodate 1 rev \n' % ())
dtsec = 21000.0
print(' TEST ------------------ S  L  0 rev \n' % ())
v1t, v2t, errorl = obu.lambertu(r1, v1, r2, 'S', 'L', 0, dtwait, dtsec, tbidu, fid)
print(' v1t \n',  (v1t))
print(' v2t \n',  (v2t))
print(' TEST ------------------ L  H  0 rev \n' % ())
v1t, v2t, errorl = obu.lambertu(r1, v1, r2, 'L', 'H', 0, dtwait, dtsec, tbiru, fid)
print(' v1t \n',  (v1t))
print(' v2t \n',  (v2t))
print(' TEST ------------------ S  L  1 rev \n' % ())

v1t, v2t, errorl = obu.lambertu(r1, v1, r2, 'S', 'L', 1, dtwait, dtsec, tbidu, fid)
print('uv1t \n',  (v1t))
print(' v2t \n',  (v2t))
print(' TEST ------------------ S  H  1 rev \n' % ())
v1t, v2t, errorl = obu.lambertu(r1, v1, r2, 'S', 'H', 1, dtwait, dtsec, tbidu, fid)
print('uv1t \n',  (v1t))
print(' v2t \n',  (v2t))
print(' TEST ------------------ L  L  1 rev \n' % ())
v1t, v2t, errorl = obu.lambertu(r1, v1, r2, 'L', 'L', 1, dtwait, dtsec, tbiru, fid)
print('uv1t \n',  (v1t))
print(' v2t \n',  (v2t))
print(' TEST ------------------ L  H  1 rev \n' % ())
v1t, v2t, errorl = obu.lambertu(r1, v1, r2, 'L', 'H', 1, dtwait, dtsec, tbiru, fid)
print('uv1t \n',  (v1t))
print(' v2t \n',  (v2t))


print('\n-------- lambertb test \n' % ())
v1t, v2t, errorl = obu.lambertb(r1, v1, r2, 's', 'l', 0, dtsec)
print(' v1dv \n',  (v1t))
print(' v2dv \n',  (v2t))

