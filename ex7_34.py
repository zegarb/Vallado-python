#     -----------------------------------------------------------------
#
#                              Ex7_3.py
#
#  this file demonstrates example 7-3 and 7-4. it also compares the gibbs and
#  herrick gibbs appraches.
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
#             7 jun 07  david vallado
#                         original
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
import spacetime_utils as stu


print('-------------------- problem ex 7-3 \n')
r1 = np.array([0.0, 0.0, 6378.137])
r2 = np.array([0.0, - 4464.696, - 5102.509])
r3 = np.array([0.0, 5740.323, 3189.068])
print(' r1 ', (r1))
print(' r2 ', (r2))
print(' r3 ', (r3))
v2g, theta, theta1, copa, errorg = obu.gibbs(r1, r2, r3)
print(' v2g km/s \n', (v2g))
print(' theta %11.6f  theta1  %11.6f  copa %11.6f \n'
      % (theta * rad2deg, theta1 * rad2deg, copa * rad2deg))
print('-------------------- problem ex 7-4 \n')
r1 = np.array([3419.85564, 6019.82602, 2784.60022])
r2 = np.array([2935.91195, 6326.18324, 2660.59584])
r3 = np.array([2434.95202, 6597.38674, 2521.52311])
print(' r1 ', (r1))
print(' r2 ', (r2))
print(' r3 ', (r3))
# just needs to be in secs - so precision isn't lost
jd1 = 0.0 / 86400.0
jd2 = (60.0 + 16.48) / 86400.0
jd3 = (120.0 + 33.04) / 86400.0
# use t1, t2, t3 as secs for greater accuracy
v2h, theta, theta1, copa, errorh = obu.hgibbs(r1, r2, r3, jd1, jd2, jd3)
print(' v2h km/s \n\n\n' , (v2h))
print(' theta %11.6f  theta1  %11.6f  copa %11.6f \n'
      % (theta * rad2deg, theta1 * rad2deg, copa * rad2deg))

### END BOOK

#  this section tests the gibbs routine

print('-------------------- compare gibbs and hgibbs \n')
r = np.zeros(3)
r[0] = 1.0
r[1] = 0.1
r[2] = 0.976
r = r * re
magr = smu.mag(r)
v = np.zeros(3)
v[0] = 0.3
v[1] = 0.1
v[2] = 0.276
v = v * velkmps
magv = smu.mag(v)
p = 1.61 * re
ecc = 0.0
incl = 34.0 * deg2rad
omega = 0.0
argp = 0.0
nu = 0.0
arglat = 0.0
truelon = 0.0
lonper = 0.0
r, v = sc.coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper)
# use kepler propagation
#        r = [-4550.4256  2220.1946  5090.3547];
#        v = [-4.975006  -5.237109  -2.11811];
#        for i = 0:50
#            [r1, v1, error] =  kepler  ( r, v, i*600.0 );
#            fprintf(1, '#4i x #11.7f  #11.7f  #11.7f  #11.7f  #11.7f  #11.7f \n', i, r1, v1 );
#          end;

# test gibbs hgibbs accuracies against two body orbits
i = 1
t2 = 0.0

magr = smu.mag(r)
magv = smu.mag(v)
rdotv = np.dot(r, v)
# -------------  find sme, alpha, and a  ------------------
sme = ((magv ** 2) * 0.5) - (mu / magr)
a = - mu / (2.0 * sme)
period = twopi * np.sqrt(np.abs(a) ** 3.0 / mu)

print(' ktr    dt sec   fraction       ang1          ang2   '
      '    coplanar        gibbs      err flg             hgibbs      '
      '     errflg')
while t2 < period * 0.45:

    if t2 < period * 0.001:
        t2 = t2 + period * 5e-06
    elif t2 < period * 0.006:
        t2 = t2 + period * 0.001
    elif t2 < period * 0.015:
        t2 = t2 + period * 0.005
    elif t2 < period * 0.2:
        t2 = t2 + period * 0.01
    else:
        t2 = t2 + period * 0.05
    t1 = 0.0
    t2 = t2
    t3 = 2.0 * t2
    r1, v1, errk = obu.kepler(r, v, t1)
    # print(errk)
    r2, v2, errk = obu.kepler(r, v, t2)
    # print(errk)
    r3, v3, errk = obu.kepler(r, v, t3)
    # print(errk)
    v2g, theta, theta1, copa, errorg = obu.gibbs(r1, r2, r3)
    # just needs to be in days - so precision isn't lost
    jd1 = t1 / 86400.0
    jd2 = t2 / 86400.0
    jd3 = t3 / 86400.0
    # use t1, t2, t3 as secs for greater accuracy
    v2h, theta, theta1, copa, errorh = obu.hgibbs(r1, r2, r3, jd1, jd2, jd3)
    #v2 = v2 / velkmps;
#v2g = v2g / velkmps;
#v2h = v2h / velkmps;
    #            fprintf(1, 'v2truth #14.7f  #14.7f  #14.7f  \n', v2(1), v2(2), v2(3) );
#            fprintf(1, 'v2g     #14.7f  #14.7f  #14.7f  \n', v2g(1), v2g(2), v2g(3) );
#            fprintf(1, 'v2h     #14.7f  #14.7f  #14.7f  \n', v2h(1), v2h(2), v2h(3) );
    #            fprintf(1, 'angles #14.7f   #14.7f  #14.7f \n', theta * rad2deg, theta1 * rad2deg, copa * rad2deg );
    tempg = v2g - v2
    temph = v2h - v2
    magtempg = smu.mag(tempg)
    magtemph = smu.mag(temph)
    #            fprintf(1, 'dg     #14.7f  #14.7f  #14.7f  \n', tempg(1), tempg(2), tempg(3) );
#            fprintf(1, 'dh     #14.7f  #14.7f  #14.7f  \n', temph(1), temph(2), temph(3) );
    print('%3i %11.3f %8.4f  %12.6f  %12.6f %10.5f dg %14.9f %12s dh %14.9f m/s %12s'
          % (i, t2, t2 / period, theta * rad2deg, theta1 * rad2deg,
             copa * rad2deg, magtempg * 1000, errorg, magtemph * 1000, errorh))
    i = i + 1



# 6 Apr 2004 19:02:28.3860     12709.670465397     9059.657416996      4937.914304833    -2.425788263     1.145863933     4.141391183
# 6 Apr 2004 19:03:28.3860     12562.052651693     9126.923951210      5185.576718081    -2.494671410     1.096292760     4.113797579
# 6 Apr 2004 19:04:28.3860     12410.326214883     9191.205376885      5431.543104340    -2.562738634     1.046363027     4.084858491

# 6 Apr 2004 19:08:28.3860     11763.346541426     9418.074295602      6396.848305457    -2.826407207     0.843386744     3.955843131
# 6 Apr 2004 19:14:28.3860     10678.648538971     9665.704898623      7780.535182624    -3.193778683     0.530988666     3.723761880

# closely spaced
r1 = np.array([12709.670465397, 9059.657416996, 4937.914304833])
r2 = np.array([12562.052651693, 9126.92395121, 5185.576718081])
r3 = np.array([12410.326214883, 9191.205376885, 5431.54310434])
print('-------------------- problem close spaced circorbit \n' % ())
print(' r1 ',  (r1))
print(' r2 ',  (r2))
print(' r3 ',  (r3))
v2g, theta, theta1, copa, errorg = obu.gibbs(r1, r2, r3)
print(' v2g km/s \n', (v2g))
print(' theta %11.6f  theta1  %11.6f  copa %11.6f \n' % (theta * rad2deg, theta1 * rad2deg, copa * rad2deg))
# just needs to be in secs - so precision isn't lost
jd1 = 0.0
jd2 = 60.0 / 86400.0

jd3 = 120.0 / 86400.0

# use t1, t2, t3 as secs for greater accuracy
v2h, theta, theta1, copa, errorh = obu.hgibbs(r1, r2, r3, jd1, jd2, jd3)
print(' v2h km/s \n', (v2h))
print(' theta %11.6f  theta1  %11.6f  copa %11.6f \n' % (theta * rad2deg, theta1 * rad2deg, copa * rad2deg))
# wider spaced
r1 = np.array([12709.670465397, 9059.657416996, 4937.914304833])
r2 = np.array([11763.346541426, 9418.074295602, 6396.848305457])
r3 = np.array([10678.648538971, 9665.704898623, 7780.535182624])
print('-------------------- problem wider spaced circorbit \n')
print(' r1 ',  (r1))
print(' r2 ',  (r2))
print(' r3 ',  (r3))
v2g, theta, theta1, copa, errorg = obu.gibbs(r1, r2, r3)
print(' v2g km/s \n', (v2g))
print(' theta %11.6f  theta1  %11.6f  copa %11.6f \n' % (theta * rad2deg, theta1 * rad2deg, copa * rad2deg))
# just needs to be in secs - so precision isn't lost
jd1 = 0.0

jd2 = 360.0 / 86400.0

jd3 = 720.0 / 86400.0

# use t1, t2, t3 as secs for greater accuracy
v2h, theta, theta1, copa, errorh = obu.hgibbs(r1, r2, r3, jd1, jd2, jd3)
print(' v2h km/s \n', (v2h))
print(' theta %11.6f  theta1  %11.6f  copa %11.6f \n' % (theta * rad2deg, theta1 * rad2deg, copa * rad2deg))
