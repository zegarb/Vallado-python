#     -----------------------------------------------------------------
#
#                              prob7_1.m
#
#  this file demonstrates problem 7-1.
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
#             9 oct 07  david vallado
#                         original
#  changes :
#             9 oct 07  david vallado
#                         original baseline
#
#     *****************************************************************

import numpy as np
import math
import orbit_utils as obu
from space_constants import *
import space_conversions as sc
import spacemath_utils as smu
import spacetime_utils as stu


tus = 86400 / au
mu = musun

# 2 McKusky Heliocentric     pg 82
#  1961 11  15  0  0  0.000    0.000000    0.000000    0.00 339.8645833 -10.89736111
#  1961 12   2  0  0  0.000    0.000000    0.000000    0.00 341.9504167 -11.77486111
#  1961 12  12  0  0  0.400    0.000000    0.000000    0.00 343.7545833 -11.93258333
# 1 Gauss Orig - Montenbruck pg 245
#  1805  9   5 24  8 42.000    0.000000    0.000000    0.00  95.9897500  22.35752222
#  1806  1  17 22  5 42.000    0.000000    0.000000    0.00 101.3112083  30.35672222
#  1806  5  23 20 20 56.400    0.000000    0.000000    0.00 121.9358333  28.04640000
# 3 Rockwell  Obj J91A01P
#  1995  3  29  8 36 41.184    0.000000    0.000000    0.00 162.8180833   6.42436111
#  1995  3  29 10  8 24.000    0.000000    0.000000    0.00 162.8065416   6.42838888
#  1995  3  29 10 59 26.880    0.000000    0.000000    0.00 162.8001250   6.43058333

jd1,jd1f = stu.jdayall(1805,9,5,24,16,5.0,'j')
rtasc1 = 95.98975 * deg2rad
decl1 = 22.35752222 * deg2rad
jd2,jd2f = stu.jdayall(1806,1,17,22,9,5.0,'j')
rtasc2 = 101.3112083 * deg2rad
decl2 = 30.35672222 * deg2rad
jd3,jd3f = stu.jdayall(1806,5,23,20,39,9.0,'j')
rtasc3 = 121.9358333 * deg2rad
decl3 = 28.0464 * deg2rad
# earth position in 1805???
rsite1 = np.transpose(np.array([- 0.9527095 * re,0.3014385 * re,0.1308476 * re]))
rsite2 = np.transpose(np.array([0.419654 * re,- 0.8163169 * re,- 0.3543437 * re]))
rsite3 = np.transpose(np.array([0.5039939 * re,0.8058414 * re,0.3497957 * re]))
# ---------------------ceres orbital observations
jd1,jd1f = stu.jdayall(1990,1,1,0,0,0.0,'j')
rtasc1 = sc.hms2rad(5,41,14.1)
decl1 = sc.dms2rad(26,26,0.0)
jd2,jd2f = stu.jdayall(1990,6,12,0,0,0.0,'j')
rtasc2 = sc.hms2rad(7,59,30.5)
decl2 = sc.dms2rad(26,58,15.0)
jd3,jd3f = stu.jdayall(1990,12,25,0,0,0.0,'j')
rtasc3 = sc.hms2rad(13,33,8.7)
decl3 = sc.dms2rad(1,0,48.0)
# earth positions at 1 jan, 12 jun, and 25 dec 1990
rsite1 = np.transpose(np.array([- 0.178961581 * re,0.887451861 * re,0.38475073 * re]))
rsite2 = np.transpose(np.array([- 0.160038443 * re,- 0.919583501 * re,- 0.398758454 * re]))
rsite3 = np.transpose(np.array([- 0.051224491 * re,0.901866361 * re,0.390962686 * re]))
# -------------------------------mars orbital observations
jd1,jd1f = stu.jdayall(1989,12,31,0,0,0.0,'j')
rtasc1 = 16.4429 * 15.0 * np.pi / 180.0
decl1 = - 21.711 * np.pi / 180.0
jd2,jd2f = stu.jdayall(1990,6,9,0,0,0.0,'j')
rtasc2 = 0.4034 * 15.0 * np.pi / 180.0
decl2 = 0.507 * np.pi / 180.0
jd3,jd3f = stu.jdayall(1990,12,26,0,0,0.0,'j')
rtasc3 = 3.6975 * 15.0 * np.pi / 180.0
decl3 = 21.975 * np.pi / 180.0
# earth positions at 1 jan, 12 jun, and 25 dec 1990
rsite1 = np.transpose(np.array([- 0.14376 * re,0.8925 * re,0.38697 * re]))
rsite2 = np.transpose(np.array([- 0.22664 * re,- 0.90768 * re,- 0.39355 * re]))
rsite3 = np.transpose(np.array([- 0.05231 * re,0.9011 * re,0.3907 * re]))
r2ans = np.array([1.08711,- 0.76841,- 0.38185])

# -------------------------------jupiter orbital observations
jd1,jd1f = stu.jdayall(1989,12,31,0,0,0.0,'j')
rtasc1 = 6.4077 * 15.0 * np.pi / 180.0
decl1 = 23.198 * np.pi / 180.0
jd2,jd2f = stu.jdayall(1990,6,9,0,0,0.0,'j')
rtasc2 = 7.045 * 15.0 * np.pi / 180.0
decl2 = 22.834 * np.pi / 180.0
jd3,jd3f = stu.jdayall(1990,12,26,0,0,0.0,'j')
rtasc3 = 9.0252 * 15.0 * np.pi / 180.0
decl3 = 17.625 * np.pi / 180.0
# earth positions at 1 jan, 12 jun, and 25 dec 1990
rsite1 = np.transpose(np.array([- 0.14376 * re,0.8925 * re,0.38697 * re]))
rsite2 = np.transpose(np.array([- 0.22664 * re,- 0.90768 * re,- 0.39355 * re]))
rsite3 = np.transpose(np.array([- 0.05231 * re,0.9011 * re,0.3907 * re]))
r2ans = np.array([- 1.74284,4.49538,1.96939])

# -------------------------------jupiter orbital observations
#     jd1=stu.jdayall(  2011,  1,   1,  0,  0,  0.000,'j');
#     rtasc1 =   sc.hms2rad(23, 49, 24.0 );
#     decl1  =   sc.dms2rad(-2, 32,  0.0 );

#     jd2=stu.jdayall(  2011,  6,  10,  0,  0,  0.000,'j');
#     rtasc2 =   sc.hms2rad( 1, 57, 24.0 );
#     decl2  =   sc.dms2rad(10, 48,  0.0 );

#     jd3=stu.jdayall(  2011,  12, 19,  0,  0,  0.000,'j');
#     rtasc3 =   sc.hms2rad( 1, 55, 12.0 );
#     decl3  =   sc.dms2rad(10, 25,  0.0 );

#     # earth positions at 1 jan, 12 jun, and 25 dec 1990
#     rsite1 = [-0.178961581*re  0.887451861*re  0.384750730*re]';
#     rsite2 = [-0.160038443*re -0.919583501*re -0.398758454*re]';
#     rsite3 = [-0.051224491*re  0.901866361*re  0.390962686*re]';
#    wrong
#      r2ans = [-1.74284  4.49538  1.96939];  # jupiter, au

#    [r2,v2] = anglesl ( decl1,decl2,decl3,rtasc1,rtasc2,rtasc3,jd1,jd2,jd3, rsite1,rsite2,rsite3, re, mu, tu );
#    processtype = 'anglesl';
r2,v2 = obu.anglesg(decl1,decl2,decl3,rtasc1,rtasc2,rtasc3,jd1, jd1f,jd2, jd2f,jd3, jd3f,rsite1,rsite2,rsite3)
print('r2 and v2 check:')
print(r2)
print(v2)
processtype = 'anglesg'

#    [r2,v2] = anglesdr ( decl1,decl2,decl3,rtasc1,rtasc2,rtasc3,jd1,jd2,jd3, rsite1,rsite2,rsite3, re, mu, tu );
#    processtype = 'anglesdr';
#    [r2,v2] = anglesgood ( decl1,decl2,decl3,rtasc1,rtasc2,rtasc3,jd1,jd2,jd3, rsite1,rsite2,rsite3 );
#    processtype = 'anglesgood';

# -------------- write out answer --------------
print('\n\ninputs: \n\n' % ())
latgc,latgd,lon,alt = sc.ecef2ll(rsite1)
print('Site obs1 %11.7f %11.7f %11.7f au  lat %11.7f lon %11.7f alt %11.7f  \n' % (rsite1[0],rsite1[1],rsite1[2], 
                                                                                   latgc, lon, alt ))
latgc,latgd,lon,alt = sc.ecef2ll(rsite2)
print('Site obs2 %11.7f %11.7f %11.7f au  lat %11.7f lon %11.7f alt %11.7f  \n' % (rsite2[0],rsite2[1],rsite2[2], 
                                                                                   latgc, lon, alt ))
latgc,latgd,lon,alt = sc.ecef2ll(rsite3)
print('Site obs3 %11.7f %11.7f %11.7f au  lat %11.7f lon %11.7f alt %11.7f  \n' % (rsite3[0],rsite3[1],rsite3[2], 
                                                                                   latgc, lon, alt ))
year,mon,day,hr,min,sec = stu.invjday(jd1,jd1f)
print('Obs#1 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n' % (year,mon,day,hr,min,sec,rtasc1 * rad2deg,decl1 * rad2deg))
year,mon,day,hr,min,sec = stu.invjday(jd2,jd2f)
print('Obs#2 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n' % (year,mon,day,hr,min,sec,rtasc2 * rad2deg,decl2 * rad2deg))
year,mon,day,hr,min,sec = stu.invjday(jd3,jd3f)
print('Obs#3 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n' % (year,mon,day,hr,min,sec,rtasc3 * rad2deg,decl3 * rad2deg))
print('\nsolution by %s \n\n' % (processtype))
print('r2     %11.7f   %11.7f  %11.7f au    %11.7f  %11.7f  %11.7f km \n' % (r2[0] / re, r2[1] / re, r2[2] / re , 
                                                                             r2[0] * au, r2[1] * au, r2[2] * au))
print('r2 ans %11.7f   %11.7f  %11.7f au    %11.7f  %11.7f  %11.7f km \n' % (r2ans[0] / re, r2ans[1] / re, r2ans[2] / re , 
                                                                             r2ans[0] * au, r2ans[1] * au, r2ans[2] * au))
print('v2     %11.7f   %11.7f  %11.7f au/tu %11.7f  %11.7f  %11.7f km/s\n' % (v2[0] / velkmps, v2[1] / velkmps, v2[2] / velkmps, 
                                                                              v2[0] * tus, v2[1] * tus, v2[2] * tus))
#    fprintf(1,'v2 ans #11.7f   #11.7f  #11.7f au/tu #11.7f  #11.7f  #11.7f km/s\n',v2ans/velkmps, v2ans);

p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper = sc.rv2coe(r2,v2,mu)

print('         p au          a au         ecc       incl deg     raan deg    argp deg     nu deg      m deg  \n' % ())
print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f \n' % (p,a,ecc,incl * rad2deg,omega * rad2deg,argp * rad2deg,nu * rad2deg,m * rad2deg))
#    [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe(r2ans,v2ans, mu);
#    fprintf(1,'         p au          a au         ecc       incl deg     raan deg    argp deg     nu deg      m deg  \n');
#    fprintf(1,'coes #11.4f #11.4f #13.9f #13.7f #11.5f #11.5f #11.5f #11.5f \n',...
#              p,a,ecc,incl*rad,omega*rad,argp*rad,nu*rad,m*rad );

rad = 180 / np.pi
au = 149597870.0

mu = 132712428000.0
e = np.array([3760.5,0.0,0.0,102.993,0.999988,0.985627,0.016676,124.167])
x = np.array([3760.5,1.3036,100.561,14.816,5.20192,0.0831123,0.048942,218.527])
incl = x(2) * deg2rad
omega = x(3) * deg2rad
lonper = x(4) * deg2rad
a = x(5) * au
n = x(6)
ecc = x(7)
truelon = x(8) * deg2rad
m = truelon - lonper
argp = lonper - omega
p = a * (1.0 - ecc ** 2)
e0,nu = obu.newtonm(ecc,m)
r2,v2 = sc.coe2rv(p,ecc,incl,omega,argp,nu,arglat,truelon,lonper,mu)
print('r2     %11.7f   %11.7f  %11.7f au    %11.7f  %11.7f  %11.7f km \n' % (r2[0] / au, r2[1] / au, r[2] / au,
                                                                             r2[0], r2[1], r2[2]))
print('v2     %11.7f   %11.7f  %11.7f au/tu %11.7f  %11.7f  %11.7f km/s\n' % (v2[0] * tus, v2[1] * tus, v2[2] * tus,
                                                                              v2[0], v2[1], v2[2]))
