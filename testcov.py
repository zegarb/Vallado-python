#
#     -----------------------------------------------------------------
#
#                               testcov
#
#  this file tests the accuracy of the covariance functions.
#
#                          companion code for
#             fundamentals of astrodyanmics and applications
#                                2013
#                          by david vallado
#
#     (h)               email davallado@gmail.com
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#            22 jun 15  david vallado
#                         many fixes for new paper
#  changes :
#             9 aug 03  david vallado
#                         fix units on n in equinoc?
#            25 jul 03  david vallado
#                         fixes for tm
#            20 jul 03  david vallado
#                         fixes to add extras for paper
#            25 jun 03  david vallado
#                         update for alternate cases
#            26 may 03  david vallado
#                         fix units and values
#            28 oct 02  david vallado
#                         fix covariance transformations
#             5 sep 02  david vallado
#                         fix cov trace input, misc
#            26 aug 02  david vallado
#                         work on partials for covariance and test cases
#            17 aug 02  david vallado
#                         work on partials for covariance
#            30 jun 02  david vallado
#                         breakout reduction, time, 2body
#            24 may 02  david vallado
#                         original baseline
#
#     *****************************************************************


import numpy as np
import orbit_utils as obu
from space_constants import au, rad2deg, re, velkmps, mu, twopi, tusec
import space_conversions as sc
import spacemath_utils as smu
import spacetime_utils as stu

small = 1e-18
rad2 = rad2deg * rad2deg
anom = 'truea'
anom = 'truen'
#anom = 'meana';
anomflt = 'latlon'
testnum = 0

ddpsi = 0.0
ddeps = 0.0
# --------------- 0. Sal test
if testnum == 0:
    reci = np.array([11074.95274,40629.74421,- 32.1123199])
    veci = np.array([- 2.940822436,0.9007122363,0.002036330819])
    aeci = np.array([0.001,0.002,0.003])
    year = 2022
    mon = 106 #Is this right? -mjc
    day = 28
    hr = 15
    min = 8
    sec = 51.655
    dut1 = 0.16236
    dat = 37
    xp = 0.0987
    yp = 0.286
    lod = 0.0
    timezone = 0
    order = 4 #Should order be 106 like the others? -mjc
    terms = 2
    cartcov = np.array([[24097166,86695628,- 5509927,- 6294.97,1752.326,17.65861],
                        [86695628,453000000.0,- 28000000.0,- 32967.4,6319.431,90.73355],
                        [- 5509927,- 28000000.0,1771703,2061.582,- 401.582,- 5.67764],
                        [- 6294.97,- 32967.4,2061.582,6949865,- 1352586,0.385006],
                        [1752.326,6319.431,- 401.582,- 1352586,263241.3,2.013476],
                        [17.65861,90.73355,- 5.67764,0.385006,2.013476,33.37338]])

# ------- navy test
if testnum == 1:
    reci = np.array([- 605.7922166, - 5870.22951108, 3493.05319896])
    veci = np.array([- 1.56825429, - 3.70234891, - 6.47948395])
    aeci = np.array([0.001, 0.002, 0.003])
    year = 2000
    mon = 12
    day = 15
    hr = 16
    min = 58
    sec = 50.208
    dut1 = 0.10597
    dat = 32
    xp = 0.0
    yp = 0.0
    lod = 0.0
    terms = 2
    timezone = 0
    order = 106
    cartcov = np.array([[8.04204e-13,- 7.41923e-13,2.02623e-12,- 3.78793e-16,1.77302e-13,- 2.48352e-13],
                        [- 7.41923e-13,2.19044e-12,3.82034e-12,- 1.1901e-15,2.83844e-13,- 3.67925e-13],
                        [2.02623e-12,3.82034e-12,1.9757e-11,- 9.35438e-15,- 1.47386e-12,3.01542e-12],
                        [- 3.78793e-16,- 1.1901e-15,- 9.35438e-15,1.56236e-17,8.86093e-17,4.56974e-16],
                        [1.77302e-13,2.83844e-13,- 1.47386e-12,8.86093e-17,9.67797e-13,7.23072e-13],
                        [- 2.48352e-13,- 3.67925e-13,3.01542e-12,4.56974e-16,7.23072e-13,1.84123e-12]])
    #  3.06012e-12  1.16922e-11  6.72154e-11  1.28647e-13  4.11455e-13 -3.89727e-12  2.20508e-05

# ------- accuracy sensitivity test
if testnum == 2:
    reci = np.array([- 605.7922166, - 5870.22951108, 3493.05319896])
    veci = np.array([- 1.56825429, - 3.70234891, - 6.47948395])
    aeci = np.array([0.001, 0.002, 0.003])
    year = 2000
    mon = 12
    day = 15
    hr = 16
    min = 58
    sec = 50.208
    dut1 = 0.10597
    dat = 32
    xp = 0.0
    yp = 0.0
    lod = 0.0
    timezone = 0
    terms = 2
    order = 106
    covntw = np.array([[1.0,0.0,0.0,0.0,0.0,0.0],
                       [0.0,10.0,0.0,0.0,0.0,0.0],
                       [0.0,0.0,1.0,0.0,0.0,0.0],
                       [0.0,0.0,0.0,1e-06,0.0,0.0],
                       [0.0,0.0,0.0,0.0,0.0001,0.0],
                       [0.0,0.0,0.0,0.0,0.0,1e-06]])
    #, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
    ut1,tut1,jdut1,jdut1frac,utc,tai,tt,ttt,jdtt,jdttfrac,tdb,ttdb,jdtdb,jdtdbfrac = \
        stu.convtime(year,mon,day,hr,min,sec,timezone,dut1,dat)
    cartstate,classstate,flstate,eqstate = \
        sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps) #latlon or radec? -mjc
    sc.printcov(covntw,'ct','m',anom)
    cartcov,tm = sc.covntw2ct(covntw,cartstate)
    sc.printcov(cartcov,'ct','m',anom)

# ------- paper test
if testnum == 3:
    reci = np.array([- 605.7922166, - 5870.22951108, 3493.05319896])
    veci = np.array([- 1.56825429, - 3.70234891, - 6.47948395])
    aeci = np.array([0.001, 0.002, 0.003])
    year = 2000
    mon = 12
    day = 15
    hr = 16
    min = 58
    sec = 50.208
    dut1 = 0.10597
    dat = 32
    xp = 0.0
    yp = 0.0
    lod = 0.0
    timezone = 0
    terms = 2
    order = 106
    cartcov = np.array([[100.0,0.01,0.01,0.0001,0.0001,0.0001],
                        [0.01,100.0,0.01,0.0001,0.0001,0.0001],
                        [0.01,0.01,100.0,0.0001,0.0001,0.0001],
                        [0.0001,0.0001,0.0001,0.0001,1e-06,1e-06],
                        [0.0001,0.0001,0.0001,1e-06,0.0001,1e-06],
                        [0.0001,0.0001,0.0001,1e-06,1e-06,0.0001]])
    classcovtrace = np.array([[6.489403,2.979136e-07,2.117582e-22,5.293956e-22,0.03465237,0.009048777],
                              [2.979136e-07,4.353986e-14,6.708473e-27,1.091152e-26,1.059842e-09,2.77089e-10],
                              [2.117582e-22,6.708473e-27,5.652526e-11,6.541195e-14,8.709806e-15,- 4.430849e-23],
                              [5.293956e-22,1.091152e-26,6.541195e-14,5.749008e-11,7.654984e-12,- 7.072071e-23],
                              [0.03465237,1.059842e-09,8.709806e-15,7.654984e-12,0.0002232401,5.830381e-05],
                              [0.009048777,2.77089e-10,- 4.430849e-23,- 7.072071e-23,5.830381e-05,1.522727e-05]])
    eqcovtrace = np.array([[5.70396e-14,2.56241e-14,2.218139e-12,6.530893e-15,- 1.258488e-18,6.777803e-18],
                           [2.56241e-14,6.348393e-14,- 1.958979e-12,7.398763e-15,- 4.223513e-18,2.274645e-17],
                           [2.218139e-12,- 1.958979e-12,3.562517e-10,2.727333e-16,2.373188e-13,- 1.278121e-12],
                           [6.530893e-15,7.398763e-15,2.727333e-16,1.256922e-15,4.930381e-32,- 1.972152e-31],
                           [- 1.258488e-18,- 4.223513e-18,2.373188e-13,4.930381e-32,2.292326e-14,- 2.061165e-17],
                           [6.777803e-18,2.274645e-17,- 1.278121e-12,- 1.972152e-31,- 2.061165e-17,2.288388e-14]])
    sc.printcov(cartcov,'ct','m',anom)
    print('Check the xxxx (mat*transpose) of cartcov \n' % ())
    print((np.transpose((cartcov * np.transpose(cartcov)))))
    #            printcov( classcovtrace,'cl','t',anom );
#            classcovt = covunits( classcovtrace,anom,'cl','m' );
#            printcov( classcovt,'cl','m',anom );
#            printcov( eqcovtrace,'eq','t',anom );
#            eqcovt = covunits( eqcovtrace,anom,'eq','m' );
#            printcov( eqcovt,'eq','m',anom );

if testnum == 4:
    reci = np.array([4364.51524926, 4748.1760294, 2430.20427647])
    veci = np.array([5.87962414, - 4.10294944, - 2.53527819])
    year = 2000
    mon = 12
    day = 14
    hr = 5
    min = 25
    sec = 3.461
    dut1 = 0.10597
    dat = 32
    xp = 0.0
    yp = 0.0
    lod = 0.0
    terms = 2
    timezone = 0
    order = 106
    # ----- in from usstratcom is lower diagonal!!!
    eqcov = np.array([[4.68914e-11,1.6009e-11,1.64731e-10,- 4.38141e-16,- 1.41195e-10,1.09999e-11,- 8.49933e-12],
                      [1.6009e-11,1.06881e-11,6.39732e-11,6.53124e-16,- 5.60554e-11,2.56099e-11,- 4.4424e-13],
                      [1.64731e-10,6.39732e-11,6.71336e-10,1.04455e-14,- 6.37507e-10,9.40554e-11,7.34847e-12],
                      [- 4.38141e-16,6.53124e-16,1.04455e-14,2.10897e-17,1.59272e-14,1.03081e-14,1.80744e-13],
                      [- 1.41195e-10,- 5.60554e-11,- 6.37507e-10,1.59272e-14,7.59844e-10,- 1.05437e-10,1.77857e-10],
                      [1.09999e-11,2.56099e-11,9.40554e-11,1.03081e-14,- 1.05437e-10,1.86472e-10,3.99631e-11],
                      [- 8.49933e-12,- 4.4424e-13,7.34847e-12,1.80744e-13,1.77857e-10,3.99631e-11,4.67207e-06]])
    sc.printcov(eqcovtrace,'eq','t',anom)
    #eqcovt = covunits(eqcovtrace,anom,'eq','m')
    #sc.printcov(eqcovt,'eq','m',anom)

if testnum == 5:
    reci = np.array([10127.26750234, 6972.89492052, 4902.05501566])
    veci = np.array([- 3.38546947, 0.81341143, - 5.70012872])
    year = 2000
    mon = 12
    day = 1
    hr = 4
    min = 39
    sec = 8.735
    dut1 = 0.11974
    dat = 32
    xp = 0.0
    yp = 0.0
    lod = 0.0
    terms = 2
    timezone = 0
    order = 106
    eqcovtrace = np.array([[1.79533e-09,9.35041e-10,- 1.99413e-09,- 1.63427e-13,1.05541e-10,- 5.77137e-10,5.44105e-10],
                           [9.35041e-10,1.04209e-09,- 5.91847e-10,- 9.78886e-14,- 4.0233e-11,- 4.39526e-10,6.10449e-10],
                           [- 1.99413e-09,- 5.91847e-10,3.63062e-09,2.62623e-13,- 3.98481e-10,9.56574e-10,2.86557e-09],
                           [- 1.63427e-13,- 9.78886e-14,2.62623e-13,5.56451e-16,5.18072e-14,- 2.60304e-14,- 2.27216e-12],
                           [1.05541e-10,- 4.0233e-11,- 3.98481e-10,5.18072e-14,1.00164e-09,2.40555e-10,- 1.10189e-10],
                           [- 5.77137e-10,- 4.39526e-10,9.56574e-10,- 2.60304e-14,2.40555e-10,2.05855e-09,- 1.97222e-11],
                           [5.44105e-10,6.10449e-10,2.86557e-09,- 2.27216e-12,- 1.10189e-10,- 1.97222e-11,1.13395e-05]])
    sc.printcov(eqcovtrace,'eq','t',anom)
    #eqcovt = covunits(eqcovtrace,anom,'eq','m')
    #sc.printcov(eqcovt,'eq','m',anom)

#         a = 6860.7631;
#         ecc = 0.0010640;
#         p = a*(1.0 - ecc^2);
#         incl = 90.0/rad; #97.65184/rad;
#         omega = 79.54701/rad;
#         argp = 83.86041/rad;
#         nu = 65.21303/rad;
#         m = 65.10238/rad;
#         arglat = 0.0;
#         truelon = 0.0;
#         lonper = 0.0;
#         [reci,veci] = coe2rv ( p,ecc,incl,omega,argp,nu,arglat,truelon,lonper );

print(' ---------------------------- begin tests ------------------------- \n' % ())
print(' ---------------------------- begin tests ------------------------- ' % ())
anomeq1 = 'mean'

anomeq1 = 'true'

anomeq2 = 'n'

anomflt = 'latlon'

#, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
ut1,tut1,jdut1,jdut1frac,utc,tai,tt,ttt,jdtt,jdttfrac,tdb,ttdb,jdtdb,jdtdbfrac = stu.convtime(year,mon,day,hr,min,sec,timezone,dut1,dat)
# --- convert the eci state into the various other state formats (classical, equinoctial, etc)
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anomeq1+anomeq2,anomflt,ddpsi,ddeps)

print('==================== do the sensitivity tests \n' % ())
print('1.  Cartesian Covariance \n' % ())
sc.printcov(cartcov,'ct','m',anomeq1+'a')
# test partials
#         # partial a wrt rx
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         magr = mag(reci);
#         delta = reci(1) * 0.00001;
#         reci(1) = reci(1) + delta;
#         magr1 = mag(reci);
#         [p,a1,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         p0 = (a-a1)/delta;
#         p1 = 2.0*a^2*reci(1) / magr^3;
#         fprintf(1,' a wrt rx #14.14f  #14.14f \n', p0, p1);

#         # partial n wrt rx
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         magr = mag(reci);
#         n = sqrt(mu/a^3);
#         recin = reci;
#         delta = recin(1) * 0.00001;
#         recin(1) = recin(1) + delta;
#         magr1 = mag(recin);
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (recin,veci);
#         n1 = sqrt(mu/a^3);
#         p0 = (n-n1)/delta;
#         p2 = -3.0*n1*a*reci(1)/magr1^3;
#         fprintf(1,' n wrt rx #14.14f  #14.14f \n', p0, p2);

#         # partial a wrt vx
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         magr = mag(veci);
#         vecin = veci;
#         delta = vecin(1) * 0.00001;
#         vecin(1) = vecin(1) + delta;
#         magv1 = mag(vecin);
#         [p,a1,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,vecin);
#         p0 = (a-a1)/delta;
#         p1 = 2.0*veci(1) / (n^2*a);
#         fprintf(1,' a wrt vx #14.14f  #14.14f \n', p0, p1);

#         # partial n wrt vx
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         magr = mag(veci);
#         n = sqrt(mu/a^3);
#         vecin = veci;
#         delta = vecin(1) * 0.00001;
#         vecin(1) = vecin(1) + delta;
#         magv1 = mag(vecin);
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,vecin);
#         n1 = sqrt(mu/a^3);
#         p0 = (n-n1)/delta;
#         p2 = -3.0*vecin(1)/(n*a^2);
#         fprintf(1,' n wrt vx #14.14f  #14.14f \n', p0, p2);

#         # partial n wrt vz
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         magr = mag(veci);
#         n = sqrt(mu/a^3);
#         vecin = veci;
#         delta = vecin(3) * 0.00001;
#         vecin(3) = vecin(3) + delta;
#         magv1 = mag(vecin);
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,vecin);
#         n1 = sqrt(mu/a^3);
#         p0 = (n-n1)/delta;
#         p2 = -3.0*vecin(3)/(n*a^2);
#         fprintf(1,' n wrt vz #14.14f  #14.14f \n', p0, p2);

#         # partial rx wrt a
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         delta = a * 0.00001;
#         a = a + delta;
#         p = a*(1-ecc^2);
#         [reci1, veci1] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);
#         p0 = (reci(1)-reci1(1))/delta;
#         p1 = reci(1) / a;
#         fprintf(1,' rx wrt a #14.14f  #14.14f \n', p0, p1);
#         # partial rx wrt n
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         n = sqrt(mu/a^3);
#         delta = n * 0.00001;
#         n = n + delta;
#         a = (mu/n^2)^(1/3);
#         p = a*(1-ecc^2);
#         [reci1, veci1] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);
#         p0 = (reci(1)-reci1(1))/delta;
#         p1 = -2*reci(1) / (3*n);
#         fprintf(1,' rx wrt n #14.14f  #14.14f \n', p0, p1);

#         # partial vx wrt a
#          [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#          delta = a * 0.00001;
#          a = a + delta;
#          p = a*(1-ecc^2);
#          [reci1, veci1] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);
#          p0 = (veci(1)-veci1(1))/delta;
#          p1 = -veci(1) / (2*a);
#          fprintf(1,' vx wrt a #14.14f  #14.14f \n', p0, p1);
#          # partial vx wrt n
#          [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#          n = sqrt(mu/a^3);
#          delta = n * 0.00001;
#          n = n + delta;
#          a = (mu/n^2)^(1/3);
#          p = a*(1-ecc^2);
#          [reci1, veci1] = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper);
#          p0 = (veci(1)-veci1(1))/delta;
#          p1 = veci(1) / (3*n);
#          fprintf(1,' vx wrt n #14.14f  #14.14f \n', p0, p1);

#         # partial vz wrt a
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         delta = a * 0.00001;
#         a1 = a + delta;
#         p = a1*(1-ecc^2);
#         [recin,vecin] = coe2rv ( p,ecc,incl,omega,argp,nu,0.0,0.0,0.0 );
#         n1 = sqrt(mu/a1^3);
#         p0 = (veci(3)-vecin(3))/delta;
#         p2 = -veci(3)/(2*a);
#         fprintf(1,' vz wrt a #14.14f  #14.14f \n', p0, p2);
#         # partial vz wrt n
#         [p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         n = sqrt(mu/a^3);
#         delta = n * 0.00001;
#         n1 = n + delta;
#         a1 = (mu/n1^2)^(1/3);
#         p = a1*(1-ecc^2);
#         [recin,vecin] = coe2rv ( p,ecc,incl,omega,argp,nu,0.0,0.0,0.0 );
#         p0 = (veci(3)-vecin(3))/delta;
#         p2 = veci(3)/(3*n);
#         fprintf(1,' vz wrt n #14.14f  #14.14f \n', p0, p2);
#         # partial lat wrt rx
#         #[p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe  (reci,veci);
#         delta = flstate(2) * 0.0000001;
#         latgc1 = flstate(2) + delta;
#         magr = flstate(5);
#             recef(1) = magr*0.001*cos(latgc1)*cos(flstate(1));  # in km
#             recef(2) = magr*0.001*cos(latgc1)*sin(flstate(1));
#             recef(3) = magr*0.001*sin(latgc1);
#             vecef = [0; 0; 0];
#             aecef = [0;0;0];
#             [recin,vecin,a] = ecef2eci(recef',vecef,aecef,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps);

#         n = sqrt(mu/a^3);
#         a1 = (mu/n1^2)^(1/3);
#         p = a1*(1-ecc^2);
#         [recin,vecin] = coe2rv ( p,ecc,incl,omega,argp,nu,0.0,0.0,0.0 );
#         p0 = (reci(1)-recin(1))/delta;
#         p1 = -magr*sin(flstate(1))*cos(flstate(2));
#         p2 = -reci(2)/sqrt(reci(1)^2 + reci(2)^2);
#         fprintf(1,' lat wrt rx #14.14f  #14.14f \n', p0, p2);
#  pause

# paper approach
# ===============================================================================================

print('===============================================================================================\n' % ())
print('1.  RSW and NTW Covariance from Cartesian #1 above ------------------- \n' % ())
covrsw,tm = sc.covct2rsw(cartcov,cartstate)
print('rsw\n' % ())
sc.printcov(covrsw,'ct','m',anom)
temt = covrsw
covntw,tm = sc.covct2ntw(cartcov,cartstate)
print('\nntw\n' % ())
sc.printcov(covntw,'ct','m',anom)
temt = covntw

print('===============================================================================================\n' % ())
# --- convert the eci state into the various other state formats (classical, equinoctial, etc)
anom = 'meana'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
print('\n2.  Classical Covariance from Cartesian #1 above (meana) ------------------- \n' % ())
classcovmeana,tmct2cl = sc.covct2cl(cartcov,cartstate,anom)
sc.printcov(classcovmeana,'cl','m',anom)
print(' \nCartesian Covariance from Classical #2 above \n' % ())
cartcovmeanarev,tmcl2ct = sc.covcl2ct(classcovmeana,classstate,anom)
sc.printcov(cartcovmeanarev,'ct','m',anom)
smu.printdiff(' cartcov - cartcovmeanarev \n',cartcov,cartcovmeanarev)
sc.printcov((tmct2cl @ tmcl2ct),'tm','m',anom)
#sc.printcov((tmcl2ct @ tmct2cl),'tm','m',anom)
#smu.printdiff( ' tmct2cl - np.linalg.inv(tmcl2ct) \n', tmct2cl, np.linalg.inv(tmcl2ct));
#rintdiff( ' tmcl2ct - np.linalg.inv(tmct2cl) \n', tmcl2ct, np.linalg.inv(tmct2cl));

# ===============================================================================================
print('\n2.  Classical Covariance from Cartesian #1 above (meann) ------------------- \n' % ())
anom = 'meann'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
classcovmeann,tmct2cl = sc.covct2cl(cartcov,cartstate,anom)
sc.printcov(classcovmeann,'cl','m',anom)
print('\n  Cartesian Covariance from Classical #2 above \n' % ())
cartcovmeannrev,tmcl2ct = sc.covcl2ct(classcovmeann,classstate,anom)
sc.printcov(cartcovmeannrev,'ct','m',anom)
print('\n' % ())
#         fprintf(1,'-------- tm cl2ct new ---------\n');
#         printcov( tmcl2ct,'tm','m','meann' );
#         fprintf(1,'-------- tm ct2cl new ---------\n');
#         printcov( tmct2cl,'tm','m','meann' );
#         fprintf(1,'\n');

smu.printdiff(' cartcov - cartcovmeannrev \n',cartcov,cartcovmeannrev)
sc.printcov(tmct2cl @ tmcl2ct,'tm','m',anom)
#smu.printdiff( ' tmct2cl - np.linalg.inv(tmcl2ct) \n', tmct2cl, np.linalg.inv(tmcl2ct));
#smu.printdiff( ' tmcl2ct - np.linalg.inv(tmct2cl) \n', tmcl2ct, np.linalg.inv(tmct2cl));


# ===============================================================================================
print('\n2.  Classical Covariance from Cartesian #1 above (truea) -------------------- \n' % ())
anom = 'truea'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
classcovtruea,tmct2cl = sc.covct2cl(cartcov,cartstate,anom)
sc.printcov(classcovtruea,'cl','m',anom)
print('\n  Cartesian Covariance from Classical #2 above \n' % ())
cartcovtruearev,tmcl2ct = sc.covcl2ct(classcovtruea,classstate,anom)
sc.printcov(cartcovtruearev,'ct','m',anom)
print('\n' % ())
tmcl2cttruea = tmcl2ct
smu.printdiff(' cartcov - cartcovtruearev \n',cartcov,cartcovtruearev)
sc.printcov(tmct2cl @ tmcl2ct,'tm','m',anom)
#smu.printdiff( ' tmct2cl - np.linalg.inv(tmcl2ct) \n', tmct2cl, np.linalg.inv(tmcl2ct));
#smu.printdiff( ' tmcl2ct - np.linalg.inv(tmct2cl) \n', tmcl2ct, np.linalg.inv(tmct2cl));


# ===============================================================================================
print('\n2.  Classical Covariance from Cartesian #1 above (truen) -------------------- \n' % ())
anom = 'truen'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
classcovtruen,tmct2cl = sc.covct2cl(cartcov,cartstate,anom)
sc.printcov(classcovtruen,'cl','m',anom)
print('  Cartesian Covariance from Classical #2 above \n' % ())
cartcovtruenrev,tmcl2ct = sc.covcl2ct(classcovtruen,classstate,anom)
sc.printcov(cartcovtruenrev,'ct','m',anom)
print('\n' % ())
tmcl2cttruen = tmcl2ct
smu.printdiff(' cartcov - cartcovtruenrev \n',cartcov,cartcovtruenrev)
sc.printcov(tmct2cl @ tmcl2ct,'tm','m',anom)
#smu.printdiff( ' tmct2cl - np.linalg.inv(tmcl2ct) \n', tmct2cl, np.linalg.inv(tmcl2ct));
#smu.printdiff( ' tmcl2ct - np.linalg.inv(tmct2cl) \n', tmcl2ct, np.linalg.inv(tmct2cl));


# ===============================================================================================
print('===============================================================================================\n' % ())
# strcat(anomeq1,anomeq2)
print('\n3.  Equinoctial Covariance from Classical #2 above (truea) \n' % ())
anom = 'truea'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
classcovtruea,tmct2cl = sc.covct2cl(cartcov,cartstate,anom)
eqcovtruea,tmcl2eq = sc.covcl2eq(classcovtruea,classstate,anom,fr)
sc.printcov(eqcovtruea,'eq','m',anom)
#eqcov = eqcovtruea; # save for later

tmct2eqtruea = tmct2cl @ tmcl2eq
print('\n4.  Classical Covariance from Equinoctial #3 above \n' % ())
classcovtruearev,tmeq2cl = sc.coveq2cl(eqcovtruea,eqstate,anom,fr)
sc.printcov(classcovtruearev,'cl','m',anom)
tmeq2cttruea = tmeq2cl @ tmcl2cttruea
smu.printdiff(' classcovtruea - classcovtruearev \n',classcovtruea,classcovtruearev)
sc.printcov(tmcl2eq @ tmeq2cl,'tm','m',anom)
#smu.printdiff( ' tmcl2eq - np.linalg.inv(tmeq2cl) \n', tmcl2eq, np.linalg.inv(tmeq2cl));
#smu.printdiff( ' tmeq2cl - np.linalg.inv(tmcl2eq) \n', tmeq2cl, np.linalg.inv(tmcl2eq));


# ===============================================================================================
# strcat(anomeq1,anomeq2)
print('\n3.  Equinoctial Covariance from Classical #2 above (truen) \n' % ())
anom = 'truen'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
classcovtruen,tmct2cl = sc.covct2cl(cartcov,cartstate,anom)
eqcovtruen,tmcl2eq = sc.covcl2eq(classcovtruen,classstate,anom,fr)
sc.printcov(eqcovtruen,'eq','m',anom)
#eqcov = eqcovtruen; # save for later

tmct2eqtruen = tmcl2eq @ tmct2cl
print('\n4.  Classical Covariance from Equinoctial #3 above \n' % ())
classcovtruenrev,tmeq2cl = sc.coveq2cl(eqcovtruen,eqstate,anom,fr)
sc.printcov(classcovtruenrev,'cl','m',anom)
tmeq2cttruen = tmcl2cttruen @ tmeq2cl
smu.printdiff(' classcovtruen - classcovtruenrev \n',classcovtruen,classcovtruenrev)
sc.printcov(tmcl2eq @ tmeq2cl,'tm','m',anom)
#smu.printdiff( ' tmcl2eq - np.linalg.inv(tmeq2cl) \n', tmcl2eq, np.linalg.inv(tmeq2cl));
#smu.printdiff( ' tmeq2cl - np.linalg.inv(tmcl2eq) \n', tmeq2cl, np.linalg.inv(tmcl2eq));


# ===============================================================================================
# strcat(anomeq1,anomeq2)
print('\n3.  Equinoctial Covariance from Classical #2 above (meana) \n' % ())
anom = 'meana'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
classcovmeana,tmct2cl = sc.covct2cl(cartcov,cartstate,anom)
eqcovmeana,tmcl2eq = sc.covcl2eq(classcovmeana,classstate,anom,fr)
sc.printcov(eqcovmeana,'eq','m',anom)
print('\n4.  Classical Covariance from Equinoctial #3 above \n' % ())
classcovmeanarev,tmeq2cl = sc.coveq2cl(eqcovmeana,eqstate,anom,fr)
sc.printcov(classcovmeanarev,'cl','m',anom)
smu.printdiff(' classcovmeana - classcovmeanarev \n',classcovmeana,classcovmeanarev)
sc.printcov(tmcl2eq @ tmeq2cl,'tm','m',anom)
#smu.printdiff( ' tmcl2eq - np.linalg.inv(tmeq2cl) \n', tmcl2eq, np.linalg.inv(tmeq2cl));
#smu.printdiff( ' tmeq2cl - np.linalg.inv(tmcl2eq) \n', tmeq2cl, np.linalg.inv(tmcl2eq));


# ===============================================================================================
# strcat(anomeq1,anomeq2)
print('\n3.  Equinoctial Covariance from Classical #2 above (meann) \n' % ())
anom = 'meann'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
classcovmeann,tmct2cl = sc.covct2cl(cartcov,cartstate,anom)
eqcovmeann,tmcl2eq = sc.covcl2eq(classcovmeann,classstate,anom,fr)
sc.printcov(eqcovmeann,'eq','m',anom)
print('\n4.  Classical Covariance from Equinoctial #3 above \n' % ())
classcovmeannrev,tmeq2cl = sc.coveq2cl(eqcovmeann,eqstate,anom,fr)
sc.printcov(classcovmeannrev,'cl','m',anom)
smu.printdiff(' classcovmeann - classcovmeannrev \n',classcovmeann,classcovmeannrev)
sc.printcov(tmcl2eq @ tmeq2cl,'tm','m',anom)
#smu.printdiff( ' tmcl2eq - np.linalg.inv(tmeq2cl) \n', tmcl2eq, np.linalg.inv(tmeq2cl));
#smu.printdiff( ' tmeq2cl - np.linalg.inv(tmcl2eq) \n', tmeq2cl, np.linalg.inv(tmcl2eq));


# ===============================================================================================
print('===============================================================================================\n' % ())
print('\n5.  Equinoctial Covariance from Cartesian #1 above (truen) \n' % ())
anom = 'truen'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
#Problem: Entire bottom row missing -mjc
eqcovDtruen,tmct2eq = sc.covct2eq(cartcov,cartstate,anom,fr)
sc.printcov(eqcovDtruen,'eq','m',anom)
smu.printdiff(' eqcovtruen - eqcovDtruen \n',eqcovtruen,eqcovDtruen)
print('\n6.  Cartesian Covariance from Equinoctial #5 above \n' % ())
# try with intermediate one
eqcovDtruenx = np.array([4.306970722e-17,1.2444284692e-14,1.4014592713e-14,- 4.7096511176e-17,- 5.1804415247e-17,2.6405645277e-17,1.2444284692e-14,6.1777913966e-12,2.6215113639e-12,- 2.2655462787e-14,- 2.3898379798e-14,- 1.5822346004e-12,1.4014592713e-14,2.6215113639e-12,6.7713072988e-12,- 1.8031837762e-14,- 1.6405199435e-14,1.4140513359e-12,- 4.7096511176e-17,- 2.2655462787e-14,- 1.8031837762e-14,2.5193704546e-12,- 2.7497854222e-13,7.2254489606e-13,- 5.1804415247e-17,- 2.3898379798e-14,- 1.6405199435e-14,- 2.7497854222e-13,2.5859832833e-12,- 2.5681926955e-12,2.640564528e-17,- 1.5822346004e-12,1.4140513359e-12,7.2254489606e-13,- 2.5681926955e-12,4.7583323083e-12])
cartcovDtruenrev,tmeq2cti = sc.coveq2ct(eqcovDtruen,eqstate,anom,fr)
sc.printcov(cartcovDtruenrev,'ct','m',anom)
smu.printdiff(' cartcov - cartcovDtruenrev \n',cartcov,cartcovDtruenrev)
sc.printcov(tmcl2eq @ tmeq2cl,'tm','m',anom)
#smu.printdiff( ' tmct2eq - np.linalg.inv(tmeq2cti) \n', tmct2eq, np.linalg.inv(tmeq2cti));
#smu.printdiff( ' tmeq2cti - np.linalg.inv(tmct2eq) \n', tmeq2cti, np.linalg.inv(tmct2eq));

# compare combined tm
print('===============================================================================================\n' % ())
# fwd
print('Combined TM')
sc.printcov(tmct2eqtruen,'tm','m',anom)
print(' \n' % ())
print('Direct function conversion')
sc.printcov(tmct2eq,'tm','m',anom)
smu.printdiff(' tmct2eqtruen - tmct2eq \n',tmct2eqtruen,tmct2eq)
print('===============================================================================================\n' % ())
# bacwd
print('Combined TM')
sc.printcov(tmeq2cttruen,'tm','m',anom)
print(' \n' % ())
print('Direct function conversion')
sc.printcov(tmeq2cti,'tm','m',anom)
smu.printdiff(' tmeq2cttruen - tmeq2cti \n',tmeq2cttruen,tmeq2cti)

# ===============================================================================================
print('===============================================================================================\n' % ())
print('\n5.  Equinoctial Covariance from Cartesian #1 above (meana) \n' % ())
anom = 'meana'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
eqcovDmeana,tmct2eq = sc.covct2eq(cartcov,cartstate,anom,fr)
sc.printcov(eqcovDmeana,'eq','m',anom)
smu.printdiff(' eqcovmeana - eqcovDmeana \n',eqcovmeana,eqcovDmeana)
print('\n6.  Cartesian Covariance from Equinoctial #5 above \n' % ())
# try with intermediate one
cartcovDmeanarev,tmeq2ct = sc.coveq2ct(eqcovDmeana,eqstate,anom,fr) #jmb changed this from tmeq2cti
sc.printcov(cartcovDmeanarev,'ct','m',anom)
smu.printdiff(' cartcov - cartcovDmeanarev \n',cartcov,cartcovDmeanarev)
#Problem: DOUBLE CHECK THIS (Close but not the same as Matlab) -mjc
sc.printcov(tmct2eq @ tmeq2ct,'tm','m',anom)
#smu.printdiff( ' tmct2eq - np.linalg.inv(tmeq2cti) \n', tmct2eq, np.linalg.inv(tmeq2cti));
#smu.printdiff( ' tmeq2cti - np.linalg.inv(tmct2eq) \n', tmeq2cti, np.linalg.inv(tmct2eq));


# ===============================================================================================
print('===============================================================================================\n' % ())
print('\n7.  Flight Covariance from Cartesian #1 above \n' % ())
anom = 'meana'
anomflt = 'latlon'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
flcov,tmct2fl = sc.covct2fl(cartcov,cartstate,anomflt,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps)
if str(anomflt) == str('latlon'):
    sc.printcov(flcov,'fl','m',anomflt)
else:
    sc.printcov(flcov,'sp','m',anomflt)
print('\n8. Cartesian Covariance from Flight #7 above \n' % ())
cartcovfltrev,tmfl2ct = sc.covfl2ct(flcov,flstate,anomflt,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps)
sc.printcov(cartcovfltrev,'ct','m',anom)
smu.printdiff(' cartcov - cartcovfltrev \n',cartcov,cartcovfltrev)
print('\n9. Classical Covariance from Flight #7 above \n' % ())
print('\n-------- tm fl2cl --------- \n' % ())
print('\n-------- SKIPPING --------- \n' % ())
'''
classcovfl,tmfl2cl = sc.covfl2cltest(flcov,flstate,anomflt,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps)
sc.printcov(classcovfl,'cl','m',anom)
smu.printdiff(' classcov - classcovfl \n',classcovmeana,classcovfl)
sc.printcov(tmct2fl @ tmfl2ct,'tm','m',anom)
'''


# ===============================================================================================
print('===============================================================================================\n' % ())
print('\n7.  Flight Covariance from Cartesian #1 above \n' % ())
anom = 'meana'
anomflt = 'radec'
cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'y',anom,anomflt,ddpsi,ddeps)
flcov,tmct2fl = sc.covct2fl(cartcov,cartstate,anomflt,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps)
if str(anomflt) == str('latlon') == 1:
    sc.printcov(flcov,'fl','m',anomflt)
else:
    sc.printcov(flcov,'sp','m',anomflt)

print('\n8. Cartesian Covariance from Flight #7 above \n' % ())
cartcovfltrev,tmfl2ct = sc.covfl2ct(flcov,flstate,anomflt,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps)
sc.printcov(cartcovfltrev,'ct','m',anom)
smu.printdiff(' cartcov - cartcovfltrev \n',cartcov,cartcovfltrev)


print('\n9. Classical Covariance from Flight #7 above \n' % ())
print('\n-------- tm fl2cl --------- \n' % ())
print('\n-------- SKIPPING --------- \n' % ())
'''
classcovfl,tmfl2cl = sc.covfl2cltest(flcov,flstate,anomflt,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps)
sc.printcov(classcovfl,'cl','m',anom)
smu.printdiff(' classcov - classcovfl \n',classcovmeana,classcovfl)
# temp testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#         fprintf(1,'2.temp testing!!  find original transformations \n');
#         [classco, tmct2clO]  = covct2clO( cartcov,cartstate,strcat(anomeq1,'a') );
#         [cartcov1, tmcl2ctO] = covcl2ctO( classco,classstate,strcat(anomeq1,'a') );

#         fprintf(1,'-------- tm cl2ct old ---------\n');
#         printcov( tmcl2ctO,'tm','m',strcat(anomeq1,'a') );
#         fprintf(1,'-------- tm ct2cl old ---------\n');
#         printcov( tmct2clO,'tm','m',strcat(anomeq1,'a') );

#         fprintf(1,'Check consistency of both approaches tmct2cl-np.linalg.inv(tmcl2ct) diff pct over 1e-18 \n');
#         fprintf(1,'-------- accuracy of tm comparing ct2cl and cl2ct --------- \n');
#         tm1 = tmct2clO;
#         tm2 = np.linalg.inv(tmcl2ctO);
#         for i=1:6
#             for j = 1:6
#                 if (abs( tm1[i,j] - tm2[i,j] ) < small) || (abs(tm1[i,j]) < small)
#                     diffm[i,j] = 0.0;
#                   else
#                     diffm[i,j] = 100.0*( (tm1[i,j]-tm2[i,j]) / tm1[i,j]);
#                 end;
#             end;
#         end;
#         fprintf(1,'pct differences if over #4e \n',small);
#         fprintf(1,'#14.4f#14.4f#14.4f#14.4f#14.4f#14.4f\n',diffm' );
# temp testing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''


diffmm = np.zeros((6,6))
print('-------- accuracy of tm comparing cl2eq and eq2cl --------- \n' % ())
tm1 = tmcl2eq
tm2 = np.linalg.inv(tmeq2cl)
for i in range(6):
    for j in range(6):
        if (np.abs(tm1[i,j] - tm2[i,j]) < small) or (np.abs(tm1[i,j]) < small):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((tm1[i,j] - tm2[i,j]) / tm1[i,j])

print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))
print('-------- accuracy of tm comparing ct2eq and eq2ct --------- \n' % ())
tm1 = tmct2eq
tm2 = np.linalg.inv(tmeq2ct)
for i in range(6):
    for j in range(6):
        if (np.abs(tm1[i,j] - tm2[i,j]) < small) or (np.abs(tm1[i,j]) < small):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((tm1[i,j] - tm2[i,j]) / tm1[i,j])

print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))
print('-------- accuracy of tm comparing ct2eq and eq2ctintermediate --------- \n' % ())
tm1 = tmct2eq
tm2 = np.linalg.inv(tmeq2cti)
for i in range(6):
    for j in range(6):
        if (np.abs(tm1[i,j] - tm2[i,j]) < small) or (np.abs(tm1[i,j]) < small):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((tm1[i,j] - tm2[i,j]) / tm1[i,j])

print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))
print('-------- accuracy of tm comparing cl2fl and fl2cl --------- \n' % ())
tm1 = tmct2fl
tm2 = np.linalg.inv(tmfl2ct)
for i in range(6):
    for j in range(6):
        if (np.abs(tm1[i,j] - tm2[i,j]) < small) or (np.abs(tm1[i,j]) < small):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((tm1[i,j] - tm2[i,j]) / tm1[i,j])

print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))
temp = tmcl2eq @ tmct2cl
temp1 = tmcl2ct @ tmeq2cl
print('\n-------- tm combined ct2cl, cl2eq --------- \n' % ())
sc.printcov(temp,'tm','m',anomeq1+anomeq2)
print('\n-------- tm ct2eq --------- \n' % ())
sc.printcov(tmct2eq,'tm','m',anomeq1+anomeq2)
print('\n-------- tm combined eq2cl, cl2ct --------- \n' % ())
sc.printcov(temp1,'tm','m',anomeq1+anomeq2)
print('\n-------- tm eq2ct --------- \n' % ())
sc.printcov(tmeq2ct,'tm','m',anomeq1+anomeq2)
print('-------- accuracy of test tm ct2eq --------- \n' % ())
tm1 = temp
tm2 = tmct2eq
for i in range(6):
    for j in range(6):
        if (np.abs(tm1[i,j] - tm2[i,j]) < small) or (np.abs(tm1[i,j]) < small):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((tm1[i,j] - tm2[i,j]) / tm1[i,j])

print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))
print('-------- accuracy of test tm eq2ct --------- \n' % ())
tm1 = temp1
tm2 = tmeq2ct
for i in range(6):
    for j in range(6):
        if (np.abs(tm1[i,j] - tm2[i,j]) < small) or (np.abs(tm1[i,j]) < small):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((tm1[i,j] - tm2[i,j]) / tm1[i,j])

print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))

#anom = 'true';
classcovrev,tm = sc.covct2cl(cartcov,cartstate,anomeq1+anomeq2)
sc.printcov(classcovrev,'cl','m',anomeq1+anomeq2)
#        [cartcov1, tm]   = covcl2ct( classco1,classstate,strcat(anomeq1,anomeq2) );
#        printcov( cartcov1,'ct','m',strcat(anomeq1,anomeq2) );


print('===================== do the eq-ct-eq tests \n' % ())
eqcov = np.array([[1e-14,1e-16,1e-16,1e-16,1e-16,1e-16],[1e-16,1e-14,1e-16,1e-16,1e-16,1e-16],[1e-16,1e-16,1e-14,1e-16,1e-16,1e-16],[1e-16,1e-16,1e-16,1e-19,1e-16,1e-16],[1e-16,1e-16,1e-16,1e-16,1e-14,1e-16],[1e-16,1e-16,1e-16,1e-16,1e-16,1e-14]])
sc.printcov(eqcov,'eq','m',anom)
classco,tm = sc.coveq2cl(eqcov,eqstate,anom,fr)
print("classco")
print(classco)
cartco,tm = sc.covcl2ct(classco,classstate,anom)
print("cartco")
print(cartco)
print("-----")
# First element (0,0) is off in classco - mjc
classco,tm = sc.covct2cl(cartco,cartstate,anom)
print("classco")
print(classco)
eqcovmeannrev,tm = sc.covcl2eq(classco,classstate,anom,fr)
sc.printcov(eqcovmeannrev,'eq','m',anom)
for i in range(6):
    for j in range(6):
        if (np.abs(eqcov[i,j] - eqcovmeannrev[i,j]) < small) or (np.abs(eqcov[i,j]) < small):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((eqcov[i,j] - eqcovmeannrev[i,j]) / eqcov[i,j])

print('-------- accuracy --------- \n' % ())
print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))
# ---------------- reset original cartcov ----------------------------
cartcov = np.array([[1.0,0.0,0.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0,0.0,0.0],[0.0,0.0,1.0,0.0,0.0,0.0],[0.0,0.0,0.0,1e-06,0.0,0.0],[0.0,0.0,0.0,0.0,1e-06,0.0],[0.0,0.0,0.0,0.0,0.0,1e-06]])
cartcov = np.array([[0.81,1e-08,1e-08,1e-08,1e-08,1e-08],[1e-08,0.81,1e-08,1e-08,1e-08,1e-08],[1e-08,1e-08,0.81,1e-08,1e-08,1e-08],[1e-08,1e-08,1e-08,1e-06,1e-08,1e-08],[1e-08,1e-08,1e-08,1e-08,1e-06,1e-08],[1e-08,1e-08,1e-08,1e-08,1e-08,1e-06]])
cartcov = np.array([[100.0,0.01,0.01,0.0001,0.0001,0.0001],[0.01,100.0,0.01,0.0001,0.0001,0.0001],[0.01,0.01,100.0,0.0001,0.0001,0.0001],[0.0001,0.0001,0.0001,0.0001,1e-06,1e-06],[0.0001,0.0001,0.0001,1e-06,0.0001,1e-06],[0.0001,0.0001,0.0001,1e-06,1e-06,0.0001]])
# -------------------------------------------------------------
# ------------- start covariance conversion tests -------------
print('anomaly is %s \n' % (anom))
print('initial cartesian covariance \n' % ())
sc.printcov(cartcov,'ct','m',anom)
print(' -------- classical covariance conversions --------- \n' % ())
print(' -----  cartesian to classical covariance  --------- \n' % ())
classcov,tm = sc.covct2cl(cartcov,cartstate,anom)
print('tm check:')
print(tm)
sc.printcov(classcov,'cl','m',anom)
print('tm ct2cl \n' % ())
print((np.transpose(tm)))
print('tm ct2cl inv \n' % ())
print((np.transpose(np.linalg.inv(tm))))
tmold = np.linalg.inv(tm)
#        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
#        fprintf(1,'#16.12f#16.12f#16.12f#16.12f#16.12f#16.12f\n',(tm * tm')' );

#        fprintf(1,'Check the xxxx(transpose*mat) of classcov \n' );
#        fprintf(1,'#16.12f#16.12f#16.12f#16.12f#16.12f#16.12f\n',orth(classcov)' );
#        fprintf(1,'Check the xxxx(mat*transpose) of cartcov \n' );
#        fprintf(1,'#16.12f#16.12f#16.12f#16.12f#16.12f#16.12f\n',orth(cartcov)' );

#        # --- now check out what the scaling should be beteween
#        for i=1:6
#          for j = 1:6
#            diffm[i,j] = classcovtrace[i,j]/classcov[i,j];
#           end;
#         end;
#        fprintf(1,'scaling factors - should be \n');
#        fprintf(1,'#16e#16e#16e#16e#16e#16e\n',diffm' );
#        for i=1:6
#           for j = 1:6
#                if (abs( classcovt[i,j]-classcov[i,j] ) < small) | (abs(classcovt[i,j]) < small)
#                    diffm[i,j] = 0.0;
#                  else
#                    diffm[i,j] = 100.0*( (classcovt[i,j]-classcov[i,j]) / classcovt[i,j]);
#                  end;
#              end;
#          end;

#        # diff just to trace output
#        fprintf(1,'-------- accuracy to trace --------- \n');
#        fprintf(1,'pct differences if over #4e \n', small);
#        fprintf(1,'#14.4f#14.4f#14.4f#14.4f#14.4f#14.4f\n',diffm' );

print('----- classical to cartesian covariance  ---------- \n' % ())
cartco,tm = sc.covcl2ct(classcov,classstate,anom)
sc.printcov(cartco,'ct','m',anom)
print('tm cl2ct \n' % ())
print((np.transpose(tm)))
#        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
#        fprintf(1,'#16e#16e#16e#16e#16e#16e\n',(tm*tm')' );

for i in range(6):
    for j in range(6):
        if (np.abs(cartcov[i,j] - cartco[i,j]) < small) or (np.abs(cartcov[i,j]) < small):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((cartcov[i,j] - cartco[i,j]) / cartcov[i,j])

print('-------- accuracy --------- \n' % ())
print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))
for i in range(6):
    for j in range(6):
        if (np.abs(tmold[i,j] - tm[i,j]) < small) or (np.abs(tmold[i,j]) < small):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((tmold[i,j] - tm[i,j]) / tmold[i,j])

print('-------- tm accuracy --------- \n' % ())
print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))

print(' -------- flight covariance conversions --------- \n' % ())
print(' -----  cartesian to flight covariance  --------- \n' % ())
flcov,tm = sc.covct2fl(cartcov,cartstate,anom,ttt,jdut1,lod,xp,yp,terms,ddpsi,ddeps)
sc.printcov(flcov,'fl','m',anom)
print('tm ct2fl \n' % ())
print((np.transpose(tm)))
print('tm ct2fl inv \n' % ())
print((np.transpose(np.linalg.inv(tm))))
tmold = np.linalg.inv(tm)
#        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
#        fprintf(1,'#16.12f#16.12f#16.12f#16.12f#16.12f#16.12f\n',(tm * tm')' );

#        fprintf(1,'Check the xxxx(transpose*mat) of flcov \n' );
#        fprintf(1,'#16.12f#16.12f#16.12f#16.12f#16.12f#16.12f\n',(flcov' * flcov)' );

print('----- flight to cartesian covariance \n' % ())
cartco,tm = sc.covfl2ct(flcov,flstate,anomflt,ttt,jdut1,lod,xp,yp,2,ddpsi,ddeps)
sc.printcov(cartco,'ct','m',anom)
print('tm  fl2ct \n' % ())
print((np.transpose(tm)))
#        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
#        fprintf(1,'#16e#16e#16e#16e#16e#16e\n',(tm*tm')' );
for i in range(6):
    for j in range(6):
        if ((np.abs(cartcov[i,j] - cartco[i,j]) < small) or (np.abs(cartcov[i,j]) < small)):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((cartcov[i,j] - cartco[i,j]) / cartcov[i,j])

print('-------- accuracy --------- \n' % ())
print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))
for i in range(6):
    for j in range(6):
        if ((np.abs(tmold[i,j] - tm[i,j]) < small) or (np.abs(tmold[i,j]) < small)):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((tmold[i,j] - tm[i,j]) / tmold[i,j])

print('-------- tm accuracy --------- \n' % ())
print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))

print(' -------- equinoctial covariance conversions --------- \n' % ())
print(' -----  classical to equinoctial covariance ---------- \n' % ())
eqcov,tm = sc.covcl2eq(classcov,classstate,anom,fr)
sc.printcov(eqcov,'eq','m',anom)
print('tm cl2eq \n' % ())
print((np.transpose(tm)))
print('tm cl2eq inv \n' % ())
print((np.transpose(np.linalg.inv(tm))))
tmold = np.linalg.inv(tm)
#        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
#        fprintf(1,'#16.12f#16.12f#16.12f#16.12f#16.12f#16.12f\n',(tm * tm')' );

#        fprintf(1,'Check the xxxx(transpose*mat) of eqcov \n' );
#        fprintf(1,'#16.12f#16.12f#16.12f#16.12f#16.12f#16.12f\n',(eqcov' * eqcov)' );

print('----- equinoctial to classical covariance \n' % ())
classco,tm = sc.coveq2cl(eqcov,eqstate,anom,fr)
sc.printcov(classco,'cl','m',anom)
print('tm eq2cl \n' % ())
print((np.transpose(tm)))
#        fprintf(1,'Check the orthogonality (mat*transpose) of tm \n' );
#        fprintf(1,'#16e#16e#16e#16e#16e#16e\n',(tm*tm')' );

for i in range(6):
    for j in range(6):
        if ((np.abs(classcov[i,j] - classco[i,j]) < small) or (np.abs(classcov[i,j]) < small)):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((classcov[i,j] - classco[i,j]) / classcov[i,j])

print('-------- accuracy --------- \n' % ())
print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))
for i in range(6):
    for j in range(6):
        if ((np.abs(tmold[i,j] - tm[i,j]) < small) or (np.abs(tmold[i,j]) < small)):
            diffmm[i,j] = 0.0
        else:
            diffmm[i,j] = 100.0 * ((tmold[i,j] - tm[i,j]) / tmold[i,j])

print('-------- accuracy tm --------- \n' % ())
print('pct differences if over %4e \n' % (small))
print((np.transpose(diffmm)))

# ------------------------ check nav tests -----------------------------------
small = 1e-18
doall = 'n'
#        fprintf(1,'year mon day hms  magr  magv  r sig  vsig  max diag  max mat \n');
# print('      ecc         incl        maxdiag        maxdiff        magr        r sig   \n' % ())
for testnumbb in np.arange(1,12+1).reshape(-1):
    #        for testnumb = 1:122
    if testnumbb == 1:
        testnumb = 3
        satnum = 107
    if testnumbb == 2:
        testnumb = 16
        satnum = 11
    if testnumbb == 3:
        testnumb = 28
        satnum = 16609
    if testnumbb == 4:
        testnumb = 42
        satnum = 20052
    if testnumbb == 5:
        testnumb = 52
        satnum = 21867
    if testnumbb == 6:
        testnumb = 64
        satnum = 23019
    if testnumbb == 7:
        testnumb = 76
        satnum = 24780
    if testnumbb == 8:
        testnumb = 87
        satnum = 25013
    if testnumbb == 9:
        testnumb = 97
        satnum = 25544
    if testnumbb == 10:
        testnumb = 101
        satnum = 25634
    if testnumbb == 11:
        testnumb = 110
        satnum = 26354
    if testnumbb == 12:
        testnumb = 113
        satnum = 26405
    ###testcove
    print('year %5i ' % (year))
    print(' mon %4i ' % (mon))
    print(' day %3i ' % (day))
    print('%3i:%2i:%8.6f\n ' % (hr,min,sec))
    print('dut1 %8.6f s' % (dut1))
    print(' dat %3i s' % (dat))
    print(' xp %8.6f "' % (xp))
    print(' yp %8.6f "' % (yp))
    print(' lod %8.6f s\n' % (lod))
    print(' ddpsi %8.6f " ddeps  %8.6f\n' % (ddpsi,ddeps))
    #        fprintf(1,'order #3i  eqeterms #31  opt #3s \n',order,eqeterms,opt );
    print('units are km and km/s and km/s2\n' % ())
    print('eci: ' , (reci))
    print(' v:\n' , (veci))
    #, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
    ut1,tut1,jdut1,jdut1frac,utc,tai,tt,ttt,jdtt,jdttfrac,tdb,ttdb,jdtdb,jdtdbfrac = stu.convtime(year,mon,day,hr,min,sec,timezone,dut1,dat)
    print('ut1 %8.6f tut1 %16.12f jdut1 %18.11f\n' % (ut1,tut1,jdut1))
    print('utc %8.6f\n' % (utc))
    print('tai %8.6f\n' % (tai))
    print('tt  %8.6f ttt  %16.12f jdtt  %18.11f\n' % (tt,ttt,jdtt))
    print('tdb %8.6f ttdb %16.12f jdtdb %18.11f\n' % (tdb,ttdb,jdtdb))
    p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper = sc.rv2coe(reci,veci)
    print('          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n' % ())
    print('coes %11.4f%11.4f%13.9f%13.7f%14.8f%14.8f%14.8f%14.8f%14.8f\n' % (p,a,ecc,incl * rad2deg,omega * rad2deg,argp * rad2deg,nu * rad2deg,m * rad2deg,arglat * rad2deg))
    if (truelon and lonper):
      print(truelon * rad2deg,lonper * rad2deg)
    cartstate,classstate,flstate,eqstate,fr = sc.setcov(reci,veci,ttt,jdut1,lod,xp,yp,terms,'n',anom,anomflt,ddpsi,ddeps)
    print('new case --------------\n' % ())
    for i in range(6):
        eqcov[i,3] = eqcov[i,3] / tusec
    for i in range(6):
        eqcov[3,i] = eqcov[3,i] / tusec
    # --------------------- write out input data --------------------------
    if doall != 'y':
        #                printcov( eqcov,'eq','m',anom );
        pass
    sc.printcov(eqcov,'eq','m',anom)
    # ------ do conversions
    classcov,tm = sc.coveq2cl(eqcov,eqstate,anom,fr)
    sc.printcov(classcov,'cl','m',anom)
    cartcov,tm = sc.covcl2ct(classcov,classstate,anom)  #jmb change
    sc.printcov(cartcov,'ct','m',anom)
    classco,tm = sc.covct2cl(cartcov,cartstate,anom)  #jmb change
    sc.printcov(classco,'cl','m',anom)
    eqcovmeannrev,tm = sc.covcl2eq(classco,classstate,anom,fr)
    if doall != 'y':
        #                printcov( eqco,'eq','m',anom );
        pass
    sc.printcov(eqcovmeannrev,'eq','m',anom)
    for i in range(6):
        for j in range(6):
            if ((np.abs(eqcov[i,j] - eqcovmeannrev[i,j]) < small) or (np.abs(eqcov[i,j]) < small)):
                diffmm[i,j] = 0.0
            else:
                diffmm[i,j] = 100.0 * ((eqcov[i,j] - eqcovmeannrev[i,j]) / eqcov[i,j])
    if doall != 'y':
        #                fprintf(1,'accuracy \n');
#                fprintf(1,'#11.4f#11.4f#11.4f#11.4f#11.4f#11.4f\n',diffm' );
        pass
    #            fprintf(1,'accuracy \n');
#            fprintf(1,'#11.4f#11.4f#11.4f#11.4f#11.4f#11.4f\n',diffm' );
    magr = np.sqrt(cartcov[0,0] + cartcov[1,1] + cartcov[2,2])
    magv = np.sqrt(cartcov[3,3] + cartcov[4,4] + cartcov[5,5])
###    magrs = np.sqrt(possigma[0] ** 2 + possigma[1] ** 2 + possigma[2] ** 2) * 1000.0
###    magvs = np.sqrt(velsigma[0] ** 2 + velsigma[1] ** 2 + velsigma[2] ** 2) * 1000.0
    if doall != 'y':
        print('position mag %11.5f velocity %11.5f \n' % (magr,magv))
###        print('position sig from msg %11.5f velocity %11.5f \n' % (magrs,magvs))
        print('-------------------------------------------------\n' % ())
    md = np.amax(np.diag(np.abs(diffmm)))
    mm = np.amax(np.amax(np.abs(diffmm)))
    if doall != 'y':
        #                fprintf(1,'#3i #5i #3i #3i #3i:#2i:#8.6f #11.5f#11.5f#11.5f#11.5f#11.5f#11.5f\n', ...
#                          testnumb,year,mon,day,hr,min,sec,magr,magv,magrs,magvs,md,mm );
#                if ( md > 0.1 ) | ( mm > 0.1 )
#                    fprintf(1,'input equinoctial covariance \n');
#                    printcov( eqcov,'eq','m',anom );
#                    fprintf(1,'calculated equinoctial covariance \n');
#                    printcov( eqco,'eq','m',anom );
#                  end
        pass
###    print('%11.6f   %11.6f   %14.6f %14.6f %14.6f %11.6f %3i \n' % (classstate(2),classstate(3) * rad2deg,md,mm,magr,magrs,testnumb))

#Consolidated verifcation
small = 1e-18
cartcov = np.array([100.0,0.01,0.01,0.0001,0.0001,0.0001,0.01,100.0,0.01,0.0001,0.0001,0.0001,0.01,0.01,100.0,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,1e-06,1e-06,0.0001,0.0001,0.0001,1e-06,0.0001,1e-06,0.0001,0.0001,0.0001,1e-06,1e-06,0.0001])
# ct2cl meana
classcov = np.array([729.98470527,3.2622542748e-05,1.9882701518e-07,- 1.5267354781e-07,0.065719274495,- 0.065719721795,3.2622542748e-05,4.7913185965e-12,1.6984420631e-14,- 1.3041848226e-14,1.9093724426e-09,- 1.9144294314e-09,1.9882701518e-07,1.6984420631e-14,1.821030447e-12,1.8574598946e-13,2.2175774749e-11,- 2.216051549e-11,- 1.5267354781e-07,- 1.3041848226e-14,1.8574598946e-13,2.0516270435e-12,- 1.6735967689e-11,1.7016422633e-11,0.065719274495,1.9093724426e-09,2.2175774749e-11,- 1.6735967689e-11,7.2061287418e-06,- 7.2039903813e-06,- 0.065719721795,- 1.9144294314e-09,- 2.216051549e-11,1.7016422633e-11,- 7.2039903813e-06,7.2018610163e-06])
tmct2clm = np.array([- 0.17683319064,- 1.7135436636,1.0196363131,- 370.38514352,- 874.40859633,- 1530.3032219,- 3.2616946703e-08,- 1.1660253726e-07,- 8.1284934173e-08,- 3.3159300134e-05,- 0.00015533153893,- 3.295424911e-05,7.316107033e-08,- 1.349752469e-08,- 9.9950269315e-09,- 0.00010964037725,2.0227611381e-05,1.4978710924e-05,1.2294251173e-07,- 2.2681729232e-08,- 1.6796005174e-08,6.627761365e-05,- 1.2227592114e-05,- 9.0546315192e-06,8.613477438e-07,- 7.8604704844e-05,0.00011232839062,- 0.041439495269,- 0.06424279752,- 0.21663961239,- 8.1498812514e-07,7.8672438414e-05,- 0.00011220665948,0.041471465416,0.064465640062,0.21650469223])
cartcov1 = np.array([100.00000044,0.0099479437413,0.010077060115,9.9969642671e-05,9.9952965602e-05,9.9850351635e-05,0.0099479482267,99.999494542,0.010306186086,9.9911319273e-05,9.9975417086e-05,9.9435118941e-05,0.01007705827,0.010306163211,100.0001513,9.9894711415e-05,9.9685411584e-05,9.9664977e-05,9.9969640504e-05,9.9911316301e-05,9.9894701574e-05,0.00010000005983,1.0001511979e-06,1.0002207062e-06,9.9952963204e-05,9.9975432157e-05,9.9685393138e-05,1.000151189e-06,0.0001000003212,1.0006353763e-06,9.9850357726e-05,9.9435120969e-05,9.9665007269e-05,1.0002207092e-06,1.0006353547e-06,0.00010000073117])
tmcl2ctm = np.array([- 0.088298081674,- 2306206.639,3435082.7677,5870229.5255,- 1409736.6,- 1411582.407,- 0.85562340421,- 3574412.879,- 633740.24247,- 605792.22533,- 3323831.848,- 3332476.496,0.50913478173,- 12053876.478,- 469289.80979,0.0,- 5830333.0764,- 5832169.9765,0.00011429153396,- 47.014983029,- 6371.951065,3702.3488924,672.11739544,673.93260795,0.00026982048565,4373.4397244,1175.564633,- 1568.25429,6524.0298553,6530.5214031,0.00047221306357,- 6250.1359462,870.51518276,0.0,- 3890.4773582,- 3885.9568039])
print('diff cartcov - cartcov1 \n' % ())
print((np.transpose(cartcov) - np.transpose(cartcov1)))
print('pctdiff cartcov - cartcov1 \n' % ())
print((100.0 * ((np.transpose(cartcov) - np.transpose(cartcov1)) / np.transpose(cartcov))))
print('diff tmct2cl - np.linalg.inv(tmcl2ct) \n' % ())
tmct2clm = np.reshape(tmct2clm,(6,6))
tmcl2ctm = np.reshape(tmcl2ctm,(6,6))
print((np.transpose(tmct2clm) - np.transpose(np.linalg.inv(tmcl2ctm))))
print('pctdiff tmct2cl - np.linalg.inv(tmcl2ct) \n' % ())
print((100.0 * ((np.transpose(tmct2clm) - np.transpose(np.linalg.inv(tmcl2ctm))) / np.transpose(tmct2clm))))
print('diff tmcl2ct - np.linalg.inv(tmct2cl) \n' % ())
print((np.transpose(tmcl2ctm) - np.transpose(np.linalg.inv(tmct2clm))))
print('pctdiff tmcl2ct - np.linalg.inv(tmct2cl) \n' % ())
print((100.0 * ((np.transpose(tmcl2ctm) - np.transpose(np.linalg.inv(tmct2clm))) / np.transpose(tmcl2ctm))))
# ct2cl truea
classcov = np.array([729.98470527,3.2622542748e-05,1.9882701518e-07,- 1.5267354781e-07,0.065719274495,- 0.065719230531,3.2622542748e-05,4.7913185965e-12,1.6984420631e-14,- 1.3041848226e-14,1.9093724426e-09,- 1.9074392492e-09,1.9882701518e-07,1.6984420631e-14,1.821030447e-12,1.8574598946e-13,2.2175774749e-11,- 2.2149481645e-11,- 1.5267354781e-07,- 1.3041848226e-14,1.8574598946e-13,2.0516270435e-12,- 1.6735967689e-11,1.700795006e-11,0.065719274495,1.9093724426e-09,2.2175774749e-11,- 1.6735967689e-11,7.2061287418e-06,- 7.2069633456e-06,- 0.065719230531,- 1.9074392492e-09,- 2.2149481645e-11,1.700795006e-11,- 7.2069633456e-06,7.2078001121e-06])
tmct2clt = np.array([- 0.17683319064,- 1.7135436636,1.0196363131,- 370.38514352,- 874.40859633,- 1530.3032219,- 3.2616946703e-08,- 1.1660253726e-07,- 8.1284934173e-08,- 3.3159300134e-05,- 0.00015533153893,- 3.295424911e-05,7.316107033e-08,- 1.349752469e-08,- 9.9950269315e-09,- 0.00010964037725,2.0227611381e-05,1.4978710924e-05,1.2294251173e-07,- 2.2681729232e-08,- 1.6796005174e-08,6.627761365e-05,- 1.2227592114e-05,- 9.0546315192e-06,8.613477438e-07,- 7.8604704844e-05,0.00011232839062,- 0.041439495269,- 0.06424279752,- 0.21663961239,- 8.7495409985e-07,7.8531006978e-05,- 0.00011245460285,0.041448320342,0.064241169378,0.21663840674])
cartcov1 = np.array([99.999999878,0.0099838874085,0.010023672447,9.9990547615e-05,9.9985265964e-05,9.9951263568e-05,0.0099838796334,99.999667381,0.010334126229,9.994238021e-05,9.9878108156e-05,9.9776738608e-05,0.01002365448,0.010334115948,99.999719738,0.00010000501669,0.00010005691598,9.9940448393e-05,9.9990549496e-05,9.9942380941e-05,0.00010000501588,0.00010000004013,1.0000408179e-06,1.0002290301e-06,9.9985272285e-05,9.9878087092e-05,0.00010005685897,1.0000408235e-06,0.00010000003267,1.0002512404e-06,9.9951259431e-05,9.9776750265e-05,9.9940480844e-05,1.0002290269e-06,1.0002512422e-06,0.00010000129048])
tmcl2ctt = np.array([- 0.088298080209,255151.93503,3435082.8085,5870229.5018,- 1409736.6043,- 1410321.5189,- 0.85562340076,2472465.6056,- 633740.25,- 605792.21528,- 3323831.8898,- 3329499.8116,0.50913478778,- 1471229.3403,- 469289.81536,0.0,- 5830333.0515,- 5826960.3806,0.00011429153431,- 1269.885894,- 6371.9510378,3702.3489388,672.11738428,673.33060821,0.00026982048904,- 7476.387043,1175.5646279,- 1568.2542948,6524.0298289,6524.6880053,0.00047221306155,801.04911367,870.51517904,0.0,- 3890.4774044,- 3882.4857291])
print('diff cartcov - cartcov1 \n' % ())
print((np.transpose(cartcov) - np.transpose(cartcov1)))
print('pctdiff cartcov - cartcov1 \n' % ())
print((100.0 * ((np.transpose(cartcov) - np.transpose(cartcov1)) / np.transpose(cartcov))))
print('diff tmct2cl - np.linalg.inv(tmcl2ct) \n' % ())
tmct2clt = np.reshape(tmcl2ctt,(6,6))
tmcl2ctt = np.reshape(tmcl2ctt,(6,6))
print((np.transpose(tmct2clt) - np.transpose(np.linalg.inv(tmcl2ctt))))
print('pctdiff tmct2cl - np.linalg.inv(tmcl2ct) \n' % ())
print((100.0 * ((np.transpose(tmct2clt) - np.transpose(np.linalg.inv(tmcl2ctt))) / np.transpose(tmct2clt))))
print('diff tmcl2ct - np.linalg.inv(tmct2cl) \n' % ())
print((np.transpose(tmcl2ctt) - np.transpose(np.linalg.inv(tmct2clt))))
print('pctdiff tmcl2ct - np.linalg.inv(tmct2cl) \n' % ())
print((100.0 * ((np.transpose(tmcl2ctt) - np.transpose(np.linalg.inv(tmct2clt))) / np.transpose(tmcl2ctt))))
# lance test

eqCovmeana = np.array([729.98470594,- 5.1232022015e-05,- 5.769702175e-05,1.938917477e-07,2.1327374042e-07,- 5.9997332907e-07,- 5.1232022015e-05,6.1778006946e-12,2.62153547e-12,- 2.2655470832e-14,- 2.3898387998e-14,4.2135772123e-12,- 5.769702175e-05,2.62153547e-12,6.7713618731e-12,- 1.8031914124e-14,- 1.6405273668e-14,- 3.6135835764e-12,1.938917477e-07,- 2.2655470832e-14,- 1.8031914124e-14,2.5193665805e-12,- 2.7497822431e-13,7.117843681e-13,2.1327374042e-07,- 2.3898387998e-14,- 1.6405273668e-14,- 2.7497822431e-13,2.5859797681e-12,- 2.5800265285e-12,- 5.9997332907e-07,4.2135772123e-12,- 3.6135835764e-12,7.117843681e-13,- 2.5800265285e-12,1.1608008196e-11])
eqCovmeana = np.reshape(eqCovmeana,(6,6))
#Debugging
# xxx = (np.linalg.svd(eqCovmeana)[1])
# print("svd:")
# print(xxx)
# print("diag:")
# print(np.diag(xxx))
s1eqA0 = np.sqrt(np.diag(np.linalg.svd(eqCovmeana)[1]))
print(s1eqA0)
eqCovmeann = np.array([0.043069707203,3.9352341277e-07,4.431823694e-07,- 1.4893213124e-09,- 1.6381982769e-09,4.608515197e-09,3.9352341277e-07,6.1778006946e-12,2.62153547e-12,- 2.2655470832e-14,- 2.3898387998e-14,4.2135772123e-12,4.431823694e-07,2.62153547e-12,6.7713618731e-12,- 1.8031914124e-14,- 1.6405273668e-14,- 3.6135835764e-12,- 1.4893213124e-09,- 2.2655470832e-14,- 1.8031914124e-14,2.5193665805e-12,- 2.7497822431e-13,7.117843681e-13,- 1.6381982769e-09,- 2.3898387998e-14,- 1.6405273668e-14,- 2.7497822431e-13,2.5859797681e-12,- 2.5800265285e-12,4.6085151971e-09,4.2135772123e-12,- 3.6135835764e-12,7.117843681e-13,- 2.5800265285e-12,1.1608008196e-11])
eqCovmeann = np.reshape(eqCovmeann,(6,6))
#np.linalg.svd(eqCovmeann)
s1eqN0 = np.sqrt(np.diag(np.linalg.svd(eqCovmeann)[1]))
a = np.array([s1eqA0[0,0],s1eqN0[0,0],s1eqA0[1,1],s1eqN0[1,1],s1eqA0[2,2],s1eqN0[2,2],s1eqA0[3,3],s1eqN0[3,3],s1eqA0[4,4],s1eqN0[4,4],s1eqA0[5,5],s1eqN0[5,5]])
# test sigmapts -------------------------------------------
cov3 = np.zeros([3,3])
cov3[0,0] = 12559.93762571587
cov3[0,1] = 12101.56371305036
cov3[0,2] = - 440.3145384949657
cov3[1,0] = 12101.56371305036
cov3[1,1] = 12017.77368889201
cov3[1,2] = 270.3798093532698
cov3[2,0] = - 440.3145384949657
cov3[2,1] = 270.3798093532698
cov3[2,2] = 4818.009967057008
r1 = np.zeros(3)
v1 = np.zeros(3)
r1[0] = - 28578179.1674123
r1[1] = 30995616.0912026
r1[2] = - 2261.22833628375
v1[0] = 33754.7320575159
v1[1] = - 2083.88260114718
v1[2] = 0.518744171831942
r1[0] = - 30762454.8061775
r1[1] = - 28804817.6521682
r1[2] = - 991451.166480117
print('cov3 starting \n' % ())
#    fprintf(1,'#20.10e #20.10e #20.10e #20.10e #20.10e #20.10e\n',cov3);
print(cov3)
#np.linalg.eig(cov3)
d, eigenaxes = np.linalg.eig(cov3)
eigenvalues = np.sqrt(d)
# print(eigenvalues) #Not in same order as matlab but same values
#eigenvalues
#     fprintf(1,'eigenaxes  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  m \n', eigenaxes(1:6) );
#     fprintf(1,'eigenaxes  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  m \n', eigenaxes(7:12) );
#     fprintf(1,'eigenaxes  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  m \n', eigenaxes(13:18) );
#     fprintf(1,'eigenaxes  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  m \n', eigenaxes(19:24) );
#     fprintf(1,'eigenaxes  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  m \n', eigenaxes(25:30) );
#     fprintf(1,'eigenaxes  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  m \n', eigenaxes(31:36) );
#     fprintf(1,'eigenvalues   #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  #16.8f  m \n', eigenvalues(1,1), eigenvalues(2,2), eigenvalues(3,3), ...
#                eigenvalues(4,4), eigenvalues(5,5), eigenvalues(6,6) );
#     fprintf(1,'mag eigenvalues  #16.8f m  \n', sqrt(eigenvalues(1,1)^2+ eigenvalues(2,2)^2+ eigenvalues(3,3)^2)  );
#     fprintf(1,'mag eigenvalues  #16.8f m/s \n', sqrt(eigenvalues(4,4)^2+ eigenvalues(5,5)^2+ eigenvalues(6,6)^2)  );

# ----- setup 12 initial states

# ---- find sigma points
# sigmapts = poscov2pts(r1,cov3)
# sigmapts
# yu,covout = remakecov(sigmapts)
# covout

# test pos/vel
r1[0] = - 30762454.8061775
r1[1] = - 28804817.6521682
r1[2] = - 991451.166480117
v1[0] = 2102.69117282211
v1[1] = - 2239.85220168736
v1[2] = - 140.227078892292
cov2 = np.zeros([6,6])
cov2[0,0] = 12559.93762571587
cov2[0,1] = 12101.56371305036
cov2[0,2] = - 440.3145384949657
cov2[0,3] = - 0.8507401236198346
cov2[0,4] = 0.9383675791981778
cov2[0,5] = - 0.0318596430999798
cov2[1,1] = 12017.77368889201
cov2[1,2] = 270.3798093532698
cov2[1,3] = - 0.8239662300032132
cov2[1,4] = 0.9321640899868708
cov2[1,5] = - 0.001327326827629336
cov2[2,2] = 4818.009967057008
cov2[2,3] = 0.02033418761460195
cov2[2,4] = 0.03077663516695039
cov2[2,5] = 0.1977541628188323
cov2[3,3] = 5.774758755889862e-05
cov2[3,4] = - 6.396031584925255e-05
cov2[3,5] = 1.079960679599204e-06
cov2[4,4] = 7.24599391355188e-05
cov2[4,5] = 1.03146660433274e-06
cov2[5,5] = 1.870413627417302e-05
cov2[1,0] = 12101.56371305036
cov2[2,0] = - 440.3145384949657
cov2[3,0] = - 0.8507401236198346
cov2[4,0] = 0.9383675791981778
cov2[5,0] = - 0.0318596430999798
cov2[2,1] = 270.3798093532698
cov2[3,1] = - 0.8239662300032132
cov2[4,1] = 0.9321640899868708
cov2[5,1] = - 0.001327326827629336
cov2[3,2] = 0.02033418761460195
cov2[4,2] = 0.03077663516695039
cov2[5,2] = 0.1977541628188323
cov2[4,3] = - 6.396031584925255e-05
cov2[5,3] = 1.079960679599204e-06
cov2[5,4] = 1.03146660433274e-06


# sigmapts = posvelcov2pts(r1,v1,cov2)
# yu,covout = remakecov(sigmapts)
# print(covout)
