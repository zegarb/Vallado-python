#     -----------------------------------------------------------------
#
#                              Ex7_2.m
#
#  this file demonstrates example 7-2.
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
import orbit_utils as obu
from space_constants import *
import space_conversions as sc
import spacemath_utils as smu
import spacetime_utils as stu
import os

#    re = 6378.137;
#    mu = 3.986004418e5;
#    tu = 86400.0
#    re = 6378.145;
#    mu = 3.986005e5;

#    re = 149597870.0;  # km in 1 au
#    mu = 1.32712428e11;
#    tu = 86400.0;

#    re = 1.0;  # 1 au
#    mu = 1.0;
#    tu = 1.0 / 58.132440906; # days in one solar tu

casenum = 0  #determines which input data is used. book uses subset of 0

#  casenum = 4;  # test

diffsites = 'n'
# typerun = 'l'; # laplace
# typerun = 'g'; # gauss
# typerun = 'd'; # doubler
# typerun = 'o'; # gooding
typerun = 'a' #run them all!

if casenum == 0:
    fn = os.path.join(os.path.dirname(__file__), "data",
                      "Sat11Access.dat")
else:
    fn = os.path.join(os.path.dirname(__file__), "data",
                      "Sat11Ex" + str(casenum) + ".dat")
filedat = np.loadtxt(fn)
print(filedat)
if casenum == 0:
    r2ans = np.array([5897.954130507, 5791.046114526, 6682.733686585])
    v2ans = np.array([-4.393910234, 4.576816355, 1.482423676])
    dut1 = - 0.609641
    dat = 35
    xp = 0.137495 * arcsec2rad
    yp = 0.342416 * arcsec2rad
    lod = 0.0
    timezone = 0
    terms = 0
    ddpsi = 0.0
    ddeps = 0.0

if casenum == 1:
    # at 8-20-07 11:50,
    r2ans = np.array([5897.954130507, 5791.046114526, 6682.733686585])
    v2ans = np.array([-4.393910234, 4.576816355, 1.482423676])
    dut1 = - 0.1639883
    dat = 33
    xp = 0.210428 * arcsec2rad
    yp = 0.286899 * arcsec2rad
    lod = 0.0
    timezone = 0
    terms = 0
    ddpsi = 0.0
    ddeps = 0.0

if casenum == 2:
    # at 8-20-12 11:48:28.000 center time,
#             year  =  2012;
#             mon   =   8;
#             day   =  20;
#             hr    =  11;
#             min   =  55;
#             sec   =  28.0000;
    dut1 = - 0.6096413
    dat = 35
    xp = 0.137495 * arcsec2rad
    yp = 0.342416 * arcsec2rad
    lod = 0.0
    timezone = 0
    terms = 0
    ddpsi = 0.0
    ddeps = 0.0

# set eop parameters for new cases
if casenum == 3:
    dut1 = 0
    dat = 34
    xp = 0.0
    yp = 0.0
    lod = 0.0
    timezone = 0
    terms = 0
    ddpsi = 0.0
    ddeps = 0.0
    obs1 = 1
    obs2 = 2
    obs3 = 3

if casenum == 4:
    dut1 = - 0.1069721
    dat = 37
    xp = 0.14847 * arcsec2rad
    yp = 0.246564 * arcsec2rad
    lod = 0.0
    timezone = 0
    terms = 0
    ddpsi = 0.0
    ddeps = 0.0
    obs1 = 1
    obs2 = 2
    obs3 = 3
    #2021 11 12 59530  0.148470  0.246564 -0.1069721  0.0001709 -0.113370 -0.007262  0.000222 -0.000052  37

match casenum:
    case 0:
        answ = '# 0  ans a = 12246.023  e = 0.2000  i = 40.00  330.000          km'
    case 1:
        answ = '# 1  old Book      2.0  ---- COEs: a= 12756.274 e= 0.20000 i= 62.3000   40.00   20.00 Topocentric case'
    case 2:
        answ = '# 2  ESC pg 282   1.5 ---- COEs:  a= 8306.247 e= 0.16419 i= 32.8780  136.53  203.95 Topocentric case -2519.36 8068.22 -2664.48 -5.9319 0.07621 2.60215'
        diffsites = 'y'
    case 3:
        answ = '# 3  ESC pg 288  2.0 ---- COEs:  a= 12756.274 e= 0.20000 i= 62.3000   40.00   20.00 Topocentric case'
        diffsites = 'n'
    case 4:
        answ = '# 4  ESC pg 289   10.0 ---- COEs:  a= 63781.37 e= 0.30000 i= 45.0000   45.00   45.00 Topocentric case'
    case 5:
        answ = '# 5  ESC pg 290    1.33 ---- COEs:  a= 8540.963 e= 0.17900 i= 32.8800  224.87  255.20 Topocentric case'
        diffsites = 'y'
    case 6:
        answ = '# 6  ESC pg 291    1.44 ---- COEs:  a= 9184.517 e= 0.23167 i= 33.2810  205.11  161.79 Topocentric case'
        diffsites = 'n'
    case 7:
        answ = '# 7  McCuskey pg 81 2.8                                          '
    case 8:
        answ = '# 8  Taff pg 232   2.67                                          '
    case 9:
        answ = '# 9  Taff pg 232   2.67 ---- COEs:  a= 26627.957 e= 0.65060 i= 64.2800  121.56  318.94 Topocentric case'
    case 10:
        answ = '# 10 Gauss Orig - Montenbruck pg 245 Heliocentric!               '
    case 11:
        answ = '# 11 Rockwell  Obj J91A01P                                       '
    case 12:
        answ = '# 12 test           2.0   ---- COEs:  a= 42163.95 e= 0.00100 i=  3.3000   40.00   20.00 Topocentric case'
    case 13:
        answ = '# 13 Der AAS 19-626 1:LEO 8597 a= 7358.3561 e= 0.001327600 i= 82.9703 263.5730 322.2972 79.6467 79.4971 41.9439 -14.5  -5504.8  4880.2  1.2167  4.8075  5.4409'
    case 14:
        answ = '# 14 Der AAS 19-626 2:MEO 28983 a= 8718.4092 e= 0.232518762 i= 105.6928 19.1724 269.5213 104.1749 77.4300 13.6962 8209  2261.3  1993.4  0.5708  -1.7319  6.4897 '
    case 15:
        answ = '# 15 Der AAS 19-626 3:SS 22824 a= 7183.4835 e= 0.263394025 i= 108.9880 307.9395 272.8097 105.1106 74.7419 17.9203 3632.38  -5828.17  2088.44  -2.12272  -1.21801  7.04160'
    case 16:
        answ = '# 16 Der AAS 19-626 4:HEO 20959 a= 26114.3997 e= 0.228037946 i= 50.5078 246.0637 207.2528 251.0742 277.0253 98.3270 16946.69  -3286.642  20413.21  0.6302  3.6057  -1.0762'
    case 17:
        answ = '# 17 Der AAS 19-626 5: HEO 36970 a= 28861.7539 e= 0.750899832 i= 26.4904 16.1794 55.8578 202.3578 270.6766 258.2157 1977.27  -37010.52  -17989.04  1.6269  1.6111  0.5452'
    case 18:
        answ = '# 18 Der AAS 19-626 6:Moly 9880 a= 27139.7858 e= 0.703343124 i= 63.6917 110.4936 250.7611 222.8043 306.6257 113.5655 -6816.554 -14641.890  23282.890  1.89019  -0.74702  -3.05218'
    case 19:
        answ = '# 19 Der AAS 19-626 7:GEO 27566 a= 42162.5679 e= 0.001093761 i= 3.4358 62.2774 204.0885 80.0126 79.8892 284.1011 40903.78  -9893.48  -2450.20  0.72866  2.98740  0.04471'
    case 20:
        answ = '# 20 Der AAS 19-626 8:GEO 27566 a= 42162.5679 e= 0.001093761 i= 3.4358 62.2774 204.0885 80.0126 79.8892 284.1011 40903.78  -9893.48  -2450.20  0.72866  2.98740  0.04471'
    case 21:
        answ = '# 21 Der AAS 19-626 9:GPS 377xx  a= 7792.4746 e= 0.001021058 i= 51.9889 320.6393 92.6220 310.1230 310.2124 42.7449 6485.933  -1110.318  4164.441  -1.707715  5.583869  4.138020'
    case 22:
        answ = '# 22 Der AAS 19-626 10:GPS 377xx  a= 7792.4746 e= 0.001021058 i= 51.9889 320.6393 92.6220 310.1230 310.2124 42.7449 6485.933  -1110.318  4164.441  -1.707715  5.583869  4.138020'
    case 23:
        answ = '# 23 Curtis pg 287  a= 10000 e= 0.1 i= 30.0000  270.0000 90.000 45.010'
    case _:
        print('Invalid option\n')

# ------ read all the data in and process
numobs = 3

obs1 = 0
obs2 = 1
obs3 = 2
yeararr = filedat[:, 2]
monarr = filedat[:, 1]
dayarr = filedat[:, 0]
hrarr = filedat[:, 3]
minarr = filedat[:, 4]
secarr = filedat[:, 5]
latarr = filedat[:, 6] * deg2rad
lonarr = filedat[:, 7] * deg2rad
altarr = filedat[:, 8]
rtascarr = filedat[:, 9] * deg2rad
declarr = filedat[:, 10] * deg2rad

#load data into x y z arrays
# for rec in filedat:
#   xxx = rec.replace("   ", ", ")
#   xxx = xxx.split(", ")
#   fd = [ddd.strip() for ddd in xxx]
#   print(fd)
#   yeararr.append(float(fd[2]))
#   monarr.append(float(fd[1]))
#   dayarr.append(float(fd[0]))
#   hrarr.append(float(fd[3]))
#   minarr.append(float(fd[4]))
#   secarr.append(float(fd[5]))
#   latarr.append(float(fd[6]) * deg2rad)
#   lonarr.append(float(fd[7]) * deg2rad)
#   altarr.append(float(fd[8]))
#   rtascarr.append(float(fd[9]) * deg2rad)
#   declarr.append(float(fd[10]) * deg2rad)

obsrecarr = []
for j in range(numobs):
    trec = {}
    trec['jd'], trec['jdf'] = \
        stu.jday(yeararr[j], monarr[j], dayarr[j], hrarr[j],
                 minarr[j], secarr[j])
    trec['latgd'] = latarr[j]
    trec['lon'] = lonarr[j]
    trec['alt'] = altarr[j]
    ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb,\
        jdtdb, jdtdbfrac = \
            stu.convtime(yeararr[j], monarr[j], dayarr[j],
                         hrarr[j], minarr[j], secarr[j], 0, dut1, dat)
    trec['rs'], trec['vs'] = obu.site(latarr[j], lonarr[j], altarr[j])
    trec['ttt'] = ttt
    trec['jdut1'] = jdut1
    trec['jdut1frac'] = jdut1frac
    trec['xp'] = xp
    trec['yp'] = yp
    trec['rtasc'] = rtascarr[j]
    trec['decl'] = declarr[j]
    obsrecarr.append(trec)

rtasc1 = obsrecarr[obs1]['rtasc']
rtasc2 = obsrecarr[obs2]['rtasc']
rtasc3 = obsrecarr[obs3]['rtasc']

decl1 = obsrecarr[obs1]['decl']
decl2 = obsrecarr[obs2]['decl']
decl3 = obsrecarr[obs3]['decl']

jd1 = obsrecarr[obs1]['jd']
jdf1 = obsrecarr[obs1]['jdf']
jd2 = obsrecarr[obs2]['jd']
jdf2 = obsrecarr[obs2]['jdf']
jd3 = obsrecarr[obs3]['jd']
jdf3 = obsrecarr[obs3]['jdf']

rs1 = obsrecarr[obs1]['rs']
vs1 = obsrecarr[obs1]['vs']
rs2 = obsrecarr[obs2]['rs']
vs2 = obsrecarr[obs2]['vs']
rs3 = obsrecarr[obs3]['rs']
vs3 = obsrecarr[obs3]['vs']

year, mon, day, hr, min, second = stu.invjday(obsrecarr[obs1]['jd'],
                                              obsrecarr[obs1]['jdf'])
utc = second
ut1 = utc + dut1
tai = utc + dat
tt = tai + 32.184
jdut1, jdut1frac = stu.jday(year, mon, day, hr, min, ut1)
jdtt, jdttfrac = stu.jday(year, mon, day, hr, min, tt)
ttt = (jdtt - 2451545.0) / 36525.0
print('year %5i ' % (year))
print('mon %4i ' % (mon))
print('day %3i ' % (day))
print('hr %3i:%2i:%8.6f\n' % (hr, min, second))
print('dut1 %8.6f s' % (dut1))
print(' dat %3i s' % (dat))
print(' xp %8.6f "' % (xp))
print(' yp %8.6f "' % (yp))
print(' lod %8.6f s\n' % (lod))
# -------------- convert each site vector from ecef to eci -----------------
a = np.array([[0], [0], [0]])

year, mon, day, hr, min, sec = stu.invjday(jd1, jdf1)
(ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb,
 ttdb, jdtdb, jdtdbfrac) = stu.convtime(year, mon, day, hr, min, sec,
                                        timezone, dut1, dat)
rsite1, vseci, aeci = sc.ecef2eci(rs1, vs1, a, ttt, jdut1 + jdut1frac, lod, xp,
                                   yp, 2, ddpsi, ddeps)
year, mon, day, hr, min, sec = stu.invjday(jd2, jdf2)
(ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb,
 jdtdb, jdtdbfrac) = stu.convtime(year, mon, day, hr, min, sec, timezone, dut1,
                                  dat)
rsite2, vseci, aeci = sc.ecef2eci(rs2, vs2, a, ttt, jdut1 + jdut1frac, lod, xp,
                                  yp, 2, ddpsi, ddeps)
year, mon, day, hr, min, sec = stu.invjday(jd3, jdf3)
(ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb,
 jdtdb, jdtdbfrac) = stu.convtime(year, mon, day, hr, min, sec, timezone, dut1,
                                  dat)
rsite3, vseci, aeci = sc.ecef2eci(rs3, vs3, a, ttt, jdut1 + jdut1frac, lod, xp,
                                  yp, 2, ddpsi, ddeps)

# ---------------------- run the angles-only routine ------------------
if typerun == 'l' or typerun == 'a':
    r2, v2 = obu.anglesl(decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1,
                         jdf1, jd2, jdf2, jd3, jdf3, diffsites, rsite1,
                         rsite2, rsite3)
    #                          rtasc3, jd1, jd2, jd3, rs1, rs2, rs3, re, mu, tu );
    processtype = 'anglesl'

# -------------- write out answer --------------
print('\n\ninputs: \n\n' % ())
latgc, latgd, lon, alt = sc.ecef2ll(rs1)

print('Site obs1 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n'
      % (rsite1[0], rsite1[1], rsite1[2], latgd * rad2deg, lon * rad2deg,
         alt * 1000))
latgc, latgd, lon, alt = sc.ecef2ll(rs2)
print('Site obs2 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n'
      % (rsite2[0], rsite2[1], rsite2[2], latgd * rad2deg, lon * rad2deg,
         alt * 1000))
latgc, latgd, lon, alt = sc.ecef2ll(rs3)
print('Site obs3 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n'
      % (rsite3[0], rsite3[1], rsite3[2], latgd * rad2deg, lon * rad2deg,
         alt * 1000))
year, mon, day, hr, min, sec = stu.invjday(jd1, jdf1)
print('obs#1 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
      % (year, mon, day, hr, min, sec, rtasc1 * rad2deg, decl1 * rad2deg))
year, mon, day, hr, min, sec = stu.invjday(jd2, jdf2)
print('obs#2 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
      % (year, mon, day, hr, min, sec, rtasc2 * rad2deg, decl2 * rad2deg))
year, mon, day, hr, min, sec = stu.invjday(jd3, jdf3)
print('Obs#3 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
      % (year, mon, day, hr, min, sec, rtasc3 * rad2deg, decl3 * rad2deg))
print('\nsolution by %s \n\n' % (processtype))
print('r2     %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n'
      % (r2[0]/re, r2[1]/re, r2[2]/re, r2[0], r2[1], r2[2]))
#        fprintf(1, 'r2 ans #11.7f   #11.7f  #11.7f er    #11.7f  #11.7f  #11.7f km \n', r2ans/re, r2ans);

print('v2     %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n'
      % (v2[0]/velkmps, v2[1]/velkmps, v2[2]/velkmps, v2[0], v2[1], v2[2]))
#        fprintf(1, 'v2 ans #11.7f   #11.7f  #11.7f er/tu #11.7f  #11.7f  #11.7f km/s\n', v2ans/velkmps, v2ans);

p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = \
    sc.rv2coeh(r2, v2, re, mu)
print('         p km          a km         ecc       incl deg     raan deg    '
      'argp deg     nu deg      m deg  \n')
print('coes p=%11.4f\n' % (p))
print('coes a=%11.4f\n' % (a))
print('coes ecc=%13.9f\n' % (ecc))
print('coes incl=%13.7f\n' % (incl * rad2deg))
print('coes omega=%11.5f\n' % (omega * rad2deg))
print('coes argp=%11.5f\n' % (argp * rad2deg))
if nu: print('coes nu=%11.5f\n' % (nu * rad2deg))
print('coes m=%11.5f\n' % (m * rad2deg))
print('%s \n' % (answ))
if typerun == 'g' or typerun == 'a':
    r2, v2 = obu.anglesg(decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1,
                         jdf1, jd2, jdf2, jd3, jdf3, rsite1, rsite2, rsite3)
    processtype = 'anglesg'

# -------------- write out answer --------------
print('\n\ninputs: \n\n' % ())
latgc, latgd, lon, alt = sc.ecef2ll(rs1)

print('Site obs1 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n'
      % (rsite1[0], rsite1[1], rsite1[2], latgd * rad2deg, lon * rad2deg,
          alt * 1000))
latgc, latgd, lon, alt = sc.ecef2ll(rs2)
print('Site obs2 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n'
      % (rsite2[0], rsite2[1], rsite2[2], latgd * rad2deg, lon * rad2deg,
          alt * 1000))
latgc, latgd, lon, alt = sc.ecef2ll(rs3)
print('Site obs3 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n'
      % (rsite3[0], rsite3[1], rsite3[2], latgd * rad2deg, lon * rad2deg,
         alt * 1000))
year, mon, day, hr, min, sec = stu.invjday(jd1, jdf1)
print('obs#1 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
      % (year, mon, day, hr, min, sec, rtasc1 * rad2deg, decl1 * rad2deg))
year, mon, day, hr, min, sec = stu.invjday(jd2, jdf2)
print('obs#2 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
      % (year, mon, day, hr, min, sec, rtasc2 * rad2deg, decl2 * rad2deg))
year, mon, day, hr, min, sec = stu.invjday(jd3, jdf3)
print('Obs#3 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
      % (year, mon, day, hr, min, sec, rtasc3 * rad2deg, decl3 * rad2deg))
print('\nsolution by %s \n\n' % (processtype))
print('r2     %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  %11.7f km \n'
      % (r2[0] / re, r2[1] / re, r2[2] / re, r2[0], r2[1], r2[2]))
#        fprintf(1, 'r2 ans #11.7f   #11.7f  #11.7f er    #11.7f  #11.7f  #11.7f km \n', r2ans/re, r2ans);

print('v2     %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  %11.7f km/s\n'
      % (v2[0] / velkmps, v2[1] / velkmps, v2[2] / velkmps, v2[0], v2[1],
         v2[2]))
#        fprintf(1, 'v2 ans #11.7f   #11.7f  #11.7f er/tu #11.7f  #11.7f  #11.7f km/s\n', v2ans/velkmps, v2ans);

p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = \
    sc.rv2coeh(r2, v2, re, mu)
print('         p km          a km         ecc       incl deg     raan deg    '
      'argp deg     nu deg      m deg  \n')
print('coes p=%11.4f\n' % (p))
print('coes a=%11.4f\n' % (a))
print('coes ecc=%13.9f\n' % (ecc))
print('coes incl=%13.7f\n' % (incl * rad2deg))
print('coes omega=%11.5f\n' % (omega * rad2deg))
print('coes argp=%11.5f\n' % (argp * rad2deg))
if nu: print('coes nu=%11.5f\n' % (nu * rad2deg))
print('coes m=%11.5f\n' % (m * rad2deg))
print('%s \n' % (answ))
if typerun == 'd' or typerun == 'a':
    r2, v2 = obu.anglesdr(decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1,
                          jdf1, jd2, jdf2, jd3, jdf3, rsite1, rsite2, rsite3,
                          re, mu)
    processtype = 'anglesdr'

# -------------- write out answer --------------
print('\n\ninputs: \n\n' % ())
latgc, latgd, lon, alt = sc.ecef2ll(rs1)

print('Site obs1 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n'
      % (rsite1[0], rsite1[1], rsite1[2], latgd * rad2deg, lon * rad2deg,
         alt * 1000))
latgc, latgd, lon, alt = sc.ecef2ll(rs2)
print('Site obs2 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n'
      % (rsite2[0], rsite2[1], rsite2[2], latgd * rad2deg, lon * rad2deg,
          alt * 1000))
latgc, latgd, lon, alt = sc.ecef2ll(rs3)
print('Site obs3 %11.7f %11.7f %11.7f km  lat %11.7f lon %11.7f alt %11.7f  \n'
      % (rsite3[0], rsite3[1], rsite3[2], latgd * rad2deg, lon * rad2deg,
          alt * 1000))
year, mon, day, hr, min, sec = stu.invjday(jd1, jdf1)
print('obs#1 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
      % (year, mon, day, hr, min, sec, rtasc1 * rad2deg, decl1 * rad2deg))
year, mon, day, hr, min, sec = stu.invjday(jd2, jdf2)
print('obs#2 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
      % (year, mon, day, hr, min, sec, rtasc2 * rad2deg, decl2 * rad2deg))
year, mon, day, hr, min, sec = stu.invjday(jd3, jdf3)
print('Obs#3 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
      % (year, mon, day, hr, min, sec, rtasc3 * rad2deg, decl3 * rad2deg))
print('\nsolution by %s \n\n' % (processtype))
if r2.all(): print('r2     %11.7f   %11.7f  %11.7f er    %11.7f  %11.7f  '
                   '%11.7f km \n'
                   % (r2[0] / re, r2[1] / re, r2[2] / re, r2[0], r2[1], r2[2]))
#        fprintf(1, 'r2 ans #11.7f   #11.7f  #11.7f er    #11.7f  #11.7f  #11.7f km \n', r2ans/re, r2ans);

if v2.all(): print('v2     %11.7f   %11.7f  %11.7f er/tu %11.7f  %11.7f  '
                   '%11.7f km/s\n'
                   % (v2[0]/velkmps, v2[1]/velkmps, v2[2]/velkmps, v2[0],
                      v2[1], v2[2]))
#        fprintf(1, 'v2 ans #11.7f   #11.7f  #11.7f er/tu #11.7f  #11.7f  #11.7f km/s\n', v2ans/velkmps, v2ans);

if r2.all() and v2.all():
  p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = \
    sc.rv2coeh(r2, v2, re, mu)
  print('         p km          a km         ecc       incl deg     '
        'raan deg    argp deg     nu deg      m deg  \n')
  if p: print('coes p=%11.4f\n' % (p))
  if a: print('coes a=%11.4f\n' % (a))
  if ecc: print('coes ecc=%13.9f\n' % (ecc))
  if incl: print('coes incl=%13.7f\n' % (incl * rad2deg))
  if omega: print('coes omega=%11.5f\n' % (omega * rad2deg))
  if argp: print('coes argp=%11.5f\n' % (argp * rad2deg))
  if nu: print('coes nu=%11.5f\n' % (nu * rad2deg))
  if m: print('coes m=%11.5f\n' % (m * rad2deg))
  print('%s \n' % (answ))
#         if typerun == 'o' || typerun == 'a'
#             [r2, v2] = anglesgood( decl1, decl2, decl3, rtasc1, rtasc2, ...
#                 rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, rsite1, rsite2, rsite3 );
#             processtype = 'anglesgood';
#         end

#         # -------------- write out answer --------------
#         fprintf(1, '\n\ninputs: \n\n');
#         [latgc, latgd, lon, alt] = ecef2ll ( rs1 ); # need to use ecef one!!
#         fprintf(1, 'Site obs1 #11.7f #11.7f #11.7f km  lat #11.7f lon #11.7f alt #11.7f  \n', rsite1[0], rsite1[1], rsite1[2], latgd * rad2deg, lon * rad2deg, alt*1000 );
#         [latgc, latgd, lon, alt] = ecef2ll ( rs2 );
#         fprintf(1, 'Site obs2 #11.7f #11.7f #11.7f km  lat #11.7f lon #11.7f alt #11.7f  \n', rsite2[0], rsite2[1], rsite2[2], latgd * rad2deg, lon * rad2deg, alt*1000 );
#         [latgc, latgd, lon, alt] = ecef2ll ( rs3 );
#         fprintf(1, 'Site obs3 #11.7f #11.7f #11.7f km  lat #11.7f lon #11.7f alt #11.7f  \n', rsite3[0], rsite3[1], rsite3[2], latgd * rad2deg, lon * rad2deg, alt*1000 );
#         [year, mon, day, hr, min, sec] = invjday ( jd1, jdf1 );
#         fprintf(1, 'obs#1 #4i #2i #2i #2i #2i #6.3f ra #11.7f de #11.7f  \n', year, mon, day, hr, min, sec, rtasc1 * rad2deg, decl1 * rad2deg );
#         [year, mon, day, hr, min, sec] = invjday ( jd2, jdf2 );
#         fprintf(1, 'obs#2 #4i #2i #2i #2i #2i #6.3f ra #11.7f de #11.7f  \n', year, mon, day, hr, min, sec, rtasc2 * rad2deg, decl2 * rad2deg );
#         [year, mon, day, hr, min, sec] = invjday ( jd3, jdf3 );
#         fprintf(1, 'Obs#3 #4i #2i #2i #2i #2i #6.3f ra #11.7f de #11.7f  \n', year, mon, day, hr, min, sec, rtasc3 * rad2deg, decl3 * rad2deg );

#         fprintf(1, '\nsolution by #s \n\n', processtype);
#         fprintf(1, 'r2     #11.7f   #11.7f  #11.7f er    #11.7f  #11.7f  #11.7f km \n', r2/re, r2);
#         #        fprintf(1, 'r2 ans #11.7f   #11.7f  #11.7f er    #11.7f  #11.7f  #11.7f km \n', r2ans/re, r2ans);

#         fprintf(1, 'v2     #11.7f   #11.7f  #11.7f er/tu #11.7f  #11.7f  #11.7f km/s\n', v2/velkmps, v2);
#         #        fprintf(1, 'v2 ans #11.7f   #11.7f  #11.7f er/tu #11.7f  #11.7f  #11.7f km/s\n', v2ans/velkmps, v2ans);

#         [p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper ] = rv2coeh (r2, v2, re, mu);
#         fprintf(1, '         p km          a km         ecc       incl deg     raan deg    argp deg     nu deg      m deg  \n');
#         fprintf(1, 'coes #11.4f #11.4f #13.9f #13.7f #11.5f #11.5f #11.5f #11.5f \n', ...
#             p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg );

#         fprintf(1, '#s \n', answ );

#         [p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper ] = rv2coeh (r2ans, v2ans, re, mu);
#         fprintf(1, '         p km          a km         ecc       incl deg     raan deg    argp deg     nu deg      m deg  \n');
#         fprintf(1, 'coes #11.4f #11.4f #13.9f #13.7f #11.5f #11.5f #11.5f #11.5f \n', ...
#             p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg );

range_ = 7550.679305
r = 1000.0 * np.array([7.721359586705, -6.293226594912, -1.333676254902])
magr = smu.mag(r)
print('r:')
print(r)
print('magr:')
print(magr)
rho = 1000.0 * np.array([4.648325989705, -2.479803066912, -5.409004552902])
magrho = smu.mag(rho)
print('rho:')
print(rho)
print('magrho:')
print(magrho)
rs = 1000.0 * np.array([3.073033597, -3.813423528, 4.075328298])
magrs = smu.mag(rs)
print('rs:')
print(rs, magrs)
print('magrs')
print(magrs)
los1 = np.array([0.6156119, -0.3284185, -0.7163541])
print('los1:')
print(los1)
print('mag los1:')
print(smu.mag(los1))
rtasc1 = 331.922 * deg2rad
rx = np.arctan(rho[1] / rho[0]) * rad2deg + 360
decl1 = -45.754 * deg2rad
dx = np.arcsin(rho[2] / smu.mag(rho)) * rad2deg
print('rtasc %11.7f  %11.7f decl %11.7f  %11.7f \n'
      % (rtasc1 * rad2deg, rx, decl1 * rad2deg, dx))
los1 = np.array([np.cos(decl1) * np.cos(rtasc1),
                 np.cos(decl1) * np.sin(rtasc1), np.sin(decl1)])
print('los1:')
print(los1)
print('mag los1:')
print(smu.mag(los1))
rtasc1x = np.arctan(rho[1] / rho[0])
los1x = np.array([np.cos(decl1) * np.cos(rtasc1x),
                  np.cos(decl1) * np.sin(rtasc1x), np.sin(decl1)])
print('los1x:')
print(los1x)
print('mag los1x:')
print(smu.mag(los1x))
urho = smu.unit(rho)
print('urho:')
print(urho)
print('mag urho:')
print(smu.mag(urho))
dotrsl2 = 2.0 * np.dot(rs, np.transpose(los1))
rhotem1 = 0.5 * (-dotrsl2 + np.sqrt(dotrsl2 * dotrsl2 - 4.0 *
                                    (magrs * magrs - magr * magr)))
r1 = rhotem1 * los1 + rs
print('r1:')
print(r1)
print('mag r1:')
print(smu.mag(r1))
rhotem1 = 0.5 * (-dotrsl2 - np.sqrt(dotrsl2 * dotrsl2 - 4.0 *
                                     (magrs * magrs - magr * magr)))
r1 = rhotem1 * los1 + rs
print('r1:')
print(r1)
print('mag r1:')
print(smu.mag(r1))
rr = smu.unit(rho)
x1 = np.array([0.0, los1[0], 0.0, rr[0]])
y1 = np.array([0.0, los1[1], 0.0, rr[1]])
z1 = np.array([0.0, los1[2], 0.0, rr[2]])
#plot3(x1, y1, z1, 'X')
mjd1 = 58896.0
rs1 = 124647954.923
mjd2 = 58897.0
rs2 = 126054577.073
mjd3 = 58898.0
rs3 = 127422565.338
mjd4 = 58899.0
rs4 = 128751471.066
smu.cubicinterp(rs1, rs2, rs3, rs4, mjd1, mjd2, mjd3, mjd4, 58897.0416667)
