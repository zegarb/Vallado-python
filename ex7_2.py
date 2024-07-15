#     -----------------------------------------------------------------
#
#                              Ex7_2
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
from space_constants import sethelp as sh
import space_conversions as sc
import spacemath_utils as smu
import spacetime_utils as stu
import os

answ = '# 0  ans a = 12246.023 km   ecc = 0.2000 deg    incl = 40.00 deg    raan = 330.000 deg'
#typerun = 'l'; # laplace
#typerun = 'g'; # gauss
#typerun = 'd'; # doubler
typerun = 'a' #run them all!


fn = os.path.join(os.path.dirname(__file__), "data",
                       "Sat11Ex" + str(0) + ".dat")
filedat = np.loadtxt(fn)

if sh.show:
    print(filedat)

#r2ans = np.array([5897.954130507, 5791.046114526, 6682.733686585])
#v2ans = np.array([-4.393910234, 4.576816355, 1.482423676])

dut1 = - 0.609641
dat = 35
xp = 0.137495 * arcsec2rad
yp = 0.342416 * arcsec2rad
lod = 0.0
timezone = 0
terms = 0
ddpsi = 0.0
ddeps = 0.0

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

if sh.show:
    print('ECEF Site Vectors:')
    print('rs1 and vs1:')
    print(rs1)
    print(vs1)
    print('rs2 and vs2:')
    print(rs2)
    print(vs2)
    print('rs3 and vs3:')
    print(rs3)
    print(vs3)
    print('\n')

year, mon, day, hr, min, second = stu.invjday(obsrecarr[obs1]['jd'],
                                              obsrecarr[obs1]['jdf'])
utc = second
ut1 = utc + dut1
tai = utc + dat
tt = tai + 32.184
jdut1, jdut1frac = stu.jday(year, mon, day, hr, min, ut1)
jdtt, jdttfrac = stu.jday(year, mon, day, hr, min, tt)
ttt = (jdtt - 2451545.0) / 36525.0
print('initial site time/data')
print('year %5i ' % (year))
print('mon %4i ' % (mon))
print('day %3i ' % (day))
print('hr %3i:%2i:%8.6f\n' % (hr, min, second))
print('dut1 %8.6f s' % (dut1))
print(' dat %3i s' % (dat))
print(' xp %8.6f "' % (xp * rad2arcsec))
print(' yp %8.6f "' % (yp * rad2arcsec))
print(' lod %8.6f s\n' % (lod))

# -------------- convert each site vector from ecef to eci -----------------
a = np.array([[0], [0], [0]])

year, mon, day, hr, min, sec = stu.invjday(jd1, jdf1)
(ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb,
 ttdb, jdtdb, jdtdbfrac) = stu.convtime(year, mon, day, hr, min, sec,
                                        timezone, dut1, dat)
rsite1, _, _ = sc.ecef2eci(rs1, vs1, a, ttt, jdut1 + jdut1frac, lod, xp,
                                   yp, 2, ddpsi, ddeps)
year, mon, day, hr, min, sec = stu.invjday(jd2, jdf2)
(ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb,
 jdtdb, jdtdbfrac) = stu.convtime(year, mon, day, hr, min, sec, timezone, dut1,
                                  dat)
rsite2, _, _ = sc.ecef2eci(rs2, vs2, a, ttt, jdut1 + jdut1frac, lod, xp,
                                  yp, 2, ddpsi, ddeps)
year, mon, day, hr, min, sec = stu.invjday(jd3, jdf3)
(ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb,
 jdtdb, jdtdbfrac) = stu.convtime(year, mon, day, hr, min, sec, timezone, dut1,
                                  dat)
rsite3, _, _ = sc.ecef2eci(rs3, vs3, a, ttt, jdut1 + jdut1frac, lod, xp,
                                  yp, 2, ddpsi, ddeps)


# ---------------------- run the angles-only routine ------------------
# ---- anglesl test ----
if typerun == 'l' or typerun == 'a':
    r2, v2 = obu.anglesl(decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1,
                         jdf1, jd2, jdf2, jd3, jdf3, rsite1)
    #                          rtasc3, jd1, jd2, jd3, rs1, rs2, rs3, re, mu, tu );
    processtype = 'anglesl'
    # -------------- write out answer --------------
    print('\n\ninputs: \n\n' % ())
    latgc, latgd, lon, alt = sc.ecef2ll(rs1)

    print('ECI Site Vectors and lat/long/alt:\n')
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

    print('\nObservation Angles:\n')
    year, mon, day, hr, min, sec = stu.invjday(jd1, jdf1)
    print('Obs#1 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
        % (year, mon, day, hr, min, sec, rtasc1 * rad2deg, decl1 * rad2deg))
    year, mon, day, hr, min, sec = stu.invjday(jd2, jdf2)
    print('Obs#2 %4i %2i %2i %2i %2i %6.3f ra %11.7f de %11.7f  \n'
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
        sc.rv2coe(r2, v2, mu)
    print('          p km          a km         ecc             incl deg      '
        'raan deg        argp deg     nu deg       m deg  \n')
    print('         %11.4f  %11.4f  %13.9f  %13.7f  %11.5f    '
        '%11.5f  %11.5f  %11.5f  \n' % (p, a, ecc, incl*rad2deg, omega*rad2deg,
                                        argp*rad2deg, nu*rad2deg, m*rad2deg))
    print('Case %s \n' % (answ))

# ---- anglesg test ----
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
        sc.rv2coe(r2, v2, mu)
    print('          p km          a km         ecc             incl deg      '
        'raan deg        argp deg     nu deg       m deg  \n')
    print('         %11.4f  %11.4f  %13.9f  %13.7f  %11.5f    '
        '%11.5f  %11.5f  %11.5f  \n' % (p, a, ecc, incl*rad2deg, omega*rad2deg,
                                        argp*rad2deg, nu*rad2deg, m*rad2deg))
    print('Case %s \n' % (answ))

# ---- anglesdr test ----
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
        sc.rv2coe(r2, v2, mu)
    print('          p km          a km         ecc             incl deg      '
        'raan deg        argp deg     nu deg       m deg  \n')
    print('         %11.4f  %11.4f  %13.9f  %13.7f  %11.5f    '
        '%11.5f  %11.5f  %11.5f  \n' % (p, a, ecc, incl*rad2deg, omega*rad2deg,
                                        argp*rad2deg, nu*rad2deg, m*rad2deg))
    print('Case %s \n' % (answ))