import math
import numpy as np

import spacetime_utils as stu
import spacemath_utils as smu
import space_conversions as sc
import orbit_utils as obu
from space_constants import *



#     -----------------------------------------------------------------
#
#                              Ex5_1.m
#
#  this file demonstrates example 5-1.
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

print('------- example 5-1 ---------')
timezone = 0

year = 2006  # need UTC that will give TDT on the 2 Apr 0 hr
mon = 4
day = 1
hr = 23
minute = 58
second = 54.816
jd, jdfrac = stu.jday(year, mon, day, hr, minute, second)
print('jd  %11.9f \n' % (jd+jdfrac))
dat = 33
dut1 = 0.2653628

ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, \
    jdtdb, jdtdbfrac = stu.convtime(year, mon, day, hr, minute, second,
                                    timezone, dut1, dat)
print('input data')
print(' year %5i ' % year)
print(' mon %4i ' % mon)
print(' day %3i ' % day)
print(' %3i:%2i:%8.6f\n ' % (hr, minute, second))
print(' dut1 %8.6f s' % dut1)
print(' dat %3i s' % dat)

print('tt  %8.6f ttt  %16.12f jdtt  %18.11f ' % (tt, ttt, jdtt))
h, m, s = stu.sec2hms(tt)
print('hms %3i %3i %8.6f \n' % (h, m, s))

rsun, rtasc, decl = obu.sun(jd + jdfrac)
print('sun  rtasc %14.6f deg decl %14.6f deg\n' % (rtasc * rad2deg, decl * rad2deg))
print('sun newTOD %11.9f%11.9f%11.9f au\n' % (rsun[0], rsun[1], rsun[2]))
print('sun newTOD %14.4f%14.4f%14.4f km\n' % (rsun[0]*au, rsun[1]*au,
                                              rsun[2]*au))

rsunaa = np.array([0.9775113, 0.1911521, 0.0828717]) * au # astronomical alm value into km
print('rs almanac  ICRF %11.9f %11.9f %11.9f km \n'
      % (rsunaa[0], rsunaa[1], rsunaa[2]))
#rs aa ICRF 146233608.380930990 28595947.006026998 12397429.803279001 km
# from jpl de430 (icrf) ephemerides are centered on jdtdb
#   2006  4  2          146233601.8528      28595943.3080      12397428.8371     149518194.6330  1.0010661
#   2006  4  2           187657.2460        286860.8848        155666.9298
print('rs jplde430 ICRF 146233601.8528      28595943.3080      12397428.8371   km \n\n')

print('now take the TOD sun vector (analytical) and move to MOD and TOD for comparison \n')
#ttt= ( jd - 2451545.0  )/ 36525.0
vmod = np.zeros((3)).T
amod = np.zeros((3)).T
reci, veci, aeci = sc.mod2eci(rsun.T, vmod, amod, ttt)
print("reci:")
print(reci)
print("veci:")
print(veci)
print("aeci:")
print(aeci)

#print('mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (db[0], db[1], db[2], smu.mag(db) ))


ddpsi=0.0
ddeps=0.0

reci, veci, aeci = sc.tod2eci(rsun.T, vmod, amod, ttt, ddpsi, ddeps)
#print('tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n' % (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci * au - rsunaa.T
#print('delta tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n\n' % (db[0], db[1], db[2], smu.mag(db) ))

#         [reci, veci, aeci] = eci2mod ( rsun.T, vmod, amod, ttt )
#         print('eci - mod  %11.9f %11.9f %11.9f %11.4f km  \n'%reci*au, smu.mag(reci)*au )
#         db = reci*au-rsunaa.T
#         print('delta eci - mod  %11.9f %11.9f %11.9f %11.4f km  \n'%db, smu.mag(db) )
#
#         [reci, veci, aeci] = eci2tod ( rsun.T, vmod, amod, ttt, ddpsi, ddeps )
#         print('eci - tod  %11.9f %11.9f %11.9f %11.4f km  \n'%reci*au, smu.mag(reci)*au )
#         db = reci*au-rsunaa.T
#         print('delta eci - tod  %11.9f %11.9f %11.9f %11.4f km  \n'%db, smu.mag(db) )

hms = sc.hms2rad(0, 44, 33.42)
dms = sc.dms2rad(4, 47, 18.3)
print('hms ast alm rtasc %11.9f decl %11.9f \n'% (hms * rad2deg, dms * rad2deg ))


# now try alamnac method
rsuna, rtasca, decla = obu.sunalmanac (jd+jdfrac)
print('\n\nsun  rtasc %14.6f deg decl %14.6f deg\n' % (rtasca * rad2deg, decla * rad2deg ))
print('sun ALM %11.9f%11.9f%11.9f au\n' % (rsuna[0], rsuna[1], rsuna[2]))
print('sun ALM %14.4f%14.4f%14.4f km\n' % (rsuna[0]*au, rsuna[1]*au, rsuna[2]*au ))

print('rs aa ICRF %11.9f %11.9f %11.9f km \n' % (rsuna[0], rsuna[1], rsuna[2]))

print('now take the TOD sun vector (analytical) and move to MOD and TOD for comparison \n')
#ttt= ( jd - 2451545.0  )/ 36525.0
vmod = np.zeros((3)).T
amod = np.zeros((3)).T
reci, veci, aeci = sc.mod2eci(rsuna.T, vmod, amod, ttt)
print("reci:")
print(reci)
print("veci:")
print(veci)
print("aeci:")
print(aeci)


# print('mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
# print('delta mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (db[0], db[1], db[2], smu.mag(db) ))

reci, veci, aeci = sc.tod2eci(rsuna.T, vmod, amod, ttt, ddpsi, ddeps)
#print('tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n' % (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n\n' % (db[0], db[1], db[2], smu.mag(db) ))

#         [reci, veci, aeci] = eci2mod ( rsuna.T, vmod, amod, ttt )
#         print('eci - mod  %11.9f %11.9f %11.9f %11.4f km  \n'%reci*au, smu.mag(reci)*au )
#         db = reci*au-rsunaa.T
#         print('delta eci - mod  %11.9f %11.9f %11.9f %11.4f km  \n'%db, smu.mag(db) )
#
#         [reci, veci, aeci] = eci2tod ( rsuna.T, vmod, amod, ttt, ddpsi, ddeps )
#         print('eci - tod  %11.9f %11.9f %11.9f %11.4f km  \n'%reci*au, smu.mag(reci)*au )
#         db = reci*au-rsunaa.T
#         print('delta eci - tod  %11.9f %11.9f %11.9f %11.4f km  \n'%db, smu.mag(db) )

hms = sc.hms2rad( 0, 44, 33.42 )
dms = sc.dms2rad( 4, 47, 18.3 )
print('hms ast alm rtasc %11.9f decl %11.9f \n'%(hms * rad2deg, dms * rad2deg) )

print('==============================================================\n')
# previous edition example
year = 1994
mon = 4
day = 1
hr = 23
minute = 58
second = 59.816
jd, jdfrac = stu.jday(year, mon, day, hr, minute, second)
print('jd  %11.9f \n'%(jd+jdfrac) )
dat = 28
dut1 = -0.0226192
ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac  \
    = stu.convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat )
print('input data \n\n')
print(' year %5i '%year)
print(' mon %4i '%mon)
print(' day %3i '%day)
print(' %3i:%2i:%8.6f\n '%(hr, minute, second) )
print(' dut1 %8.6f s'%dut1)
print(' dat %3i s'%dat)

print('tt  %8.6f ttt  %16.12f jdtt  %18.11f '%(tt, ttt, jdtt) )
h, m, s = stu.sec2hms( tt )
print('hms %3i %3i %8.6f \n'%(h, m, s))

rsun, rtasc, decl = obu.sun ( jd+jdfrac )
print('sun  rtasc %14.6f deg decl %14.6f deg\n'%(rtasc * rad2deg, decl * rad2deg) )
print('sun ICRS %11.9f%11.9f%11.9f au\n' % (rsun[0], rsun[1], rsun[2]))
print('sun ICRS %14.4f%14.4f%14.4f km\n' % (rsun[0]*au, rsun[1]*au, rsun[2]*au ))

rsunaa = np.array([0.9772766, 0.1922635, 0.0833613])*au # astronomical alm value into km
print('rs almanac MOD %11.9f %11.9f %11.9f km \n' % (rsunaa[0], rsunaa[1], rsunaa[2]))

print('now take the TOD sun vector (analytical) and move to MOD and TOD for comparison \n')
#ttt= ( jd - 2451545.0  )/ 36525.0
vmod = np.zeros((3)).T
amod = np.zeros((3)).T
reci, veci, aeci = sc.mod2eci  ( rsun.T, vmod, amod, ttt )
print("reci:")
print(reci)
print("veci:")
print(veci)
print("aeci:")
print(aeci)


#print('mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (db[0], db[1], db[2], smu.mag(db) ))

reci, veci, aeci = sc.tod2eci  ( rsun.T, vmod, amod, ttt, ddpsi, ddeps )
#print('tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n\n' % (db[0], db[1], db[2],  smu.mag(db) ))

#         [reci, veci, aeci] = eci2mod ( rsun.T, vmod, amod, ttt )
#         print('eci - mod  %11.9f %11.9f %11.9f %11.4f km  \n'%reci*au, smu.mag(reci)*au )
#         db = reci*au-rsunaa.T
#         print('delta eci - mod  %11.9f %11.9f %11.9f %11.4f km  \n'%db, smu.mag(db) )
#
#         [reci, veci, aeci] = eci2tod ( rsun.T, vmod, amod, ttt, ddpsi, ddeps )
#         print('eci - tod  %11.9f %11.9f %11.9f %11.4f km  \n'%reci*au, smu.mag(reci)*au )
#         db = reci*au-rsunaa.T
#         print('delta eci - tod  %11.9f %11.9f %11.9f %11.4f km  \n'%db, smu.mag(db) )

# now try alamnac method
rsuna, rtasca, decla = obu.sunalmanac ( jd+jdfrac )
print('\n\nsun  rtasc %14.6f deg decl %14.6f deg\n'%(rtasca * rad2deg, decla * rad2deg) )
print('sun ALM %11.9f%11.9f%11.9f au\n' % (rsuna[0], rsuna[1], rsuna[2]))
print('sun ALM %14.4f%14.4f%14.4f km\n' % (rsuna[0]*au, rsuna[1]*au, rsuna[2]*au ))

print('rs aa ICRF %11.9f %11.9f %11.9f km \n' % (rsunaa[0], rsunaa[1], rsunaa[2]))

print('now take the TOD sun vector (analytical) and move to MOD and TOD for comparison \n')
#ttt= ( jd - 2451545.0  )/ 36525.0
vmod = np.zeros((3)).T
amod = np.zeros((3)).T
reci, veci, aeci = sc.mod2eci  ( rsuna.T, vmod, amod, ttt )
print("reci:")
print(reci)
print("veci:")
print(veci)
print("aeci:")
print(aeci)


#print('mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (db[0], db[1], db[2], smu.mag(db) ))

reci, veci, aeci = sc.tod2eci  ( rsuna.T, vmod, amod, ttt, ddpsi, ddeps )
#print('tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n' % (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n' % (db[0], db[1], db[2], smu.mag(db) ))

print('==============================================================\n')
# another example tdt = 29+32.184 secs less than 4/2 at 0 hrs
year = 1995
mon = 4
day = 1
hr = 23
minute = 58
second = 58.816
jd, jdfrac = stu.jday(year, mon, day, hr, minute, second)
print('jd  %11.9f \n'%(jd+jdfrac) )
dat = 29
dut1 = 0.1535663
ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac  \
    = stu.convtime ( year, mon, day, hr, minute, second, timezone, dut1, dat )
print('input data \n\n')
print(' year %5i '%year)
print(' mon %4i '%mon)
print(' day %3i '%day)
print(' %3i:%2i:%8.6f\n '%(hr, minute, second) )
print(' dut1 %8.6f s'%dut1)
print(' dat %3i s'%dat)

print('ut1 %8.6f tut1 %16.12f jdut1 %18.11f '%(ut1, tut1, jdut1) )
h, m, s = stu.sec2hms( ut1 )
print('hms %3i %3i %8.6f \n'%(h, m, s))
print('utc %8.6f '%utc )
h, m, s = stu.sec2hms( utc )
print('hms %3i %3i %8.6f \n'%(h, m, s))
print('tai %8.6f'%tai )
h, m, s = stu.sec2hms( tai )
print('hms %3i %3i %8.6f \n'%(h, m, s))
print('tt  %8.6f ttt  %16.12f jdtt  %18.11f '%(tt, ttt, jdtt) )
h, m, s = stu.sec2hms( tt )
print('hms %3i %3i %8.6f \n'%(h, m, s))
print('tdb %8.6f ttdb %16.12f jdtdb %18.11f\n'%(tdb, ttdb, jdtdb) )

rsun, rtasc, decl = obu.sun ( jd+jdfrac )
print('sun  rtasc %14.6f deg decl %14.6f deg\n'%(rtasc * rad2deg, decl * rad2deg) )
print('sun ICRS %11.9f%11.9f%11.9f au\n' % (rsun[0], rsun[1], rsun[2]))
print('sun ICRS %14.4f%14.4f%14.4f km\n' % (rsun[0]*au, rsun[1]*au, rsun[2]*au ))

rsunaa = np.array([0.9781158, 0.1884327, 0.0816997])*au # astronomical alm value into km
print('rs almanac MOD %11.9f %11.9f %11.9f km \n' % (rsunaa[0], rsunaa[1], rsunaa[2]))

#ttt= ( jd - 2451545.0  )/ 36525.0
vmod = np.zeros((3)).T
amod = np.zeros((3)).T
reci, veci, aeci = sc.mod2eci  ( rsun.T, vmod, amod, ttt )
print("reci:")
print(reci)
print("veci:")
print(veci)
print("aeci:")
print(aeci)


#print('mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (db[0], db[1], db[2], smu.mag(db) ))

reci, veci, aeci = sc.tod2eci  ( rsun.T, vmod, amod, ttt, ddpsi, ddeps )
#print('tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n\n' % (db[0], db[1], db[2],  smu.mag(db) ))

#         [reci, veci, aeci] = eci2mod ( rsun.T, vmod, amod, ttt )
#         print('eci - mod  %11.9f %11.9f %11.9f %11.4f km  \n'%reci*au, smu.mag(reci)*au )
#         db = reci*au-rsunaa.T
#         print('delta eci - mod  %11.9f %11.9f %11.9f %11.4f km  \n'%db, smu.mag(db) )
#
#         [reci, veci, aeci] = eci2tod ( rsun.T, vmod, amod, ttt, ddpsi, ddeps )
#         print('eci - tod  %11.9f %11.9f %11.9f %11.4f km  \n'%reci*au, smu.mag(reci)*au )
#         db = reci*au-rsunaa.T
#         print('delta eci - tod  %11.9f %11.9f %11.9f %11.4f km  \n'%db, smu.mag(db) )

# now try alamnac method
rsuna, rtasca, decla = obu.sunalmanac ( jd+jdfrac )
print('\n\nsun  rtasc %14.6f deg decl %14.6f deg\n'%(rtasca * rad2deg, decla * rad2deg) )
print('sun ALM %11.9f%11.9f%11.9f au\n' % (rsuna[0], rsuna[1], rsuna[2]))
print('sun ALM %14.4f%14.4f%14.4f km\n' % (rsuna[0]*au, rsuna[1]*au, rsuna[2]*au ))

print('rs aa ICRF %11.9f %11.9f %11.9f km \n' % (rsunaa[0], rsunaa[1], rsunaa[2]))

print('now take the TOD sun vector (analytical) and move to MOD and TOD for comparison \n')
#ttt= ( jd - 2451545.0  )/ 36525.0
vmod = np.zeros((3)).T
amod = np.zeros((3)).T
reci, veci, aeci = sc.mod2eci  ( rsuna.T, vmod, amod, ttt )
print("reci:")
print(reci)
print("veci:")
print(veci)
print("aeci:")
print(aeci)


#print('mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta mod - eci  %11.9f %11.9f %11.9f %11.4f km  \n'% (db[0], db[1], db[2], smu.mag(db) ))

reci, veci, aeci = sc.tod2eci  ( rsuna.T, vmod, amod, ttt, ddpsi, ddeps )
#print('tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n' % (reci[0]*au, reci[1]*au, reci[2]*au, smu.mag(reci)*au ))
db = reci*au-rsunaa.T
#print('delta tod - eci  %11.9f %11.9f %11.9f %11.4f km  \n' % (db[0], db[1], db[2], smu.mag(db) ))


