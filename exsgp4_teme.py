import os
import numpy as np
import math
import orbit_utils as obu
from space_constants import *
import space_conversions as sc
import spacemath_utils as smu
import spacetime_utils as stu


#
#  exsgp4_teme
#
#  tests the teme conversions including of date and of epoch
#

# ---------- test the teme to pef conversion ------------
recef = np.array([[- 1033.479383],[7901.2952754],[6380.3565958]])
vecef = np.array([[- 3.22563652],[- 2.87245145],[5.531924446]])
aecef = np.array([[1],[1],[1]])
year = 2004
mon = 4
day = 6
hr = 7
min = 51
sec = 28.386009
dut1 = - 0.4399619
dat = 32
xp = - 0.140682 * arcsec2rad

yp = 0.333309 * arcsec2rad
lod = 0.0015563

ddpsi = - 0.052195 * arcsec2rad

ddeps = - 0.003875 * arcsec2rad
ddx = - 0.000205 * arcsec2rad

ddy = - 0.000136 * arcsec2rad
timezone = 0
order = 106
eqeterms = 2

print('test program for reduction functions \n\n')
print('input data \n\n')
print('year %5i ' % (year))
print(' mon %4i ' % (mon))
print(' day %3i ' % (day))
print('%3i:%2i:%8.6f\n ' % (hr,min,sec))
print('dut1 %8.6f s' % (dut1))
print(' dat %3i s' % (dat))
print(' xp %8.6f "' % (xp * rad2arcsec))
print(' yp %8.6f "' % (yp * rad2arcsec))
print(' lod %8.6f s\n' % (lod))
print(' ddpsi %8.6f " ddeps  %8.6f\n' % (ddpsi * rad2arcsec,ddeps * rad2arcsec))
print(' ddx %8.6f " ddy  %8.6f\n' % (ddx * rad2arcsec,ddy * rad2arcsec))
print('order %3i  eqeterms %3i  \n' % (order,eqeterms))
print('units are km and km/s \n')
timezone = 0
# -------- convtime    - convert time from utc to all the others
print('convtime results\n')
# , tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
ut1,tut1,jdut1,jdut1frac,utc,tai,tt,ttt,jdtt,jdttfrac,tdb,ttdb,jdtdb,jdtdbfrac = stu.convtime(year,mon,day,hr,min,sec,timezone,dut1,dat)
print('ut1 %8.6f tut1 %16.12f jdut1 %18.11f \n' % (ut1,tut1,jdut1))
# ---------------------- teme transformations ---------------------
# -------- ecef2teme    - transform ecef to teme vectors
rteme,vteme,ateme = sc.ecef2teme(recef,vecef,aecef,ttt,jdut1 + jdut1frac,lod,xp,yp,eqeterms)
print('\n\n start from ecef \n')
print(f'ecef-teme\n rteme \n{rteme}')
print(f' vteme \n{vteme}')
# -------- teme2ecef    - transform teme to ecef vectors
recef1,vecef1,aecef1 = sc.teme2ecef(rteme,vteme,ateme,ttt,jdut1 + jdut1frac,lod,xp,yp,eqeterms)
print(f'teme-ecef\n recef \n{recef1}')
print(f' vecef \n{vecef1}\n')
dr = 1000 * (recef - recef1)

print('diff in ecef  %14.7f %14.7f %14.7f %14.7f \n' % (dr[0,0],dr[1,0],dr[2,0],smu.mag(dr)))
# ---------------------- eci transformations ----------------------
# -------- teme2eci    - transform teme to eci vectors
print(f'\n\n start from teme \n')
reci,veci,aeci = sc.teme2eci(rteme,vteme,ateme,ttt,ddpsi,ddeps)
print(f'teme-eci\n reci \n{reci}')
print(f' veci \n{veci}\n')
# -------- eci2teme    - transform eci to teme vectors
rteme1,vteme1,ateme1 = sc.eci2teme(reci,veci,aeci,ttt,ddpsi,ddeps)
print(f'ecef-teme\n rteme \n{rteme1}')
print(f' vteme \n{vteme1}\n')
dr = 1000 * (rteme - rteme1)

print('diff in teme  %14.7f %14.7f %14.7f %14.7f \n' % (dr[0,0],dr[1,0],dr[2,0],smu.mag(dr)))
# check standard ecef to eci transformation
recig,vecig,aecig = sc.ecef2eci(recef,vecef,aecef,ttt,jdut1 + jdut1frac,lod,xp,yp,2,ddpsi,ddeps)
print(f'GCRF 2 w corr IAU-76/FK5   \n{recig}')
print(f' v \n{vecig}\n')
dr = 1000 * (reci - recig)

print('diff in eci  %14.7f %14.7f %14.7f %14.7f \n' % (dr[0,0],dr[1,0],dr[2,0],smu.mag(dr)))

# -------------------- teme of epoch/date --------------------
print('\n\n =================== now do teme of date and of epoch example  ================= \n')
#       //typerun = 'c' compare 1 year of full satcat data
#       //typerun = 'v' verification run, requires modified elm file with
#       //              start stop and delta times
rad = 180.0 / np.pi
opsmode = 'a'

whichconst = 72
eqeterms = 2
longstr1 = '1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753'
longstr2 = '2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667     0.0      4320.0       720.00000'
#       // convert the char string to sgp4 elements
#       // includes initialization of sgp4
typerun = 'v'
startmfe,stopmfe,deltamin,satrec = sc.twoline2rv(longstr1,longstr2,typerun,'e',opsmode,whichconst)
print('%11.7f %11.7f %11.7f \n' % (startmfe,stopmfe,deltamin))
print(' %d\n' % (satrec['satnum']))
#      // call the propagator to get the initial state vector value
satrec,ro,vo = obu.sgp4(satrec,0.0)
print(' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n' % (satrec['t'],ro[0],ro[1],ro[2],vo[0],vo[1],vo[2]))
tsince = startmfe
#  // check so the first value isn't written twice
if (np.abs(tsince) > 1e-08):
    tsince = tsince - deltamin

# eop data
# date      MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
#2000 06 27 51722  0.109021  0.285481  0.2052872  0.0004465 -0.053678 -0.006320 -0.000025  0.000028  32
year,mon,day,hr,min,sec = stu.invjday(satrec['jdsatepoch'],satrec['jdsatepochf'])
dut1 = 0.2052872
dat = 32
xp = 0.109021 * arcsec2rad

yp = 0.285481 * arcsec2rad
lod = 0.0004465

ddpsi = - 0.053678 * arcsec2rad

ddeps = - 0.00632 * arcsec2rad
ddx = - 2.5e-05 * arcsec2rad

ddy = 2.8e-05 * arcsec2rad
#, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
ut1,tut1,jdut1,jdut1frac,utc,tai,tt,ttt,jdtt,jdttfrac,tdb,ttdb,jdtdb,jdtdbfrac = stu.convtime(year,mon,day,hr,min,sec,timezone,dut1,dat)
print('year %5i ' % (year))
print(' mon %4i ' % (mon))
print(' day %3i ' % (day))
print('%3i:%2i:%8.6f\n' % (hr,min,sec))
# for teme of epoch, find the transformation matrix at t0
prec,psia,wa,ea,xa = obu.precess(ttt,'80')
deltapsie,trueepse,meanepse,omegae,nute = obu.nutation(ttt,ddpsi,ddeps)

eqeg = deltapsie * np.cos(meanepse)
eqeg = np.fmod(eqeg,2.0 * np.pi)
eqe = np.zeros((3,3))
eqe[0,0] = np.cos(eqeg)
eqe[0,1] = np.sin(eqeg)
eqe[0,2] = 0.0
eqe[1,0] = - np.sin(eqeg)
eqe[1,1] = np.cos(eqeg)
eqe[1,2] = 0.0
eqe[2,0] = 0.0
eqe[2,1] = 0.0
eqe[2,2] = 1.0
tme = prec * nute * np.transpose(eqe)
# // loop to perform the propagation
# while ((tsince < stopmfe) && (satrec.error == 0))
while ((tsince < stopmfe)):

    tsince = tsince + deltamin
    if (tsince > stopmfe):
        tsince = stopmfe
    satrec,ro,vo = obu.sgp4(satrec,tsince)
    if (satrec['error'] == 0):
        if ((typerun != 'e') and (typerun != 'd')):
            jd = satrec['jdsatepoch']
            jdf = satrec['jdsatepochf'] + tsince / 1440.0
            year,mon,day,hr,min,sec = stu.invjday(jd,jdf)
            print(' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f  %5i%3i%3i %2i:%2i:%9.6f\n' % (satrec['t'],ro[0],ro[1],ro[2],vo[0],vo[1],vo[2],year,mon,day,hr,min,sec))
        else:
            print(' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f' % (tsince,ro[0],ro[1],ro[2],vo[0],vo[1],vo[2]))
            p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper = sc.rv2coe(ro,vo)
            print(' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f \n' % (a,ecc,incl * rad,omega * rad,argp * rad,nu * rad,m * rad))


rteme = np.array([[ro[0]],[ro[1]],[ro[2]]])
vteme = np.array([[vo[0]],[vo[1]],[vo[2]]])
print(' rteme \n', (rteme))
print(' vteme \n', (vteme))
# for teme of date, find the transformation matrix at t=t0 + 3days
#2000 06 30 51725  0.110329  0.281805  0.2042651  0.0001678 -0.054522 -0.006209  0.000002  0.000016  32
dut1 = 0.2042651
dat = 32
xp = 0.110329 * arcsec2rad

yp = 0.281805 * arcsec2rad
lod = 0.0001678

ddpsi = - 0.054522 * arcsec2rad

ddeps = - 0.006209 * arcsec2rad
ddx = 2e-06 * arcsec2rad

ddy = 1.6e-05 * arcsec2rad
#, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac
ut1,tut1,jdut1,jdut1frac,utc,tai,tt,ttt,jdtt,jdttfrac,tdb,ttdb,jdtdb,jdtdbfrac = stu.convtime(year,mon,day,hr,min,sec,timezone,dut1,dat)
# ------ teme of date
prec,psia,wa,ea,xa = obu.precess(ttt,'80')
deltapsi,trueeps,meaneps,omega,nut = obu.nutation(ttt,ddpsi,ddeps)

eqeg = deltapsi * np.cos(meaneps)
eqeg = np.fmod(eqeg,2.0 * np.pi)
eqe[0,0] = np.cos(eqeg)
eqe[0,2] = np.sin(eqeg)
eqe[0,2] = 0.0
eqe[1,0] = - np.sin(eqeg)
eqe[1,1] = np.cos(eqeg)
eqe[1,2] = 0.0
eqe[2,0] = 0.0
eqe[2,1] = 0.0
eqe[2,2] = 1.0
tm = prec * nut * np.transpose(eqe)
print(f'of date transformation teme - eci   \n{tm} \n')
# -------- teme2ecef    - transform teme to ecef vectors
recefd,vecefd,aecefd = sc.teme2ecef(rteme,vteme,ateme,ttt,jdut1 + jdut1frac,lod,xp,yp,eqeterms)
print('teme-ecef of date\n recefd \n', (recefd))
print(' vecefd \n', (vecefd))
# -------- teme2eci    - transform teme to eci vectors
print('\n\n start from teme \n')
recid,vecid,aecid = sc.teme2eci(rteme,vteme,ateme,ttt,ddpsi,ddeps)
print('teme-eci of date\n recid \n', (recid))
print(' vecid \n' % (vecid))
# convert ecef to eci standard
# -------- ecef2eci    - transform ecef to eci vectors
reci,veci,aeci = sc.ecef2eci(recefd,vecefd,aecefd,ttt,jdut1 + jdut1frac,lod,xp,yp,2,ddpsi,ddeps)
print(f'ecef-eci fk5/iau76\n reci {reci}')
print(f' veci {veci}\n')
# ------ teme of epoch
# -------- teme2eci    - transform teme to eci vectors using epoch tm matrix
print(f'of epoch transformation teme - eci  \n{tme} \n')
recie = tme * rteme
vecie = tme * vteme
aecie = tme * ateme
print('teme-eci of epoch\n recie \n', (recie))
print(' vecie \n', (vecie))
# ---- now convert back to ecef using std techniques
recefe,vecefe,aecefe = sc.eci2ecef(recie,vecie,aecie,ttt,jdut1 + jdut1frac,lod,xp,yp,eqeterms,ddpsi,ddeps)
print('teme-ecef of epoch\n recefe \n',(recefe))
print(' vecefe \n', (vecefe))
dr = recid - recie
dv = vecid - vecie
print('eci diffs\n dr (m) \n', (dr * 1000.0))
print(' dv \n', (dv * 1000.0))
#print('%11.7f \n' % (smu.mag(dr)))
dr = recefd - recefe
dv = vecefd - vecefe
print('ecef diffs\n dr (m) \n', (dr * 1000.0))
print(' dv \n', (dv * 1000.0))
#print('%11.7f \n' % (smu.mag(dr)))
