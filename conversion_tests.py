from space_conversions import *
import numpy as np
from pprint import pprint as pp

latgc = math.pi*0.20
latgd = gc2gd(latgc)
print(latgd, latgc)
latgc = gd2gc(latgd)
print(latgd, latgc)

apin = 250
print(apin)
kpout = ap2kp(apin)
print(kpout)
apout = kp2ap(kpout)
print(apout)

# LEO test
recef = np.array([[-1033.4793830],  [7901.2952754],  [6380.3565958]])
vecef = np.array([[-3.225636520],  [-2.872451450],   [5.531924446]])
aecef = np.array([[0.001],[0.002],[0.003]])

year=2004
mon = 4
day = 6
hr =  7
min= 51
sec= 28.386009

dut1 = -0.4399619  # sec
dat  = 32         # sec
xp   = -0.140682 * arcsec2rad  # " to rad
yp   =  0.333309 * arcsec2rad
lod  =  0.0015563
ddpsi = -0.052195 * arcsec2rad  # " to rad
ddeps = -0.003875 * arcsec2rad
ddx = -0.000205 * arcsec2rad  # " to rad
ddy = -0.000136 * arcsec2rad
order = 106
terms = 2
eqeterms = 2
timezone=0
opt = 'c' # specify the iau00 cio approach

print('input data \n\n')
print(' year %5i '% year)
print(' mon %4i '% mon)
print(' day %3i '% day)
print(' %3i:%2i:%8.6f\n'% (hr, min, sec))
print(' dut1 %8.6f s'% dut1)
print(' dat %3i s'% dat)
print(' xp %8.6f "'% (xp * rad2arcsec))
print(' yp %8.6f "'% (yp * rad2arcsec))
print(' lod %8.6f s\n'% lod)
print(' ddpsi %8.6f " ddeps  %8.6f\n'% (ddpsi * rad2arcsec, ddeps * rad2arcsec))
print(' ddx   %8.6f " ddy    %8.6f\n'% (ddx * rad2arcsec, ddy * rad2arcsec))
print(' order %3i  eqeterms %3i  opt %3s \n'% (order, eqeterms, opt))
print('units are km and km/s and km/s2\n')

# -------- convtime    - convert time from utc to all the others
print('convtime results\n')
ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, \
    ttdb, jdtdb, jdtdbfrac \
    = stu.convtime(year, mon, day, hr, min, sec, timezone, dut1, dat)
print('ut1 %8.6f tut1 %16.12f jdut1 %18.11f '% (ut1, tut1, jdut1+jdut1frac))

print('hms2rad, rad2hms')
hms = hms2rad(hr, min, sec)
print(hms)
hr, min, sec = rad2hms(hms)
print('%f hr, %f min, %f sec\n' % (hr, min, sec))


print('input vectors:')
pp(recef)
pp(vecef)
pp(aecef)


#jdut1 2453101.82740678312


nu = lon2nu (jdut1, 7.020438698 * deg2rad, 0.070273056 * deg2rad, 19.90450011 * deg2rad,
                352.5056022 * deg2rad)
print("nu is ", nu)
lon = nu2lon(jdut1,nu,0.070273056 * deg2rad, 19.90450011 * deg2rad, 352.5056022 * deg2rad)
print("lon is ", lon)

#--------------------------------ecef2--------------------------------------------------------------
rpef, vpef, apef = ecef2pef(recef, vecef, aecef, '80', xp, yp, ttt)
print('ecef2pef 80 returned: ')
pp(rpef)
pp(vpef)
pp(apef)

recef, vecef, aecef = pef2ecef(rpef, vpef, apef, '80', xp, yp, ttt)
print('pef2ecef 80 returned: ')
pp(recef)
pp(vecef)
pp(aecef)

# technically tirs conversions and not pef for 2000 theory.
# don't know if it should be given its own function name -zeg
rpef, vpef, apef = ecef2pef(recef, vecef, aecef, '01', xp, yp, ttt)
print('ecef2pef 01 returned: ')
pp(rpef)
pp(vpef)
pp(apef)

recef, vecef, aecef = pef2ecef(rpef, vpef, apef, '01', xp, yp, ttt)
print('pef2ecef 01 returned: ')
pp(recef)
pp(vecef)
pp(aecef)

# aI, aJ ecef values slightly off
rtod, vtod, atod = ecef2tod(recef, vecef, aecef, ttt,
                            jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps)
print('ecef2tod returned: ')
pp(rtod)
pp(vtod)
pp(atod)
recef, vecef, aecef = tod2ecef(rtod, vtod, atod, ttt, jdut1 + jdut1frac, lod,
                               xp, yp, 2, ddpsi, ddeps)
print('tod2ecef returned: ')
pp(recef)
pp(vecef)
pp(aecef)

# aI, aJ ecef values slightly off
rmod, vmod, amod = ecef2mod(recef, vecef, aecef, ttt,
                            jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps)
print('ecef2mod returned: ')
pp(rmod)
pp(vmod)
pp(amod)
recef, vecef, aecef = mod2ecef(rmod, vmod, amod, ttt, jdut1 + jdut1frac, lod,
                               xp, yp, 2, ddpsi, ddeps)
print('mod2ecef returned: ')
pp(recef)
pp(vecef)
pp(aecef)


recig, vecig, aecig = ecef2eci(recef, vecef, aecef, ttt,
                                jdut1+jdut1frac, lod, xp, yp, 2,
                                ddpsi, ddeps)
print('ecef2eci returned: ')
pp(recig)
pp(vecig)
pp(aecig)
recef, vecef, aecef = eci2ecef(recig, vecig, aecig, ttt,
                                jdut1+jdut1frac, lod, xp, yp, 2,
                                ddpsi, ddeps)
print('eci2ecef returned:')
pp(recef)
pp(vecef)
pp(aecef)

rcirs,vcirs,acirs = ecef2cirsiau06(recef,vecef,aecef,ttt,jdut1,lod,xp,yp,'a',
                                   ddx, ddy)
print('ecef2cirsiau06 a returned: ')
pp(rcirs)
pp(vcirs)
pp(acirs)

recef,vecef,aecef = cirs2ecefiau06(rcirs,vcirs,acirs,ttt,jdut1,lod,xp,yp,'a',
                                   ddx, ddy)
print('cirs2ecefiau06 a returned: ')
pp(recef)
pp(vecef)
pp(aecef)

rcirs,vcirs,acirs = ecef2cirsiau06(recef,vecef,aecef,ttt,jdut1,lod,xp,yp,'b',
                                   ddx, ddy)
print('ecef2cirsiau06 b returned: ')
pp(rcirs)
pp(vcirs)
pp(acirs)

recef,vecef,aecef = cirs2ecefiau06(rcirs,vcirs,acirs,ttt,jdut1,lod,xp,yp,'b',
                                   ddx, ddy)
print('cirs2ecefiau06 b returned: ')
pp(recef)
pp(vecef)
pp(aecef)

rcirs,vcirs,acirs = ecef2cirsiau06(recef,vecef,aecef,ttt,jdut1,lod,xp,yp,'c',
                                   ddx, ddy)
print('ecef2cirsiau06 c returned: ')
pp(rcirs)
pp(vcirs)
pp(acirs)

recef,vecef,aecef = cirs2ecefiau06(rcirs,vcirs,acirs,ttt,jdut1,lod,xp,yp,'c',
                                   ddx, ddy)
print('cirs2ecefiau06 c returned: ')
pp(recef)
pp(vecef)
pp(aecef)

reci,veci,aeci = ecef2eciiau06(recef,vecef,aecef,ttt,jdut1,
                                lod,xp,yp,'a', ddx, ddy)
print('ecef2eciiau06 a returned: ')
pp(reci)
pp(veci)
pp(aeci)

recef,vecef,aecef = eci2ecefiau06(reci,veci,aeci,ttt,jdut1,
                                    lod,xp,yp,'a', ddx, ddy)
print('eci2ecefiau06 a returned: ')
pp(recef)
pp(vecef)
pp(aecef)

reci,veci,aeci = ecef2eciiau06(recef,vecef,aecef,ttt,jdut1,
                                lod,xp,yp,'b', ddx, ddy)
print('ecef2eciiau06 b returned: ')
pp(reci)
pp(veci)
pp(aeci)

recef,vecef,aecef = eci2ecefiau06(reci,veci,aeci,ttt,jdut1,
                                    lod,xp,yp,'b', ddx, ddy)
print('eci2ecefiau06 b returned: ')
pp(recef)
pp(vecef)
pp(aecef)

reci,veci,aeci = ecef2eciiau06(recef,vecef,aecef,ttt,jdut1,
                                lod,xp,yp,'c', ddx, ddy)
print('ecef2eciiau06 c returned: ')
pp(reci)
pp(veci)
pp(aeci)

recef,vecef,aecef = eci2ecefiau06(reci,veci,aeci,ttt,jdut1,
                                    lod,xp,yp,'c', ddx, ddy)
print('eci2ecefiau06 c returned: ')
pp(recef)
pp(vecef)
pp(aecef)

rteme, vteme, ateme = ecef2teme(recef, vecef, aecef, ttt,
                                jdut1+jdut1frac, lod, xp, yp, eqeterms)
print('ecef2teme returned: ')
pp(rteme)
pp(vteme)
pp(ateme)

recef, vecef, aecef = teme2ecef(rteme, vteme, ateme, ttt,
                                jdut1+jdut1frac, lod, xp, yp, eqeterms)
print('teme2ecef returned: ')
pp(recef)
pp(vecef)
pp(aecef)

#--------------------------------eci2--------------------------------------------------------------
rteme, vteme, ateme = eci2teme(recig, vecig, aecig, ttt, ddpsi, ddeps)
print('eci2teme returned:')
pp(rteme)
pp(vteme)
pp(ateme)

recig, vecig, aecig, = teme2eci(rteme, vteme, ateme,
                                tt, ddpsi, ddeps)
print('teme2eci returned:')
pp(recig)
pp(vecig)
pp(aecig)

rtod, vtod, atod = eci2tod(recig, vecig, aecig, ttt, ddpsi, ddeps)
print('eci2tod returned:')
pp(rtod)
pp(vtod)
pp(atod)

recig, vecig, aecig = tod2eci(rtod, vtod, atod, ttt, ddpsi, ddeps)
print('tod2eci returned:')
pp(recig)
pp(vecig)
pp(aecig)

rmod, vmod, amod = eci2mod(recig, vecig, aecig, ttt)
print('eci2mod returned:')
pp(rmod)
pp(vmod)
pp(amod)

recig, vecig, aecig = mod2eci(rmod, vmod, amod, ttt)
print('mod2eci returned:')
pp(recig)
pp(vecig)
pp(aecig)

rpef, vpef, apef = eci2pef(recig, vecig, aecig, ttt, jdut1, lod, terms, ddpsi,
                           ddeps)
print('eci2pef returned: ')
pp(rpef)
pp(vpef)
pp(apef)

recig, vecig, aecig = pef2eci(rpef, vpef, apef, ttt, jdut1, lod, eqeterms,
                              ddpsi, ddeps)
print('pef2eci returned: ')
pp(recig)
pp(vecig)
pp(aecig)

rtirs, vtirs, atirs = eci2tirsiau06(recig, vecig, aecig, "a", ttt, jdut1, lod,
                                    0, 0)
print("eci2tirsiau06 a returned:")
pp(rtirs)
pp(vtirs)
pp(atirs)

recig, vecig, aecig = tirs2eciiau06(rtirs, vtirs, atirs, "a", ttt, jdut1, lod,
                                 0, 0)
print("tirs2eciiau06 a returned:")
pp(recig)
pp(vecig)
pp(aecig)

rtirs, vtirs, atirs = eci2tirsiau06(recig, vecig, aecig, "b", ttt, jdut1, lod,
                                    0, 0)
print("eci2tirsiau06 b returned:")
pp(rtirs)
pp(vtirs)
pp(atirs)

recig, vecig, aecig = tirs2eciiau06(rtirs, vtirs, atirs, "b", ttt, jdut1, lod,
                                 0, 0)
print("tirs2eciiau06 b returned:")
pp(recig)
pp(vecig)
pp(aecig)

rtirs, vtirs, atirs = eci2tirsiau06(recig, vecig, aecig, "c", ttt, jdut1, lod,
                                    0, 0)
print("eci2tirsiau06 c returned:")
pp(rtirs)
pp(vtirs)
pp(atirs)

recig, vecig, aecig = tirs2eciiau06(rtirs, vtirs, atirs, "c", ttt, jdut1, lod,
                                 0, 0)
print("tirs2eciiau06 c returned:")
pp(recig)
pp(vecig)
pp(aecig)

rcirs,vcirs,acirs = eci2cirsiau06(recig,vecig,aecig,ttt,'a', ddx, ddy)
print('eci2cirsiau06 a returned: ')
pp(rcirs)
pp(vcirs)
pp(acirs)

recig,vecig,aecig = cirs2eciiau06(rcirs,vcirs,acirs,ttt,'a', ddx, ddy)
print('cirs2eciiau06 a returned: ')
pp(recig)
pp(vecig)
pp(aecig)

rcirs,vcirs,acirs = eci2cirsiau06(recig,vecig,aecig,ttt,'b', ddx, ddy)
print('eci2cirsiau06 b returned: ')
pp(rcirs)
pp(vcirs)
pp(acirs)

recig,vecig,aecig = cirs2eciiau06(rcirs,vcirs,acirs,ttt,'b', ddx, ddy)
print('cirs2eciiau06 b returned: ')
pp(recig)
pp(vecig)
pp(aecig)

rcirs,vcirs,acirs = eci2cirsiau06(recig,vecig,aecig,ttt,'c', ddx, ddy)
print('eci2cirsiau06 c returned: ')
pp(rcirs)
pp(vcirs)
pp(acirs)

recig,vecig,aecig = cirs2eciiau06(rcirs,vcirs,acirs,ttt,'c', ddx, ddy)
print('cirs2eciiau06 c returned: ')
pp(recig)
pp(vecig)
pp(aecig)

print("----------------------------------------")

a = 6860.7631
ecc = 0.75
p = a*(1.0 - ecc**2)
incl = 15 * deg2rad #97.65184/rad
omega = 79.54701 * deg2rad
argp = 83.86041 * deg2rad
nu = 65.21303 * deg2rad
arglat = 0.0
truelon = 0.0
lonper = 0.0

print('p km=%f  a km=%f  ecc=%f  incl deg=%f  raan deg=%f argp deg=%f  nu deg=%f'% \
        (p,a,ecc,incl * rad2deg,omega * rad2deg,argp * rad2deg,nu * rad2deg))
print('incl rad=%f  raan rad=%f  argp rad=%f  nu rad=%f'% \
        (incl, omega, argp, nu))
reci,veci = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper)
print("coe2rv returned: ", reci, veci)

#    reci = np.array([[11074.95274], [40629.74421], [-32.1123199]])
#    veci = np.array([[-2.940822436], [0.9007122363], [0.002036330819]])
#    reci = np.array([11074.95274, 40629.74421, -32.1123199])
#    veci = np.array([-2.940822436, 0.9007122363, 0.002036330819])
#    reci = np.array([6524.834000000,  6862.875000000, 6448.296000000])
#    veci = np.array([4.9013270000,    5.5337560000,   -1.9763410000])
p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = \
    rv2coe(reci, veci)
print('p km=%f  a km=%f  ecc=%f  incl =%f  raan=%f  argp=%f  nu=%f  m=%f '% \
        (p,a,ecc,incl,omega,argp,nu,m))
print('p km=%f  a km=%f  ecc=%f  incl deg=%f  raan deg=%f  argp deg=%f  nu deg=%f  m deg=%f '% \
        (p,a,ecc,incl * rad2deg,omega*rad2deg,argp*rad2deg,nu*rad2deg,m*rad2deg))
print('     arglat   truelon    lonper ',\
        arglat,truelon,lonper)
#            arglat * rad2deg,truelon * rad2deg,lonper * rad2deg)

#---------------------------rv2-----------------------------------------------

#alt, terms needs init
alt = 0.0
terms = 2
rho, trtasc, tdecl, drho, dtrtasc, dtdecl = rv2tradec(reci, veci, latgd, lon,
                                                        alt, ttt, jdut1, lod,
                                                        xp, yp, terms, ddpsi,
                                                        ddeps)
print("rv2tradec")
print('rho %f trtasc %f tdecl %f drho %f dtrtasc %f dtdecl %f' %
        (rho, trtasc, tdecl, drho, dtrtasc, dtdecl))

#   tradc2rv()

rrsw, vrsw, transmat = rv2rsw(reci, veci)
print('rv2rsw:')
print(rrsw)
print(vrsw)
print(transmat)

#rsw2rv()

rntw, vntw, transmat = rv2ntw(reci, veci)
print('rv2ntw:')
print(rntw)
print(vntw)
print(transmat)

#ntw2rv()

lon,latgc,rtasc,decl,fpa,az,magr,magv = rv2flt(reci, veci, ttt, jdut1, lod,
                                                xp, yp, terms, ddpsi, ddeps,)
print('rv2flt')
print('lon %f, latgc %f, rtasc %f, decl %f, fpa %f, az %f, magr %f,\
        magv %f' % (lon, latgc, rtasc, decl, fpa, az, magr, magv))

reci,veci = flt2rv(magr, magv, latgc, lon, fpa, az, ttt, jdut1, lod, xp,
                    yp, terms, ddpsi, ddeps)

print('flt2rv:')
print(reci)
print(veci)

a,n,af,ag,chi,psi,meanlonM,meanlonNu,fr = rv2eq(reci, veci)
print('rv2eq:')
print('a %f, n %f, af %f, ag %f, chi %f, psi %f, meanlonM %f, meanlonNu\
        %f, fr %f' % (a,n,af,ag,chi,psi,meanlonM,meanlonNu,fr))

#r = reci, v = veci?
r, v = eq2rv(a, af, ag, chi, psi, meanlonM, fr)
print('eq2rv:')
print(r)
print(v)

rmag,vmag,rtasc,decl,fpav,az = rv2adbar(r, v)
print('rv2adbar:')
print('rmag %f, vmag %f, rtasc %f, decl %f, fpav %f, az %f' % \
        (rmag, vmag, rtasc, decl, fpav, az))

r, v = adbar2rv(rmag, vmag, rtasc, decl, fpav, az)
print('adbar2rv:')
print(r)
print(v)

# rv2coeS was deleted and replaced by rv2coe
# p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper = rv2coeS(r, v)
# print('rv2coeS:')
# print(f'p {p}, a {a}, ecc {ecc}, incl {incl}, omega {omega}, argp {argp}, nu {nu}, arglat {arglat},\
#         truelon {truelon}, lonper {lonper}') # % (p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper))

# rv2coeh was combined with rv2coe
# p,a,ecc,incl,raan,argp,nu,m,arglat,truelon,lonper = rv2coeh(r, v, re, mu)
# print('rv2coeh:')
# print(f'p {p}, a {a}, ecc {ecc}, incl {incl}, omega {omega}, argp {argp}, nu {nu}, arglat {arglat},\
#         truelon {truelon}, lonper {lonper}') # % (p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper))


p,a,ecc,incl,raan,argp,nu,m,arglat,truelon,lonper = rv2coe(r, v)
print('rv2coe:')
print(f'p {p}, a {a}, ecc {ecc}, incl {incl}, omega {omega}, argp {argp}, nu {nu}, arglat {arglat},\
truelon {truelon}, lonper {lonper}') # % (p,a,ecc,incl,omega,argp,nu,m,arglat,truelon,lonper))

# coe2rvh was combined with coe2rv
r, v = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper, mu)
print('coe2rv:')
print(r)
print(v)

rr,rtasc,decl,drr,drtasc,ddecl = rv2radec(r, v)
print('rv2radec:')
print ('rr %f, rtasc %f, decl %f, drr %f, drtasc %f, ddecl %f' \
        % (rr,rtasc,decl,drr,drtasc,ddecl))

r, v = radec2rv(rr, rtasc, decl, drr, drtasc, ddecl)
print('radec2rv')
print(r)
print(v)

rho, trtasc, tdecl, drho, dtrtasc, dtdecl = rv2tradec(reci, veci, latgd,
                                                        lon, alt, ttt, jdut1,
                                                        lod, xp, yp, terms,
                                                        ddpsi, ddeps)
print('rv2tradec:')
print('rho %f, trtasc %f, tdecl %f, drho %f, dtrtasc %f, dtdecl %f' %\
        (rho, trtasc, tdecl, drho, dtrtasc, dtdecl))

#rseci, vseci needs init
#reci, veci = tradec2rv(rho, trtasc, tdecl, drho, dtrtasc, dtdecl, rseci,
#                       vseci, lod)
#print('tradec2rv:')
#print(reci)
#print(veci)

rho, az, el, drho, daz, del_ = rv2razel(reci, veci, latgd, lon, alt, ttt,
                                        jdut1, lod, xp, yp, terms, ddpsi,
                                        ddeps)
print('rv2razel:')
print(f'rho {rho}, az {az}, drho {drho}, daz {daz}, del {del_}')

reci, veci = razel2rv(rho, az, el, drho, daz, del_, latgd, lon, alt, ttt,
                        jdut1, lod, xp, yp, terms, ddpsi, ddeps)
print('razel2rv:')
print(reci)
print(veci)

rr, ecllon, ecllat, drr, decllon, decllat = rv2ell(reci, veci)
print('rv2ell')
print('rr %f, ecllon %f, ecllat %f, drr %f, decllon %f, decllat %f' %\
        (rr, ecllon, ecllat, drr, decllon, decllat))

rijk, vijk = ell2rv(rr, ecllon, ecllat, drr, decllon, decllat)
print('ell2rv:')
print(rijk)
print(vijk)

rhosez, drhosez = raz2sez(rho, az, el, drho, daz, del_)
print('raz2sez:')
print(rhosez)
print(drhosez)

rho, az, el, drho, daz, del_ = sez2raz(rhosez, drhosez)
print('sez2raz:')
print('rho %f, az %f, el %f, drho %f, daz %f, del_ %f' %\
        (rho, az, el, drho, daz, del_))


#---- hilleqcm tests

print('hill to ecqm and ecqm to hill conversion test\n')
x = 500
y = 10
z = 250
dx = 200
dy = 0.01
dz = 20
rinteqcm = np.array([x, y, z])
vinteqcm = np.array([dx, dy, dz])
rtgteci = np.array([-605.7904308, -5870.230407, 3493.052004])
vtgteci = np.array([-1.568251615, -3.702348353, -6.479484915])

rintecix, vintecix = hilleqcm2eci(rtgteci, vtgteci, rinteqcm, vinteqcm)
pp(rintecix)
pp(vintecix)
rinteqcm, vinteqcm = eci2hilleqcm(rtgteci, vtgteci, rintecix, vintecix)
pp(rinteqcm)
pp(vinteqcm)


