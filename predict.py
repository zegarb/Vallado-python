import numpy as np
from space_constants import *
import space_conversions as sc
import spacetime_utils as stu
import orbit_utils as obu
import spacemath_utils as smu

#
#
#   prediction problem
#
#

testnum = 3
if testnum == 3:
    a = 6768.3568
    #p = 8634.2349
    #ecc = 0.002
    ecc = 0.0005770
    p = a * (1 - ecc**2)
    #incl = 50.0 * deg2rad
    incl = 51.6190 *deg2rad
    #omega = 115.0 * deg2rad
    omega = 13.334 *deg2rad
    #argp = 90.0 * deg2rad
    argp = 102.5680 * deg2rad
    #nu = 0.0 * deg2rad
    m = 257.5950 * deg2rad
    e0, nu = smu.newtonm(ecc,m)
    arglat = 0.0 * deg2rad
    truelon = 0.0 * deg2rad
    lonper = 0.0 * deg2rad
    reci,veci = sc.coe2rv(p,ecc,incl,omega,argp,nu,arglat,truelon,lonper)
    print(reci)
    print(veci)
    #reci = np.array([6585.038266,1568.184321,9.116355])
    #veci = np.array([- 1.1157766,4.6316816,6.0149576])
    #aeci = np.array([0.001,0.002,0.003])
    year = 1997
    mon = 4
    day = 1
    hr = 21 #UTC (EST = 16)
    min = 36
    sec = 0.0
    dut1 = - 0.2913774
    dat = 32
    xp = - 0.19108
    #year 1997.25 from IERS
    #xp = - 0.189116
    yp = 0.329624
    #year 1997.25 from IERS
    #yp = 0.331998
    lod = 0.0
    terms = 0
    timezone = 0

print('year %5i ' % (year))
print('mon %4i ' % (mon))
print('day %3i ' % (day))
print('hr %3i:%2i:%8.6f\n' % (hr,min,sec))
print('dut1 %8.6f s' % (dut1))
print(' dat %3i s' % (dat))
print(' xp %8.6f "' % (xp))
print(' yp %8.6f "' % (yp))
print(' lod %8.6f s\n' % (lod))
# -------- convtime    - convert time from utc to all the others
#            fprintf(1,'convtime results\n');
#            fprintf(1,'ut1 #8.6f tut1 #16.12f jdut1 #18.11f\n',ut1,tut1,jdut1 );
#            fprintf(1,'utc #8.6f\n',utc );
#            fprintf(1,'tai #8.6f\n',tai );
#            fprintf(1,'tt  #8.6f ttt  #16.12f jdtt  #18.11f\n',tt,ttt,jdtt );
#            fprintf(1,'tdb #8.6f ttdb #16.12f jdtdb #18.11f\n',tdb,ttdb,jdtdb );

# ---- perform prediction
latgd = 42.38 * deg2rad
lon = - 71.13 * deg2rad
alt = 0.024 #km (24m)
dtsec = 120.0
# timezone = 0 Assumming Epoch time (UTC)
jdepoch,jdepochf = stu.jday(year,mon,day,hr,min,sec)

def predict(reci, veci, latgd, lon, alt, dtsec, jdepoch, jdepochf, dut1, dat, xp, yp, lod, terms=0, timezone=0):
    ndot = 0.0
    nddot = 0.0
    rho = 0.0
    az = 0.0
    el = 0.0
    vis = 'radar sun'

    conv = math.pi / (180.0*3600.0)

    rsecef,vsecef = obu.site(latgd,lon,alt)
    print('site ecef :\n',rsecef,vsecef)
    print()

    year, mon, day, hr, min, sec = stu.invjday(jdepoch,jdepochf)

    for i in range(0,120):
        #                [reci1,veci1,error] =  kepler  ( reci,veci, i*dtsec );
    #                reci = reci';
    #                veci = veci';
        reci1,veci1 = obu.pkepler(reci,veci,i * dtsec,ndot,nddot)
        ut1,tut1,jdut1,jdut1frac,utc,tai,tt,ttt,jdtt,jdttfrac,tdb,ttdb,jdtdb,jdtdbfrac \
            = stu.convtime(year,mon,day,hr,min,sec + i * dtsec,timezone,dut1,dat)
         # -------------------- convert eci to ecef --------------------
        a = np.array([[0],[0],[0]])
        # What are ddpsi and ddeps supposed to be? ttt? lod?
        #if i == 106:
        #    reci1 = [-2811.27691,3486.2632,5069.5763]
        #    veci1 = [-6.859691, -2.964792, -1.764721]

        recef,vecef,aecef = sc.eci2ecef(reci1,veci1,a,ttt,jdut1 + jdut1frac,lod,xp,yp,terms,0.000014*conv,-0.000181*conv)
        print(f'Julian Time: {jdut1 + jdut1frac}')
        print(f'reci1 {i} x {reci1} {veci1}')
        print(f'recef {i} x {recef} {vecef} \n')
        # ------- find ecef range vector from site to satellite -------
        rhoecef = recef - rsecef
          # ------------- convert to sez for calculations ---------------
        tempvec = smu.rot3(rhoecef,lon)
        rhosez = smu.rot2(tempvec,halfpi - latgd)

        if i == 106:
            print(jdut1 + jdut1frac)
            print(f'reci1 {i} x {reci1} {veci1}')
            y,m,d,h,mn,s = stu.invjday(jdut1,jdut1frac - dut1 / 86400.0)
            print(f'{y} {m} {d:.0f} {h:.0f}:{mn:02.0f} {s} \n')

        if i == 106:
            print(f'recef {i} x {recef} {vecef} \n')


        if i == 106:
            print(f'rhosez {i} x {rhosez} \n')

        rho,az,el,drho,daz,del_ = sc.rv2razel(reci1,veci1,latgd,lon,alt,ttt,jdut1 + jdut1frac,lod,xp,yp,terms,0.0,0.0)
        # fprintf(1,'rvraz #14.7f#14.7f#14.7f#14.7f#14.7f#14.7f\n',rho,az * rad2deg,el * rad2deg,drho,daz * rad2deg,del * rad2deg );
        if az < 0.0:
            az = az + twopi
        if rhosez[2] > 0.0:
            rsun,rtasc,decl = obu.sun(jdtt + jdttfrac)
            if i == 106:
                print(f'rsun{i} {rsun} \n')
                print(f'rsun{i} {rsun*au} \n')
            rsun = rsun * au
            rseci,vseci,aeci = sc.ecef2eci(rsecef,vsecef,a,ttt,jdut1 + jdut1frac,lod,xp,yp,2, 0.0, 0.0)
            if i == 106:
                print(f'rseci {i} x {rseci} {vseci} \n')
            if np.dot(rsun,rseci) > 0.0:
                vis = 'radar sun'
            else:
                rxr = np.cross(rsun,reci1)
                magrxr = smu.mag(rxr)
                magr = smu.mag(reci1)
                magrsun = smu.mag(rsun)
                zet = np.arcsin(magrxr / (magrsun * magr))
                dist = smu.mag(reci1) * np.cos(zet - halfpi)
                if i == 106:
                    print('zet  %11.7f dist %11.7f  \n' % (zet * rad2deg,dist))
                if dist > re:
                    vis = 'visible'
                else:
                    vis = 'radar night'
        else:
            vis = 'not visible'
        y,m,d,h,mn,s = stu.invjday(jdut1,jdut1frac - dut1 / 86400.0)
        print('%5i %3i %3i %2i:%2i %6.3f %12s %11.7f  %11.7f  %11.7f  \n' % (y,m,d,h,mn,s,vis,rho,az * rad2deg,el * rad2deg))


    return jdut1, jdut1frac, rho, az, el, vis

predict(reci, veci, latgd, lon, alt, dtsec, jdepoch, jdepochf, dut1, dat, xp, yp, terms, timezone=0)
