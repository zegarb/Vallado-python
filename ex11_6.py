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
    #reci = np.array([6585.038266,1568.184321,9.116355])
    #veci = np.array([- 1.1157766,4.6316816,6.0149576])
    #aeci = np.array([0.001,0.002,0.003])
    # April 1, 1997, 4:36 pm ET
    year1 = 1997
    mon1 = 4
    day1 = 1
    hr1 = 21 #UTC (ET = 16)
    min1 = 36
    sec1 = 0.0
    dut1 = - 0.2913774
    # April 2, 1997 at 1:08:0.00 UTC
    year2 = year1
    mon2 = mon1
    day2 = 2
    hr2 = 1
    min2 = 8
    sec2 = 0.0
    dat = 32
    # Missing arcsec2rad conversion in matlab code
    xp = - 0.19108 * arcsec2rad
    #year1 1997.25 from IERS
    #xp = - 0.189116 * arcsec2rad
    yp = 0.329624 * arcsec2rad
    #year1 1997.25 from IERS
    #yp = 0.331998 * arcsec2rad
    lod = 0.0
    terms = 0
    timezone = 0

print('\nreci', end='')
print (reci)
print('veci', end='')
print(veci)
print('year1 %5i ' % (year1))
print('mon1 %4i ' % (mon1))
print('day1 %3i ' % (day1))
print('hr1 %3i:%2i:%8.6f UTC\n' % (hr1,min1,sec1))
print('dut1 %8.6f s' % (dut1))
print(' dat %3i s' % (dat))
print(' xp %8.6f "' % (xp*rad2arcsec))
print(' yp %8.6f "' % (yp*rad2arcsec))
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
jdepoch,jdepochf = stu.jday(year1,mon1,day1,hr1,min1,sec1)
jend, jendf = stu.jday(year2,mon2,day2,hr2,min2,sec2)

# 106 iterations gets the 1:08:0.00 UTC example (pg 913)
jdutend, rho, az, el, vis = obu.predict(reci, veci, jdepoch+jdepochf,jend+jendf,
                                        latgd, lon, alt, dtsec, dut1, dat, xp, yp)

y,m,d,h,mn,s = stu.invjday(jdutend - (dut1 / 86400.0))
print('%5i %3i %3i %2i:%2i %6.3f %12s %11.7f  %11.7f  %11.7f  \n' % \
              (y,m,d,h,mn,s,vis,rho,az * rad2deg,el * rad2deg))

jend, jendf = stu.jday(year2,mon2,day2,hr2,min2+26,sec2)


# more iterations to complete the table (Vallado 4th edition pg 914, Table 11-4)
jdutend, rho, az, el, vis = obu.predict(reci, veci, jdepoch+jdepochf, jend+jendf,
                                        latgd, lon, alt, dtsec, dut1, dat, xp, yp)
y,m,d,h,mn,s = stu.invjday(jdutend - (dut1 / 86400.0))
print('%5i %3i %3i %2i:%2i %6.3f %12s %11.7f  %11.7f  %11.7f  \n' % \
              (y,m,d,h,mn,s,vis,rho,az * rad2deg,el * rad2deg))

