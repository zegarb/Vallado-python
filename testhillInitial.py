#
# test the hills conversion programs
# this does the intiial checkouts
#  use testhillSTK for te longer runs with perturbed values, etc
#
# Uses:
#  EQCM_to_ECI_RTN_
#      f_ECI_to_RTN_
#      newtonnu, newtone
#      inverselliptic2
#      elliptic12
#  ECI_to_EQCM_RTN_
#
#
#
#
# STK needs to be running (testhillsc)
#
#
#

import numpy as np
from space_constants import *
import os
import spacemath_utils as smu
import space_conversions as sc
print('------------------------------------------------ initial accuracy ----------------------------- \n')
# now test the ability to convert eci - hills and back
casenumo = 1

casetest = 7

ang_step = 1e-08
# case where you read in the int and tgt ephemerides
ropt = 'm'
outfilehill = open(os.path.join(os.path.dirname(__file__), 'testoutput',
                   f'thillc{casenumo}n{casetest}.out'), 'wt')
# outfilehill = open(strcat('d:/STKFiles Educational Files/Hills/thillc', int2str(casenumo),'n', int2str(casetest),'.out'),'wt')

# -----------------------------------------------------------------------------------------------
# -------------------------------- do initial accuracy checks fwd, back, ------------------------
# -----------------------------------------------------------------------------------------------
# --- this one also sets the case for the rest of the program here...
for ktr in range(1, 12):
    # set last one to be the case of interest
    if ktr < 11:
        casenum = ktr
    else:
        casenum = casenumo
    rtgteci = np.zeros(3)
    vtgteci = np.zeros(3)
    match casenum:
        case 1:
            rtgteci[0] = 6378.137 + 500.0
            rtgteci[1] = 0.0
            rtgteci[2] = 0.0
            magrt = smu.mag(rtgteci)
            vtgteci[0] = 0.0
            vtgteci[1] = np.sqrt(mu / magrt)
            vtgteci[2] = 0.0
            circ = 'y'
            rtgteci = rtgteci.T
            vtgteci = vtgteci.T
            #  previous doesn't work in STK92 - fixed in 9.3 for smaller inclination
            a = 6378.137 + 500.0
            ecc = 0.0
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, 0.0, 0.001 * deg2rad, 0.0, 0.0,
                                         0.0, 0.0, 0.0, 0.0)
            circ = 'y'
            print('rtgt = [%20.13f %20.13f %20.13f]; '
                  '\n vtgt = [%20.13f %20.13f %20.13f]; '
                  '\n' % (rtgteci[0], rtgteci[1], rtgteci[2],
                          vtgteci[0], vtgteci[1], vtgteci[2]))
        case 2:
            a = 26500.0
            ecc = 0.0
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, 0.0, 0.001 * deg2rad, 0.0, 0.0,
                                         0.0, 0.0, 0.0, 0.0)
            circ = 'y'
        case 3:
            a = 26500.0
            ecc = 0.0
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, 0.0, 55.0 * deg2rad, 0.0, 0.0, 0.0,
                                          0.0, 0.0, 0.0)
            circ = 'y'
        case 4:
            a = 42164.0
            ecc = 0.0
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, 0.0, 0.001 * deg2rad, 0.0, 0.0,
                                         0.0, 0.0, 0.0, 0.0)
            circ = 'y'
        case 5:
            rtgteci = np.array([-605.7904308,- 5870.230407, 3493.052004])
            vtgteci = np.array([- 1.568251615,- 3.702348353,- 6.479484915])
            circ = 'n'
        case 6:
            rtgteci = np.array([- 761.24075519, 22699.40899449, 13644.176032])
            vtgteci = np.array([- 2.630715586, 1.396719096,- 2.491891])
            circ = 'n'
        case 7:
            rtgteci = np.array([- 40588.150362,- 11462.167028, 27.147649])
            vtgteci = np.array([0.834787457,- 2.958305691,- 0.001173016])
            circ = 'n'
        case 8:
            rtgteci = np.array([- 5551.898646,- 2563.049696, 3257.756165])
            vtgteci = np.array([2.149073,- 7.539457,- 2.185709])
            circ = 'n'
        case 9:
            rtgteci = np.array([9668.14551571, 6999.07240705, 4041.43303674])
            vtgteci = np.array([- 3.65292306, 0.64966519,- 5.82123516])
            circ = 'n'
        case 10:
            a = 8164.7188
            ecc = 0.08
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, ecc, 32.861 * deg2rad, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            circ = 'n'
    x = 10.0
    y = 10.0
    z = 10.0
    xd = 0.01
    yd = 0.01
    zd = 0.01
    hro1 = np.zeros(3)
    hvo1 = np.zeros(3)
    hro1[0] = x / 1000.0
    hro1[1] = y / 1000.0
    hro1[2] = z / 1000.0
    hvo1[0] = xd / 1000.0
    hvo1[1] = yd / 1000.0
    hvo1[2] = zd / 1000.0
    rintecix, vintecix = sc.hilleqcm2eci(rtgteci, vtgteci, hro1, hvo1)
    rhillx, vhillx = sc.eci2hilleqcm(rtgteci, vtgteci, rintecix, vintecix)
    dr = rhillx - hro1
    dv = vhillx - hvo1
    mdr = smu.mag(dr * 1000000)
    mdv = smu.mag(dv * 1000000)
    print('dr m         %20.13f %20.13f %20.13f %15.8f mm %20.13f %20.13f %20.13f %15.8f mm/s \n'
          % (dr[0] * 1000, dr[1] * 1000, dr[2] * 1000, mdr, dv[0] * 1000,
             dv[1] * 1000, dv[2] * 1000, mdv))
    p, a, ecc, incl, omega, argpx1, nux1, m, arglat, truelon, lonper = sc.rv2coe(rtgteci, vtgteci)
    print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n'
          % (p, a, ecc, incl * rad2deg, omega * rad2deg, argpx1 * rad2deg,
             nux1 * rad2deg, m * rad2deg, arglat * rad2deg, truelon * rad2deg,
             lonper * rad2deg))
    p, a, ecc, incl, omega, argpx1, nux1, m, arglat, truelon, lonper = sc.rv2coe(rintecix, vintecix)
    print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n'
          % (p, a, ecc, incl * rad2deg, omega * rad2deg, argpx1 * rad2deg,
             nux1 * rad2deg, m * rad2deg, arglat * rad2deg, truelon * rad2deg,
             lonper * rad2deg))

# -----------------------------------------------------------------------------------------------
# ------------------- check various positions to determine if the veocity is correct ------------
# -----------------------------------------------------------------------------------------------
for ktr in range(1, 7):
    if ktr == 1:
        x = 100.0
    else:
        x = 0.0
    if ktr == 2:
        y = 100.0
    else:
        y = 0.0
    if ktr == 3:
        z = 100.0
    else:
        z = 0.0
    if ktr == 4:
        xd = 0.01
    else:
        xd = 0.0
    if ktr == 5:
        yd = 0.01
    else:
        yd = 0.0
    if ktr == 6:
        zd = 0.01
    else:
        zd = 0.0
    hro1[0] = x / 1000.0
    hro1[1] = y / 1000.0
    hro1[2] = z / 1000.0
    hvo1[0] = xd / 1000.0
    hvo1[1] = yd / 1000.0
    hvo1[2] = zd / 1000.0
    rinteci, vinteci = sc.hilleqcm2eci(rtgteci, vtgteci, hro1, hvo1)
    print('hillsin       %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
          % (hro1[0] * 1000, hro1[1] * 1000, hro1[2] * 1000,
             hvo1[0] * 1000, hvo1[1] * 1000, hvo1[2] * 1000))
    mmagrti = smu.mag(rtgteci)
    mmagvti = smu.mag(vtgteci)
    print('rtgteci        %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f %15.8f %15.8f \n'
          % (rtgteci[0], rtgteci[1], rtgteci[2], vtgteci[0], vtgteci[1],
             vtgteci[2], mmagrti, mmagvti))
    mmagri = smu.mag(rinteci)
    mmagvi = smu.mag(vinteci)
    print('rinteci        %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f  %15.8f %15.8f \n'
          % (rinteci[0], rinteci[1], rinteci[2], vinteci[0],
             vinteci[1], vinteci[2], mmagri, mmagvi))

# -----------------------------------------------------------------------------------------------
# -------------------------- perturb each one, or all of the vector components ------------------
# -----------------------------------------------------------------------------------------------
print('\n\n------------------------------- perturb each one ---------------------------------------- \n')
# initialize
x = 0.0

y = 0.0

z = 0.0

xd = 0.0

yd = 0.0

zd = 0.0

fid = 2
for ktr in range(1, 7):
    if ktr == 1:
        x = 100.0
    else:
        x = 0.0
    if ktr == 2:
        y = 100.0
    else:
        y = 0.0
    if ktr == 3:
        z = 100.0
    else:
        z = 0.0
    if ktr == 4:
        xd = 1
    else:
        xd = 0.0
    if ktr == 5:
        yd = 1
    else:
        yd = 0.0
    if ktr == 6:
        zd = 1
    else:
        zd = 0.0
    hro1[0] = x / 1000.0
    hro1[1] = y / 1000.0
    hro1[2] = z / 1000.0
    hvo1[0] = xd / 1000.0
    hvo1[1] = yd / 1000.0
    hvo1[2] = zd / 1000.0

    # reset this!! so m for hills call!!
    hro = np.zeros(3)
    hvo = np.zeros(3)
    hro[0] = x
    hro[1] = y
    hro[2] = z
    hvo[0] = xd
    hvo[1] = yd
    hvo[2] = zd
    rinteci, vinteci = sc.hilleqcm2eci(rtgteci, vtgteci, hro1, hvo1)
    rhill2, vhill2 = sc.eci2hilleqcm(rtgteci, vtgteci, rinteci, vinteci)
    rinteci1, vinteci1 = sc.hilleqcm2eci(rtgteci, vtgteci, rhill2, vhill2)
    rhill, vhill = sc.eci2hilleqcm(rtgteci, vtgteci, rinteci1, vinteci1)
    rhill1, vhill1 = sc.eci2hilleqcm(rtgteci, vtgteci, rinteci, vinteci)

    # --------------------- write out various transformation methods --- just at ktr = 8
    if (ktr <= 8):
        print('\nhills in ell   %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (hro[0], hro[1], hro[2], hvo[0], hvo[1], hvo[2]))
        print('\nhills in       %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (hro1[0] * 1000, hro1[1] * 1000, hro1[2] * 1000,
                 hvo1[0] * 1000, hvo1[1] * 1000, hvo1[2] * 1000))
        mmagrti = smu.mag(rtgteci)
        mmagvti = smu.mag(vtgteci)
        print('rtgteci        %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n'
              % (rtgteci[0], rtgteci[1], rtgteci[2], mmagrti, vtgteci[0],
                 vtgteci[1], vtgteci[2], mmagvti))
        print('rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (rhill2[0] * 1000, rhill2[1] * 1000, rhill2[2] * 1000,
                 vhill2[0] * 1000, vhill2[1] * 1000, vhill2[2] * 1000))
        mmagrti = smu.mag(rinteci)
        mmagvti = smu.mag(vinteci)
        print('rinteci        %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n'
              % (rinteci[0], rinteci[1], rinteci[2], mmagrti, vinteci[0],
                 vinteci[1], vinteci[2], mmagvti))
        mmagrti = smu.mag(rinteci1)
        mmagvti = smu.mag(vinteci1)
        print('rinteci1 eq %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n'
              % (rinteci1[0], rinteci1[1], rinteci1[2], mmagrti,
                 vinteci1[0], vinteci1[1], vinteci1[2], mmagvti))
        print('hills hill  %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (rhill[0] * 1000, rhill[1] * 1000, rhill[2] * 1000,
                 vhill[0] * 1000, vhill[1] * 1000, vhill[2] * 1000))
        print('hills1 eq   %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (rhill1[0] * 1000, rhill1[1] * 1000, rhill1[2] * 1000,
                 vhill1[0] * 1000, vhill1[1] * 1000, vhill1[2] * 1000))
        print('==== s codes ====\n')
        #---            [rinteci, vinteci]   = EQCM_to_ECI_NTW_( rtgteci*1000, vtgteci*1000, hro', hvo' );      # do in m~!!!find int pos eqcm
        #---            [rhill2, vhill2]          = ECI_to_EQCM_NTW_( rtgteci*1000, vtgteci*1000, rinteci, vinteci ); # find how much this would be for true hills
        rinteci, vinteci = sc.hilleqcm2eci(rtgteci, vtgteci, hro, hvo)
        rhill2, vhill2 = sc.eci2hilleqcm(rtgteci, vtgteci, rinteci, vinteci)
        mmagrti = smu.mag(rinteci)
        mmagvti = smu.mag(vinteci)
        print('rinteci hil %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n'
              % (rinteci[0] * 0.001, rinteci[1] * 0.001,
                 rinteci[2] * 0.001, mmagrti * 0.001, vinteci[0] * 0.001,
                 vinteci[1] * 0.001, vinteci[2] * 0.001, mmagvti * 0.001))
        print('rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (rhill2[0], rhill2[1], rhill2[2], vhill2[0], vhill2[1], vhill2[2]))

outfilehill.close()
