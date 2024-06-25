#
# test the hills conversion programs
# this does the intiial checkouts
#  use testhillSTK for te longer runs with perturbed values, etc
#
# Uses:
#  EQCM_to_ECI_RTN_sal
#      f_ECI_to_RTN_sal
#      newtonnu, newtone
#      inverselliptic2
#      elliptic12
#  ECI_to_EQCM_RTN_sal
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
print('------------------------------------------------ initial accuracy ----------------------------- \n' % ())
# now test the ability to convert eci - hills and back
constmath
constastro
casenumo = 1

casetest = 7

ang_step = 1e-08
# case where you read in the int and tgt ephemerides
ropt = 'm'
outfilehill = open(os.path.join(os.path.dirname(__file__), 'testoutput',
                   'thillc' + casenumo + 'n' + casetest + '.out'), 'wt')
# outfilehill = open(strcat('d:/STKFiles Educational Files/Hills/thillc', int2str(casenumo),'n', int2str(casetest),'.out'),'wt')

# -----------------------------------------------------------------------------------------------
# -------------------------------- do initial accuracy checks fwd, back, ------------------------
# -----------------------------------------------------------------------------------------------
# --- this one also sets the case for the rest of the program here...
for ktr in range(1, 11):
    # set last one to be the case of interest
    if ktr < 10:
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
            rtgteci = np.array([-605.7904308,- 5870.230407, 3493.052004]).T
            vtgteci = np.array([- 1.568251615,- 3.702348353,- 6.479484915]).T
            circ = 'n'
        case 6:
            rtgteci = np.array([- 761.24075519, 22699.40899449, 13644.176032]).T
            vtgteci = np.array([- 2.630715586, 1.396719096,- 2.491891]).T
            circ = 'n'
        case 7:
            rtgteci = np.array([- 40588.150362,- 11462.167028, 27.147649]).T
            vtgteci = np.array([0.834787457,- 2.958305691,- 0.001173016]).T
            circ = 'n'
        case 8:
            rtgteci = np.array([- 5551.898646,- 2563.049696, 3257.756165]).T
            vtgteci = np.array([2.149073,- 7.539457,- 2.185709]).T
            circ = 'n'
        case 9:
            rtgteci = np.array([9668.14551571, 6999.07240705, 4041.43303674]).T
            vtgteci = np.array([- 3.65292306, 0.64966519,- 5.82123516]).T
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
    hrokm = hro1.T
    hvo1[0] = xd / 1000.0
    hvo1[1] = yd / 1000.0
    hvo1[2] = zd / 1000.0
    hvokm = hvo1.T
    rintecix, vintecix = EQCM_to_ECI_RTN_sal(rtgteci, vtgteci, hrokm, hvokm)
    rhillx, vhillx = ECI_to_EQCM_RTN_sal(rtgteci, vtgteci, rintecix, vintecix, outfilehill)
    dr = rhillx - hrokm
    dv = vhillx - hvokm
    mdr = sc.mag(dr * 1000000)
    mdv = sc.mag(dv * 1000000)
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
    hrokm = hro1.T
    hvo1[0] = xd / 1000.0
    hvo1[1] = yd / 1000.0
    hvo1[2] = zd / 1000.0
    hvokm = hvo1.T
    rintecisal, vintecisal = EQCM_to_ECI_RTN_sal(rtgteci, vtgteci, hrokm, hvokm)
    print('hillsin       %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
          % (hrokm[0] * 1000, hrokm[1] * 1000, hrokm[2] * 1000,
             hvokm[0] * 1000, hvokm[1] * 1000, hvokm[2] * 1000))
    mmagrti = sc.mag(rtgteci)
    mmagvti = sc.mag(vtgteci)
    print('rtgteci        %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f %15.8f %15.8f \n'
          % (rtgteci[0], rtgteci[1], rtgteci[2], vtgteci[0], vtgteci[1],
             vtgteci[2], mmagrti, mmagvti))
    mmagri = sc.mag(rintecisal)
    mmagvi = sc.mag(vintecisal)
    print('rinteci        %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f  %15.8f %15.8f \n'
          % (rintecisal[0], rintecisal[1], rintecisal[2], vintecisal[0],
             vintecisal[1], vintecisal[2], mmagri, mmagvi))

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
    hrokm = hro1.T
    hvo1[0] = xd / 1000.0
    hvo1[1] = yd / 1000.0
    hvo1[2] = zd / 1000.0
    hvokm = hvo1.T
    # reset this!! so m for hills call!!
    hro = np.zeros(3)
    hvo = np.zeros(3)
    hro[0] = x
    hro[1] = y
    hro[2] = z
    hvo[0] = xd
    hvo[1] = yd
    hvo[2] = zd
    rintecisal, vintecisal = EQCM_to_ECI_RTN_sal(rtgteci, vtgteci,
                                                 hrokm, hvokm)
    rhill2, vhill2 = ECI_to_EQCM_RTN_sal(rtgteci, vtgteci, rintecisal,
                                         vintecisal, outfilehill)
    rintecisal1, vintecisal1 = EQCM_to_ECI_RTN_sal(rtgteci, vtgteci,
                                                   rhill2, vhill2)
    rhillsal, vhillsal = ECI_to_EQCM_RTN_sal(rtgteci, vtgteci, rintecisal1,
                                             vintecisal1, outfilehill)
    rhillsal1, vhillsal1 = ECI_to_EQCM_RTN_sal(rtgteci, vtgteci, rintecisal,
                                               vintecisal, outfilehill)
    # ------------ set the vectors for future use in program
    rinteci = np.zeros(3)
    vinteci = np.zeros(3)
    rinteci[0] = rintecisal[0]
    rinteci[1] = rintecisal[1]
    rinteci[2] = rintecisal[2]
    vinteci[0] = vintecisal[0]
    vinteci[1] = vintecisal[1]
    vinteci[2] = vintecisal[2]
    # set initial values to compare and make sure STK doesn't have an error
    rintecio = rinteci.copy()
    vintecio = vinteci.copy()
    rtgtecio = rtgteci.copy()
    vtgtecio = vtgteci.copy()
    # --------------------- write out various transformation methods --- just at ktr = 8
    if (ktr <= 8):
        print('\nhills in ell   %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (hro[0], hro[1], hro[2], hvo[0], hvo[1], hvo[2]))
        print('\nhills in       %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (hrokm[0] * 1000, hrokm[1] * 1000, hrokm[2] * 1000,
                 hvokm[0] * 1000, hvokm[1] * 1000, hvokm[2] * 1000))
        mmagrti = sc.mag(rtgteci)
        mmagvti = sc.mag(vtgteci)
        print('rtgteci        %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n'
              % (rtgteci[0], rtgteci[1], rtgteci[2], mmagrti, vtgteci[0],
                 vtgteci[1], vtgteci[2], mmagvti))
        print('rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (rhill2[0] * 1000, rhill2[1] * 1000, rhill2[2] * 1000,
                 vhill2[0] * 1000, vhill2[1] * 1000, vhill2[2] * 1000))
        mmagrti = sc.mag(rinteci)
        mmagvti = sc.mag(vinteci)
        print('rinteci        %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n'
              % (rinteci[0], rinteci[1], rinteci[2], mmagrti, vinteci[0],
                 vinteci[1], vinteci[2], mmagvti))
        mmagrti = sc.mag(rintecisal)
        mmagvti = sc.mag(vintecisal)
        print('rintecisal hil %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n'
              % (rintecisal[0], rintecisal[1], rintecisal[2], mmagrti,
                 vintecisal[0], vintecisal[1], vintecisal[2], mmagvti))
        mmagrti = sc.mag(rintecisal1)
        mmagvti = sc.mag(vintecisal1)
        print('rintecisal1 eq %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n'
              % (rintecisal1[0], rintecisal1[1], rintecisal1[2], mmagrti,
                 vintecisal1[0], vintecisal1[1], vintecisal1[2], mmagvti))
        print('hillssal hill  %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (rhillsal[0] * 1000, rhillsal[1] * 1000, rhillsal[2] * 1000,
                 vhillsal[0] * 1000, vhillsal[1] * 1000, vhillsal[2] * 1000))
        print('hillssal1 eq   %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (rhillsal1[0] * 1000, rhillsal1[1] * 1000, rhillsal1[2] * 1000,
                 vhillsal1[0] * 1000, vhillsal1[1] * 1000, vhillsal1[2] * 1000))
        print('==== sals codes ====\n')
        #---            [rintecisal, vintecisal]   = EQCM_to_ECI_NTW_sal( rtgteci*1000, vtgteci*1000, hro', hvo' );      # do in m~!!!find int pos eqcm
        #---            [rhill2, vhill2]          = ECI_to_EQCM_NTW_sal( rtgteci*1000, vtgteci*1000, rintecisal, vintecisal ); # find how much this would be for true hills
        rintecisal, vintecisal = EQCM_to_ECI_RTN_sal(rtgteci, vtgteci,
                                                     hro.T, hvo.T)
        rhill2, vhill2 = ECI_to_EQCM_RTN_sal(rtgteci, vtgteci, rintecisal,
                                             vintecisal, outfilehill)
        mmagrti = sc.mag(rintecisal)
        mmagvti = sc.mag(vintecisal)
        print('rintecisal hil %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n'
              % (rintecisal[0] * 0.001, rintecisal[1] * 0.001,
                 rintecisal[2] * 0.001, mmagrti * 0.001, vintecisal[0] * 0.001,
                 vintecisal[1] * 0.001, vintecisal[2] * 0.001, mmagvti * 0.001))
        print('rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (rhill2[0], rhill2[1], rhill2[2], vhill2[0], vhill2[1], vhill2[2]))

outfilehill.close()
