#
# test the hills conversion programs
#
# Uses:
#  hilleqcm2eci
#  eci2rsw
#  newtonnu, newtone
#  inverselliptic2
#  elliptic12
#  eci2hilleqcm
#
# STK needs to be running (testhillsc)

import numpy as np
import math
import os
import spacemath_utils as smu
from space_constants import *
import space_conversions as sc
import orbit_utils as obu
print('------------------------------------------------ initial accuracy ----------------------------- \n')
# now test the ability to convert eci - hills and back
casenumo = 1

casetest = 7

juststk = 1

ang_step = 1e-08
# case where you read in the int and tgt ephemerides
ropt = 'm'

# clear the contents out between runs w deletes existing data in the file
# outfile = open('d:/STKFiles Educational Files/Hills/testhill.out','wt')
outdir = os.path.join(os.path.dirname(__file__), 'testoutput')
outfile = open(os.path.join(outdir, 'testhill.out'), 'wt')

outfiledet = open(os.path.join(outdir,
                               f'testhilldetc{casenumo}n{casetest}.out'), 'wt')

outfiledethp = open(os.path.join(outdir,
                                 f'testhilldethpc{casenumo}n{casetest}.out'), 'wt')

outfiledetH = open(os.path.join(outdir,
                                f'testhilldetHc{casenumo}n{casetest}.out'), 'wt')

outfiledethpH = open(os.path.join(outdir,
                                  f'testhilldethpHc{casenumo}n{casetest}.out'), 'wt')

outfilehill = open(os.path.join(outdir,
                                f'thillc{casenumo}n{casetest}.out'), 'wt')

# ------------------------- target satellite info ----------------------------
# print('=== be sure to have d:/STKFiles Educational Files/Hills/testhilla.sc open ===\n')
# ---------------------------------------------------------------------
# ------------------------ do initial accuracy checks fwd, back, ------
# ---------------------------------------------------------------------
# --- this one also sets the case for the rest of the program here...

rtgteci = np.zeros(3)
vtgteci = np.zeros(3)
hro1 = np.zeros(3)
hvo1 = np.zeros(3)
hro = np.zeros(3)
hvo = np.zeros(3)
rinteci = np.zeros(3)
vinteci = np.zeros(3)

for ktr in range(1, 12):
    # set last one to be the case of interest
    if ktr < 11:
        casenum = ktr
    else:
        casenum = casenumo

    match casenum:
        case 1:
            rtgteci[0] = (6378.137 + 500.0)
            rtgteci[1] = 0.0
            rtgteci[2] = 0.0
            magrt = smu.mag(rtgteci)
            vtgteci[0] = 0.0
            vtgteci[1] = np.sqrt(mu / magrt)
            vtgteci[2] = 0.0
            circ = 'y'
            #  previous doesn't work in STK92 - fixed in 9.3 for smaller inclination
            a = 6378.137 + 500.0
            ecc = 0.0
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, 0.0, 0.001 * deg2rad,
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            circ = 'y'
            print('rtgt = [%20.13f %20.13f %20.13f]; \n '
                  'vtgt = [%20.13f %20.13f %20.13f]; \n'
                   % (rtgteci[0], rtgteci[1], rtgteci[2], vtgteci[0],
                      vtgteci[1], vtgteci[2]))
        case 2:
            a = 26500.0
            ecc = 0.0
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, 0.0, 0.001 * deg2rad,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            circ = 'y'
        case 3:
            a = 26500.0
            ecc = 0.0
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, 0.0, 55.0 * deg2rad,
                                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            circ = 'y'
        case 4:
            a = 42164.0
            ecc = 0.0
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, 0.0, 0.001 * deg2rad,
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            circ = 'y'
        case 5:
            rtgteci = np.array([-605.7904308, -5870.230407, 3493.052004])
            vtgteci = np.array([-1.568251615, -3.702348353, -6.479484915])
            circ = 'n'
        case 6:
            rtgteci = np.array([-761.24075519, 22699.40899449,
                                13644.176032])
            vtgteci = np.array([-2.630715586, 1.396719096, -2.491891])
            circ = 'n'
        case 7:
            rtgteci = np.array([-40588.150362, -11462.167028, 27.147649])
            vtgteci = np.array([0.834787457, -2.958305691, -0.001173016])
            circ = 'n'
        case 8:
            rtgteci = np.array([- 5551.898646, -2563.049696, 3257.756165])
            vtgteci = np.array([2.149073, -7.539457, -2.185709])
            circ = 'n'
        case 9:
            rtgteci = np.array([9668.14551571, 6999.07240705, 4041.43303674])
            vtgteci = np.array([- 3.65292306, 0.64966519, -5.82123516])
            circ = 'n'
        case 10:
            a = 8164.7188
            ecc = 0.08
            p = a * (1.0 - ecc * ecc)
            rtgteci, vtgteci = sc.coe2rv(p, ecc, 32.861 * deg2rad,
                                         0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            circ = 'n'
    x = 10.0
    y = 10.0
    z = 10.0
    xd = 0.01
    yd = 0.01
    zd = 0.01
    hro1[0] = x / 1000.0
    hro1[1] = y / 1000.0
    hro1[2] = z / 1000.0
    hvo1[0] = xd / 1000.0
    hvo1[1] = yd / 1000.0
    hvo1[2] = zd / 1000.0
    #--- [rintecix, vintecix]   = eqcm2eci( rtgteci, vtgteci, hrokm, hvokm );     # find int pos eqcm
    #--- [rhillx, vhillx]      = eci2eqcm( rtgteci, vtgteci, rintecix, vintecix ); # find how much this would be for true hills

    rintecix, vintecix = sc.hilleqcm2eci(rtgteci, vtgteci, hro1, hvo1)
    rhillx, vhillx = sc.eci2hilleqcm(rtgteci, vtgteci, rintecix, vintecix)
    #        [rintecix, vintecix]   = eqcm2eci( rtgteci, vtgteci, hrokm, hvokm );     # find int pos eqcm
    #        [rhillx, vhillx]      = eci2eqcm( rtgteci, vtgteci, rintecix, vintecix ); # find how much this would be for true hills
    dr = (rhillx - hro1) * 1000
    dv = (vhillx - hvo1) * 1000
    mdr = smu.mag(dr * 1000)
    mdv = smu.mag(dv * 1000)
    print('dr m         %20.13f %20.13f %20.13f %15.8f '
          'mm %20.13f %20.13f %20.13f %15.8f mm/s \n'
          % (dr[0], dr[1], dr[2], mdr, dv[0], dv[1], dv[2], mdv))
    p, a, ecc, incl, omega, argpx1, nux1, m, arglat, truelon, lonper = \
        sc.rv2coe(rtgteci, vtgteci)
    print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f '
          '%11.5f %11.5f %11.5f\n'
          % (p, a, ecc, incl * rad2deg, omega * rad2deg, argpx1 * rad2deg,
             nux1 * rad2deg, m * rad2deg, arglat * rad2deg, truelon * rad2deg,
             lonper * rad2deg))
    p, a, ecc, incl, omega, argpx1, nux1, m, arglat, truelon, lonper = \
        sc.rv2coe(rintecix, vintecix)
    print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f '
          '%11.5f %11.5f %11.5f\n'
          % (p, a, ecc, incl * rad2deg, omega * rad2deg, argpx1 * rad2deg,
             nux1 * rad2deg, m * rad2deg, arglat * rad2deg, truelon * rad2deg,
             lonper * rad2deg))

if (juststk != 0):
    # -----------------------------------------------------------------------------------------------
    # ------------------- check various positions to determine if the velocity is correct ------------
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
        #---        [rinteci, vinteci]   = eqcm2eci( rtgteci, vtgteci, hrokm, hvokm );              # find int pos eqcm
        rinteci, vinteci = sc.hilleqcm2eci(rtgteci, vtgteci, hro1, hvo1)
        print('hillsin  %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n'
              % (hro1[0] * 1000, hro1[1] * 1000, hro1[2] * 1000,
                 hvo1[0] * 1000, hvo1[1] * 1000, hvo1[2] * 1000))
        mmagrti = smu.mag(rtgteci)
        mmagvti = smu.mag(vtgteci)
        print('rtgteci  %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f '
              '%15.8f %15.8f \n'
              % (rtgteci[0], rtgteci[1], rtgteci[2],
                 vtgteci[0], vtgteci[1], vtgteci[2], mmagrti, mmagvti))
        mmagri = smu.mag(rinteci)
        mmagvi = smu.mag(vinteci)
        print('rinteci  %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f  '
              '%15.8f %15.8f \n'
              % (rinteci[0], rinteci[1], rinteci[2],
                 vinteci[0], vinteci[1], vinteci[2], mmagri, mmagvi))

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
        hro[0] = x
        hro[1] = y
        hro[2] = z
        hvo[0] = xd
        hvo[1] = yd
        hvo[2] = zd
        #---        [rinteci, vinteci]   = eqcm2eci( rtgteci, vtgteci, hrokm, hvokm );              # find int pos eqcm
        #---        [rhill2, vhill2]          = eci2eqcm( rtgteci, vtgteci, rinteci, vinteci ); # find how much this would be for true hills
        #---        [rinteci1, vinteci1] = eqcm2eci( rtgteci, vtgteci, rhill2, vhill2 );            # find int that would have found this new hills
        #---        [rhill, vhill]      = eci2eqcm( rtgteci, vtgteci, rinteci1, vinteci1 );  #
        #---        [rhill1, vhill1]    = eci2eqcm( rtgteci, vtgteci, rinteci, vinteci );    #
        rinteci, vinteci = sc.hilleqcm2eci(rtgteci, vtgteci, hro1,
                                                     hvo1)
        rhill2, vhill2 = sc.eci2hilleqcm(rtgteci, vtgteci, rinteci,
                                             vinteci)
        rinteci1, vinteci1 = sc.hilleqcm2eci(rtgteci, vtgteci,
                                                       rhill2, vhill2)
        rhill, vhill = sc.eci2hilleqcm(rtgteci, vtgteci, rinteci1,
                                                 vinteci1)
        rhill1, vhill1 = sc.eci2hilleqcm(rtgteci, vtgteci,
                                                   rinteci, vinteci)
        # hro = rhill2*1000;  # m, used for hills propagtion
        # hvo = vhill2*1000;

        # --------------------- write out various transformation methods --- just at ktr = 8
        if (ktr <= 8):
            print('\nhills in ell   %20.13f %20.13f %20.13f %20.13f %20.13f '
                  '%20.13f \n'
                  % (hro[0], hro[1], hro[2], hvo[0], hvo[1], hvo[2]))
            print('\nhills in       %20.13f %20.13f %20.13f %20.13f %20.13f '
                  '%20.13f \n'
                  % (hro1[0] * 1000, hro1[1] * 1000, hro1[2] * 1000,
                     hvo1[0] * 1000, hvo1[1] * 1000, hvo1[2] * 1000))
            mmagrti = smu.mag(rtgteci)
            mmagvti = smu.mag(vtgteci)
            print('rtgteci        %20.13f %20.13f %20.13f %15.8f %20.13f '
                  '%20.13f %20.13f %15.8f \n'
                  % (rtgteci[0], rtgteci[1], rtgteci[2], mmagrti, vtgteci[0],
                     vtgteci[1], vtgteci[2], mmagvti))
            print('rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f '
                  '%20.13f \n'
                  % (rhill2[0] * 1000, rhill2[1] * 1000, rhill2[2] * 1000,
                     vhill2[0] * 1000, vhill2[1] * 1000, vhill2[2] * 1000))
            mmagrti = smu.mag(rinteci)
            mmagvti = smu.mag(vinteci)
            print('rinteci        %20.13f %20.13f %20.13f %15.8f %20.13f '
                  '%20.13f %20.13f %15.8f \n'
                  % (rinteci[0], rinteci[1], rinteci[2], mmagrti, vinteci[0],
                     vinteci[1], vinteci[2], mmagvti))
            mmagrti = smu.mag(rinteci1)
            mmagvti = smu.mag(vinteci1)
            print('rinteci1 eq %20.13f %20.13f %20.13f %15.8f %20.13f '
                  '%20.13f %20.13f %15.8f \n'
                  % (rinteci1[0], rinteci1[1], rinteci1[2], mmagrti,
                     vinteci1[0], vinteci1[1], vinteci1[2], mmagvti))
            print('hills hill  %20.13f %20.13f %20.13f %20.13f %20.13f '
                  '%20.13f \n'
                  % (rhill[0] * 1000, rhill[1] * 1000,
                     rhill[2] * 1000, vhill[0] * 1000,
                     vhill[1] * 1000, vhill[2] * 1000))
            print('hills1 eq   %20.13f %20.13f %20.13f %20.13f %20.13f '
                  '%20.13f \n'
                  % (rhill1[0] * 1000, rhill1[1] * 1000,
                     rhill1[2] * 1000, vhill1[0] * 1000,
                     vhill1[1] * 1000, vhill1[2] * 1000))
            print('==== s codes ====\n' % ())
            #---            [rinteci, vinteci]   = EQCM_to_ECI_NTW_( rtgteci*1000, vtgteci*1000, hro', hvo' );      # do in m~!!!find int pos eqcm
            #---            [rhill2, vhill2]          = ECI_to_EQCM_NTW_( rtgteci*1000, vtgteci*1000, rinteci, vinteci ); # find how much this would be for true hills
            rinteci, vinteci = sc.hilleqcm2eci(rtgteci, vtgteci,
                                                         hro, hvo)
            rhill2, vhill2 = sc.eci2hilleqcm(rtgteci, vtgteci, rinteci,
                                                 vinteci)
            mmagrti = smu.mag(rinteci)
            mmagvti = smu.mag(vinteci)
            print('rinteci hil %20.13f %20.13f %20.13f %15.8f %20.13f '
                  '%20.13f %20.13f %15.8f \n'
                  % (rinteci[0] * 0.001, rinteci[1] * 0.001,
                     rinteci[2] * 0.001, mmagrti * 0.001,
                     vinteci[0] * 0.001, vinteci[1] * 0.001,
                     vinteci[2] * 0.001, mmagvti * 0.001))
            print('rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f '
                  '%20.13f \n' % (rhill2[0], rhill2[1], rhill2[2], vhill2[0],
                                  vhill2[1], vhill2[2]))


# -----------------------------------------------------------------------------------------------
# -------------------------------------- now do it in STK ---------------------------------------
# -----------------------------------------------------------------------------------------------
print('\n\n-------------------------------------------------- now in STK -------------------------------------------- \n')
ktr8str = '10.0  0.01'
ktr9str = '100.0  0.01'
ktr10str = '1000.0  0.01'
ktr11str = '10.0  0.1'
ktr12str = '100.0  0.1'
ktr13str = '1000.0  0.1'
ktr14str = '10.0  1.0'
ktr15str = '100.0  1.0'
ktr16str = '1000.0  1.0'
first = 0

for ktr in range(casetest, casetest + 1):
    if ktr == 7:
        x = 100000.0
        y = 2348.0
        z = - 345.0
        xd = 1
        yd = 2.34
        zd = 0.3829
        currktrstr = 'misc'
    if ktr == 8:
        x = 10.0
        y = 10.0
        z = 10.0
        xd = 0.01
        yd = 0.01
        zd = 0.01
        #             x = 1000.0;  # in m
        #             y = 1000.0;#100.0;  # in m
        #             z = 1000.0;#100.0;  # in m
        #             xd = 0.1000;#0.10;
        #             yd = 0.1000;#0.10;
        #             zd = 0.1000;#0.10;
        currktrstr = '10.0  0.01'
    if ktr == 9:
        x = 100.0
        y = 100.0
        z = 100.0
        xd = 0.01
        yd = 0.01
        zd = 0.01
        currktrstr = '100.0  0.01'
    if ktr == 10:
        x = 1000.0
        y = 1000.0
        z = 1000.0
        xd = 0.01
        yd = 0.01
        zd = 0.01
        currktrstr = '1000.0  0.01'
    if ktr == 11:
        x = 10.0
        y = 10.0
        z = 10.0
        xd = 0.1
        yd = 0.1
        zd = 0.1
        scalef = x * y * z / 1000000.0
        currktrstr = '10.0  0.1'
    if ktr == 12:
        x = 100.0
        y = 100.0
        z = 100.0
        xd = 0.1
        yd = 0.1
        zd = 0.1
        currktrstr = '100.0  0.1'
    if ktr == 13:
        x = 1000.0
        y = 1000.0
        z = 1000.0
        xd = 0.1
        yd = 0.1
        zd = 0.1
        currktrstr = '1000.0  0.1'
    if ktr == 14:
        x = 10.0
        y = 10.0
        z = 10.0
        xd = 1.0
        yd = 1.0
        zd = 1.0
        currktrstr = '10.0  1.0'
    if ktr == 15:
        x = 100.0
        y = 0.0
        z = 0.0
        xd = 0.0
        yd = 0.0
        zd = 0.0
        currktrstr = '100.0  1.0'
    if ktr == 16:
        x = 10000.0
        y = 10000.0
        z = 10000.0
        xd = 1.0
        yd = 1.0
        zd = 1.0
        currktrstr = '1000.0  1.0'
    hro1[0] = x / 1000.0
    hro1[1] = y / 1000.0
    hro1[2] = z / 1000.0
    hvo1[0] = xd / 1000.0
    hvo1[1] = yd / 1000.0
    hvo1[2] = zd / 1000.0
    #        [rinteci, vinteci]   = eqcm2eci( rtgteci, vtgteci, hrokm, hvokm );              # find int pos eqcm
    #old way        [rinteci, vinteci]   = hill2eci( rtgteci, vtgteci, hrokm, hvokm );              # find int pos eqcm
    rinteci, vinteci = sc.hilleqcm2eci(rtgteci, vtgteci, hro1, hvo1)
    #        pause;
    #---        [rhill2, vhill2]          = eci2eqcm( rtgteci, vtgteci, rinteci, vinteci ); # find how much this would be for true hills
    rhill2, vhill2 = sc.eci2hilleqcm(rtgteci, vtgteci, rinteci, vinteci)

    #        if ktr >= 8  # cases for paper
    # get the vectors to STK, produce .e files, and then process the
    # differences for excel
    if first == 0:
        outfile.write(' ============================ Case for book %i \n'
                      % (casenumo))
        mmagrti = smu.mag(rinteci)
        mmagvti = smu.mag(vinteci)
        outfile.write('rint      %20.13f %20.13f %20.13f %15.8f %20.13f '
                      '%20.13f %20.13f %15.8f \n'
                      % (rinteci[0], rinteci[1], rinteci[2], mmagrti,
                         vinteci[0], vinteci[1], vinteci[2], mmagvti))
        mmagrti = smu.mag(rtgteci)
        mmagvti = smu.mag(vtgteci)
        outfile.write('rtgt      %20.13f %20.13f %20.13f %15.8f %20.13f '
                      '%20.13f %20.13f %15.8f \n'
                      % (rtgteci[0], rtgteci[1], rtgteci[2], mmagrti,
                         vtgteci[0], vtgteci[1], vtgteci[2], mmagvti))
        outfile.write('rhill         %20.13f %20.13f %20.13f   %20.13f %20.13f'
                      ' %20.13f   \n'
                      % (rhill2[0], rhill2[1], rhill2[2], vhill2[0],
                         vhill2[1], vhill2[2]))
        outfiledet.write(' ============================ Case for book %i \n'
                         % (casenumo))
        mmagrti = smu.mag(rinteci)
        mmagvti = smu.mag(vinteci)
        outfiledet.write('rint, ,     %20.13f, %20.13f, %20.13f, %20.13f, '
                         '%20.13f, %20.13f, %15.8f %15.8f \n'
                         % (rinteci[0], rinteci[1], rinteci[2], vinteci[0],
                            vinteci[1], vinteci[2], mmagrti, mmagvti))
        mmagrti = smu.mag(rtgteci)
        mmagvti = smu.mag(vtgteci)
        outfiledet.write('rtgt,  ,    %20.13f, %20.13f, %20.13f, %20.13f, '
                         '%20.13f, %20.13f, %15.8f %15.8f \n'
                         % (rtgteci[0], rtgteci[1], rtgteci[2], vtgteci[0],
                            vtgteci[1], vtgteci[2], mmagrti, mmagvti))
        outfiledet.write('rhill ,  ,      %20.13f, %20.13f, %20.13f,   '
                         '%20.13f, %20.13f, %20.13f   \n'
                         % (rhill2[0], rhill2[1], rhill2[2], vhill2[0],
                            vhill2[1], vhill2[2]))
        outfiledethp.write(' ============================ Case for book %i \n'
                           % (casenumo))
        mmagrti = smu.mag(rinteci)
        mmagvti = smu.mag(vinteci)
        outfiledethp.write('rint,  ,    %20.13f, %20.13f, %20.13f, %20.13f, '
                           '%20.13f, %20.13f, %15.8f %15.8f \n'
                           % (rinteci[0], rinteci[1], rinteci[2], vinteci[0],
                              vinteci[1], vinteci[2], mmagrti, mmagvti))
        mmagrti = smu.mag(rtgteci)
        mmagvti = smu.mag(vtgteci)
        outfiledethp.write('rtgt,  ,    %20.13f, %20.13f, %20.13f, %20.13f, '
                           '%20.13f, %20.13f, %15.8f %15.8f \n'
                           % (rtgteci[0], rtgteci[1], rtgteci[2], vtgteci[0],
                              vtgteci[1], vtgteci[2], mmagrti, mmagvti))
        outfiledethp.write('rhill ,  ,      %20.13f, %20.13f, %20.13f,   '
                           '%20.13f, %20.13f, %20.13f   \n'
                           % (rhill2[0], rhill2[1], rhill2[2], vhill2[0],
                              vhill2[1], vhill2[2]))
        first = 1
    else:
        outfile.write(' \n')
        outfiledet.write(' \n')
        outfiledethp.write(' \n')

#     # commenting out stk stuff for now. -zeg
#     # ---- Attach to running STK (stkeducationalfiles\testhilla.sc)
#     app = actxGetRunningServer('STK.application')
#     root = app.get('Personality2')
#     scen = root.CurrentScenario
#     # ---------------------------------- set these up one time ----------------------------------
#     dtsec = 180.0
#     numsteps = 1440
#     # note that the .e files go for 4 days, or 5760 min
#     #          fprintf(1,'start setting stk sats  - please wait... \n' );
#     StartT = '1 jun 2011 00:00:00.000'
#     StopT = '5 jun 2011 00:00:00.000'
#     sat = scen.Children.Item('InterceptorHPOP')
#     #sat.SetPropagatorType( 'ePropagatorHPOP' );
#     sat.Propagator.StartTime = StartT
#     sat.Propagator.StopTime = StopT
#     sat.Propagator.Step = dtsec
#     sat.Propagator.InitialState.Representation.AssignCartesian('eCoordinateSystemICRF', rinteci[0], rinteci[1], rinteci[2], vinteci[0], vinteci[1], vinteci[2])
#     sat.Propagator.Propagate
#     # Export ephemeris
#     ephtool = sat.ExportTools.GetEphemerisStkExportTool
#     ephtool.CoordinateSystem = 'eStkEphemCoordinateSystemICRF'
#     ephtool.Export('d:\STKFiles Educational Files\Hills\inthpop.txt')
#     sat = scen.Children.Item('Interceptor')
#     sat.SetPropagatorType('ePropagatorTwoBody')
#     sat.Propagator.StartTime = StartT
#     sat.Propagator.StopTime = StopT
#     sat.Propagator.Step = dtsec
#     sat.Propagator.InitialState.Representation.AssignCartesian('eCoordinateSystemICRF', rinteci[0], rinteci[1], rinteci[2], vinteci[0], vinteci[1], vinteci[2])
#     sat.Propagator.Propagate
#     # Export ephemeris
#     ephtool = sat.ExportTools.GetEphemerisStkExportTool
#     ephtool.CoordinateSystem = 'eStkEphemCoordinateSystemICRF'
#     ephtool.Export('d:\STKFiles Educational Files\Hills\int.txt')
#     sat = scen.Children.Item('TargetHPOP')
#     #sat.SetPropagatorType( 'ePropagatorHPOP' );
#     sat.Propagator.StartTime = StartT
#     sat.Propagator.StopTime = StopT
#     sat.Propagator.Step = dtsec
#     sat.Propagator.InitialState.Representation.AssignCartesian('eCoordinateSystemICRF', rtgteci[0], rtgteci[1], rtgteci[2], vtgteci[0], vtgteci[1], vtgteci[2])
#     sat.Propagator.Propagate
#     # Export ephemeris
#     ephtool = sat.ExportTools.GetEphemerisStkExportTool
#     ephtool.CoordinateSystem = 'eStkEphemCoordinateSystemICRF'
#     ephtool.Export('d:\STKFiles Educational Files\Hills\tgthpop.txt')
#     sat = scen.Children.Item('Target')
#     #sat.SetPropagatorType( 'ePropagatorHPOP' );
#     sat.Propagator.StartTime = StartT
#     sat.Propagator.StopTime = StopT
#     sat.Propagator.Step = dtsec
#     sat.Propagator.InitialState.Representation.AssignCartesian('eCoordinateSystemICRF', rtgteci[0], rtgteci[1], rtgteci[2], vtgteci[0], vtgteci[1], vtgteci[2])
#     sat.Propagator.Propagate
#     # Export ephemeris
#     ephtool = sat.ExportTools.GetEphemerisStkExportTool
#     ephtool.CoordinateSystem = 'eStkEphemCoordinateSystemICRF'
#     ephtool.Export('d:\STKFiles Educational Files\Hills\tgt.txt')
#     print('case %i, set maually \n' % (casenumo))
#     #     pause;
#     #          fprintf(1,'done with setting stk sats  \n' );
#     # goto compiled program!!
#     print('\n\n------------------------------- done in stk, go to compiled program? ---------------------------------------- \n')
#     res1 = readdote('d:\STKFiles Educational Files\Hills\int.txt')
#     res2 = readdote('d:\STKFiles Educational Files\Hills\inthpop.txt')
#     res3 = readdote('d:\STKFiles Educational Files\Hills\tgt.txt')
#     res4 = readdote('d:\STKFiles Educational Files\Hills\tgthpop.txt')
#     # Check if stk mixes up state for low ecc, etc
#     #             if (mag(rintecio-res1.pos(:, 1)'*0.001) > 0.001)
#     #                 fprintf(1,'error in rinteci, set maually #20.13f #20.13f #20.13f  \n', rintecio, res1.pos(:, 1)*0.001 );
#     #                 pause;
#     #                 #    sat = scen.Children.Item( 'Interceptor' );  # Get current satellite
#     #                 #    sat.Propagator.Propagate;
#     #                 [res1] = readdote ( 'd:\STK Educational Files\Hills\int.txt' );  # all these files are in meters
#     #             end
#     #             if (mag(rintecio-res2.pos(:, 1)'*0.001) > 0.001)
#     #                 fprintf(1,'error in H rinteci, set maually #20.13f #20.13f #20.13f \n', rintecio, res2.pos(:, 1)*0.001 );
#     #             end
#     #             if (mag(rtgtecio-res3.pos(:, 1)*0.001) > 0.001)
#     #                 fprintf(1,'error in rtgteci, set maually #20.13f #20.13f #20.13f \n', rtgtecio, res3.pos(:, 1)*0.001 );
#     #             end
#     #             if (mag(rtgtecio-res4.pos(:, 1)*0.001) > 0.001)
#     #                 fprintf(1,'error in H rtgteci, set maually #20.13f #20.13f #20.13f  \n', rtgtecio, res4.pos(:, 1)*0.001 );
#     #             end
#     #             if (mag(vintecio-res1.vel(:, 1)'*0.001) > 0.1)
#     #                 fprintf(1,'error in vinteci, set maually #20.13f #20.13f #20.13f \n', vintecio, res1.vel(:, 1)*0.001 );
#     #                 pause;
#     #                 #    sat = scen.Children.Item( 'Interceptor' );  # Get current satellite
#     #                 #    sat.Propagator.Propagate;
#     #                 [res1] = readdote ( 'd:\STK Educational Files\Hills\int.txt' );  # all these files are in meters
#     #             end
#     #             if (mag(vintecio-res2.vel(:, 1)'*0.001) > 0.1)
#     #                 fprintf(1,'error in H vinteci, set maually #20.13f #20.13f #20.13f \n', vintecio, res2.vel(:, 1)*0.001 );
#     #             end
#     #             if (mag(vtgtecio-res3.vel(:, 1)*0.001) > 0.1)
#     #                 fprintf(1,'error in vtgteci, set maually #20.13f #20.13f #20.13f \n', vtgtecio, res3.vel(:, 1)*0.001 );
#     #             end
#     #             if (mag(vtgtecio-res4.vel(:, 1)*0.001) > 0.1)
#     #                 fprintf(1,'error in H vtgteci, set maually #20.13f #20.13f #20.13f \n', vtgtecio, res4.vel(:, 1)*0.001 );
#     #             end
#     print('start stepping through times and finding diffs  \n')
#     # setup initial hills displacement from various sources (above)
#     kk = 1
#     # use the linear approximation only
#     hrost = hrokm
#     hvost = hvokm
#     ri = res3.pos[:, kk]
#     # xxx m doesn't work
#     #            hrost = hro;  # m
#     #            hvost = hvo;
#     # use the two-body values
#     #            [hrost, hvost] = eci2eqcm( res3.pos(:, kk), res3.vel(:, kk), res1.pos(:, kk), res1.vel(:, kk) ); # two-body
#     #             hrost = hrost / 1000;
#     #             hvost = hvost / 1000;
#     # units must be kmxxx hpop
#     #            ri = res3.pos(:, kk);  # must be the target!!!!
#     # use the hpop values
#     #            [hrost, hvost] = eci2eqcm(res4.pos(:, kk), res4.vel(:, kk), res2.pos(:, kk), res2.vel(:, kk) );
#     #             hrost = hrost / 1000;
#     #             hvost = hvost / 1000;
#     # units must be kmxxx hpop
#     #            ri = res4.pos(:, kk);  # must be the target!!!!
#     magralt = (smu.mag(ri) - 6378137.0) * 0.001
#     # select various initial hill's estimates

#     for kk in range(numsteps):
#         #load data into x y z arrays
#         rinteci2 = res1.pos[:, kk] / 1000.0
#         vinteci2 = res1.vel[:, kk] / 1000.0
#         ttime = res1.t[:, kk]
#         rintecih = res2.pos[:, kk] / 1000.0
#         vintecih = res2.vel[:, kk] / 1000.0
#         rtgteci2 = res3.pos[:, kk] / 1000.0
#         vtgteci2 = res3.vel[:, kk] / 1000.0
#         rtgtecih = res4.pos[:, kk] / 1000.0
#         vtgtecih = res4.vel[:, kk] / 1000.0
#         # do conversions at each step from the numerical ephemerides
#         #---                [rhill2body, vhill2body] = eci2eqcm( rtgteci2*0.001, vtgteci2*0.001, rinteci2*0.001, vinteci2*0.001 ); # in km
#         #---                [rhillhpop, vhillhpop]   = eci2eqcm( rtgtecih*0.001, vtgtecih*0.001, rintecih*0.001, vintecih*0.001 );
#         rhill2body, vhill2body = ECI_to_EQCM_RTN_(rtgteci2, vtgteci2, rinteci2, vinteci2, outfilehill)
#         rhillhpop, vhillhpop = ECI_to_EQCM_RTN_(rtgtecih, vtgtecih, rintecih, vintecih, outfilehill)
#         #rhill2body
#         # xxx
#         # if , then change 1 more asignments below and hrost above...
#         #               [rhill2body, vhill2body] = ECI_to_EQCM_NTW_( rtgteci2, vtgteci2, rinteci2, vinteci2 ); # in m
#         #               [rhillhpop, vhillhpop]   = ECI_to_EQCM_NTW_( rtgtecih, vtgtecih, rintecih, vintecih );
#         rhill2body = rhill2body.T
#         vhill2body = vhill2body.T
#         rhillhpop = rhillhpop.T
#         vhillhpop = vhillhpop.T
#         # use Hills with various initial conditions
#         rhill, vhill = obu.hillsr(hrost, hvost, magralt, ttime)
#         # store rhill for later plotting against numerical versions
#         rhillarr = np.zeros((kk, 3))
#         rhillhpoparr = np.zeros((kk, 3))
#         rhill2barr = np.zeros((kk, 3))
#         vhillarr = np.zeros((kk, 3))
#         vhillhpoparr = np.zeros((kk, 3))
#         vhill2barr = np.zeros((kk, 3))

#         rhillarr[kk, 0] = rhill[0]
#         rhillarr[kk, 1] = rhill[1]
#         rhillarr[kk, 2] = rhill[2]
#         rhillhpoparr[kk, 0] = rhillhpop[0]
#         rhillhpoparr[kk, 1] = rhillhpop[1]
#         rhillhpoparr[kk, 2] = rhillhpop[2]
#         rhill2barr[kk, 0] = rhill2body[0]
#         rhill2barr[kk, 1] = rhill2body[1]
#         rhill2barr[kk, 2] = rhill2body[2]
#         vhillarr[kk, 0] = vhill[0]
#         vhillarr[kk, 1] = vhill[1]
#         vhillarr[kk, 2] = vhill[2]
#         vhillhpoparr[kk, 0] = vhillhpop[0]
#         vhillhpoparr[kk, 1] = vhillhpop[1]
#         vhillhpoparr[kk, 2] = vhillhpop[2]
#         vhill2barr[kk, 0] = vhill2body[0]
#         vhill2barr[kk, 1] = vhill2body[1]
#         vhill2barr[kk, 2] = vhill2body[2]
#         dr = rhill2body - rhill
#         dv = vhill2body - vhill
#         drh = rhillhpop - rhill
#         dvh = vhillhpop - vhill
#         if (math.remainder(kk, 50) == 0) or (kk == 1):
#             #if (kk>798) && (kk <801)
#             print('%6i,%12.3f  h2b %6.3f %6.3f %6.3f, h %6.3f %6.3f %6.3f, dr %6.3f %6.3f %6.3f, int %11.3f %11.3f %11.3f, tgt %11.3f %11.3f %11.3f \n' % (kk, ttime, rhill2body[0], rhill2body[1], rhill2body[2], rhill[0], rhill[1], rhill[2], dr[0], dr[1], dr[2], rinteci2[0], rinteci2[1], rinteci2[2], rtgteci2[0], rtgteci2[1], rtgteci2[2]))
#             #     dbstop in testhilldav at 771
#         #     1,   0.000  h2b  0.100  0.000  0.000, h  0.100  0.000  0.000, int 6878237.000       0.000       0.000, tgt 6878137.000       0.000       0.000
#         #     2, 180.000  h2b  0.102 -0.020  0.000, h  0.106 -0.001  0.000, int 6742198.558 1361223.730      23.758, tgt 6742094.602 1361223.334      23.758
#         #     3, 360.000  h2b  0.108 -0.042 -0.000, h  0.124 -0.006  0.000, int 6339464.620 2668602.690      46.576, tgt 6339348.963 2668599.503      46.576
#         #     4, 540.000  h2b  0.117 -0.067 -0.000, h  0.152 -0.021  0.000, int 5685966.422 3870422.267      67.552, tgt 5685831.863 3870411.426      67.551
#         #     5, 720.000  h2b  0.130 -0.096 -0.000, h  0.190 -0.049  0.000, int 4807554.692 4919143.853      85.855, tgt 4807395.079 4919117.922      85.855
#         # xxx
#         dr = dr * 1000
#         dv = dv * 1000
#         drh = drh * 1000
#         dvh = dvh * 1000
#         #                 drpp1 = rinteci2 - rintecih; # overall difference in vectors, 2body and hpop
#         #                 drpp2 = rtgteci2 - rtgtecih;
#         #                 dvpp1 = vinteci2 - vintecih; # overall difference in vectors, 2body and hpop
#         #                 dvpp2 = vtgteci2 - vtgtecih;
#         #                 drp(kk, 3) = mag(drpp1);
#         #                 drp(kk, 4) = mag(drpp2);
#         #                 dvp(kk, 3) = mag(dvpp1);
#         #                 dvp(kk, 4) = mag(dvpp2);
#         dra = np.zeros((kk, 4))
#         drha = np.zeros((kk, 4))
#         dva = np.zeros((kk, 4))
#         dvha = np.zeros((kk, 4))

#         dra[kk, 0] = dr[0]
#         dra[kk, 1] = dr[1]
#         dra[kk, 2] = dr[2]
#         dra[kk, 3] = smu.mag(dr)
#         drha[kk, 0] = drh[0]
#         drha[kk, 1] = drh[1]
#         drha[kk, 2] = drh[2]
#         drha[kk, 3] = smu.mag(drh)
#         dva[kk, 0] = dv[0]
#         dva[kk, 1] = dv[1]
#         dva[kk, 2] = dv[2]
#         dva[kk, 3] = smu.mag(dv)
#         dvha[kk, 0] = dvh[0]
#         dvha[kk, 1] = dvh[1]
#         dvha[kk, 2] = dvh[2]
#         dvha[kk, 3] = smu.mag(dvh)
#         outfiledet.write('diff 2body, %6i, %20.13f, %20.13f, %20.13f, %20.13f,'
#                          ' %20.13f, %20.13f \n'
#                          % (ttime, dr[0], dr[1], dr[2], dv[0], dv[1], dv[2]))
#         outfiledethp.write('diff hpop,  %6i, %20.13f, %20.13f, %20.13f, '
#                            '%20.13f, %20.13f, %20.13f \n'
#                            % (ttime, drh[0], drh[1], drh[2], dvh[0], dvh[1],
#                               dvh[2]))
#         outfiledetH.write('Hill hill,  %6i, %20.13f, %20.13f, %20.13f, '
#                           '%20.13f, %20.13f, %20.13f \n'
#                            % (ttime, rhill[0], rhill[1], rhill[2], vhill[0],
#                               vhill[1], vhill[2]))
#         outfiledethpH.write('Hill hpop,  %6i, %20.13f, %20.13f, %20.13f, '
#                             '%20.13f, %20.13f, %20.13f \n'
#                             % (ttime, rhillhpop[0], rhillhpop[1], rhillhpop[2],
#                                vhillhpop[0], vhillhpop[1], vhillhpop[2]))
#         #   [rhill2body, vhill2body]      = ECI_to_EQCM_RTN_( rtgteci2, vtgteci2, rinteci2, vinteci2 ); # find how much this would be for true hills
#         rx = np.zeros(3)
#         vx = np.zeros(3)
#         rx[0] = rhill2body[0]
#         rx[1] = rhill2body[1]
#         rx[2] = rhill2body[2]
#         vx[0] = rhill2body[0]
#         vx[1] = rhill2body[1]
#         vx[2] = rhill2body[2]
#         rintecisl, vintecisl = EQCM_to_ECI_RTN_(rtgteci2, vtgteci2, rx, vx)
#         outfilehill.write(' %6i %20.13f %20.13f %20.13f %20.13f %20.13f '
#                           '%20.13f \n'
#                           % (ttime, rintecisl[0], rintecisl[1], rintecisl[2],
#                              vintecisl[0], vintecisl[1], vintecisl[2]))

# outfiledet.close()
# outfiledethp.close()
# outfiledetH.close()
# outfiledethpH.close()
# outfilehill.close()
# #pause;

# #         res1.t(1:numsteps) = res1.t(1:numsteps) / 86400.0;
# #         #res1.t(numsteps)

# #         #    figure[0];
# #         #    plot(res1.t(1:numsteps), dra(:, 4),'b-','LineWidth', 1.5); # magnitude
# #         #    hold on;


# #         if ktr == 8 || ktr == 1
# #             color1 = 'b-';
# #         end;
# #         if ktr == 9 || ktr == 2
# #             color1 = 'c-';
# #         end;
# #         if ktr == 10 || ktr == 3
# #             color1 = 'r-';
# #         end;
# #         if ktr == 11
# #             color1 = 'm-';
# #         end;
# #         if ktr == 12
# #             color1 = 'g-';
# #         end;
# #         if ktr == 13
# #             color1 = 'y-';
# #         end;
# #         #    plot(res1.t(1:numsteps), dra(:, 1),'b-'); #
# #         #    plot(res1.t(1:numsteps), dra(:, 2),'b-.'); #
# #         #    plot(res1.t(1:numsteps), dra(:, 3),'b-.'); #
# #         figure[0];
# #         #plot(res1.t(1:numsteps), dra(:, 4), color1,'LineWidth', 1.5); # magnitude
# #         #hold on;
# #         nummpoints = 10;
# #         mytemp = reshape( dra(:, 4), numsteps/nummpoints, nummpoints );
# #         mytime = reshape( res1.t(1:numsteps), numsteps/nummpoints, nummpoints );
# #         errorbar( mean(mytime), mean(mytemp), std(mytemp), color1 )
# #         legend( ktr8str, ktr9str, ktr10str, ktr11str, ktr12str, ktr13str, 'Location', 'NorthWest' )

# #         fprintf(1,'kep  case #6.2f #6.2f  #10.3f #10.3f #10.3f #10.3f #10.3f #10.3f \n', ...
# #             x, xd, mean(dra(1:240, 4)), mean(dra(241:480, 4)), mean(dra(481:720, 4)), mean(dra(721:960, 4)), mean(dra(961:1200, 4)), mean(dra(1201:1440, 4)) );
# #         fprintf(1,'kep  case #6.2f #6.2f  #10.3f #10.3f #10.3f #10.3f #10.3f #10.3f \n', ...
# #             x, xd, std(dra(1:240, 4)), std(dra(241:480, 4)), std(dra(481:720, 4)), std(dra(721:960, 4)), std(dra(961:1200, 4)), std(dra(1201:1440, 4)) );
# #         p = polyfit( mean(mytime), mean(mytemp), 1); # find the polynomial fit
# #         fprintf(outfile,'kep  #s  p = #f (t) + #f \n', currktrstr, p[0], p[1] )
# #         if ktr == 8
# #             # xlabel('time (days)');
# #             # ylabel('error (m)');
# #             xlabel(gcf,'String','time (days)','FontName','Times New Roman','FontWeight','normal','FontSize', 9);
# #             ylabel(gcf,'String','error (m)','FontName','Times New Roman','FontWeight','normal','FontSize', 9);
# #             #        set(gca,'YLim',[1.0, 1000.0]);
# #             # set(gca,'YScale','log');
# #             set(gca,'XMinorTick','on');
# #             set(gca,'YMinorTick','on');
# #             hold on;
# #             drawnow;
# #         end
# #         #     # ---- for first graph
# #         #     figure;
# #         #     plot(res1.t(1:numsteps), dra(:, 1),'r-'); #
# #         #     hold on
# #         #     plot(res1.t(1:numsteps), dra(:, 2),'c-'); #
# #         #     plot(res1.t(1:numsteps), dra(:, 3),'m-'); #
# #         #     plot(res1.t(1:numsteps), dra(:, 4),'b-'); #
# #         #     numpoints = 10;
# #         #     mytemp = reshape( dra(:, 4), numsteps/nummpoints, nummpoints );
# #         #     mytime = reshape( res1.t(1:numsteps), numsteps/nummpoints, nummpoints );
# #         #     errorbar( mean(mytime), mean(mytemp), std(mytemp), 'g-','LineWidth', 1.5 )
# #         #     legend( 'Radial', 'Along-track','Cross-track','Magnitude', 'Errorbars', 'Location', 'NorthWest' )

# #         # check the tgt 2 tgt differences
# #         #mytemp = reshape( drp(:, 4), numsteps/nummpoints, nummpoints );
# #         #errorbar( mean(mytime), mean(mytemp), std(mytemp), color1 );



# #         #    plot(res1.t(1:numsteps), drha(:, 1),'m-'); #
# #         #    plot(res1.t(1:numsteps), drha(:, 2),'m-.'); #
# #         #    plot(res1.t(1:numsteps), drha(:, 3),'m-.'); #
# #         figure[1];
# #         #plot(res1.t(1:numsteps), drha(:, 4), color1,'LineWidth', 1.5); # magnitude
# #         #hold on;
# #         mytemp = reshape( drha(:, 4), numsteps/nummpoints, nummpoints );
# #         mytime = reshape( res1.t(1:numsteps), numsteps/nummpoints, nummpoints );
# #         errorbar( mean(mytime), mean(mytemp), std(mytemp), color1 )
# #         legend( ktr8str, ktr9str, ktr10str, ktr11str, ktr12str, ktr13str, 'Location', 'NorthWest' )

# #         fprintf(1,'hpop case #6.2f #6.2f  #10.3f #10.3f #10.3f #10.3f #10.3f #10.3f \n', ...
# #             x, xd, mean(drha(1:240, 4)), mean(drha(241:480, 4)), mean(drha(481:720, 4)), mean(drha(721:960, 4)), mean(drha(961:1200, 4)), mean(drha(1201:1440, 4)) );
# #         fprintf(1,'hpop case #6.2f #6.2f  #10.3f #10.3f #10.3f #10.3f #10.3f #10.3f \n', ...
# #             x, xd, std(drha(1:240, 4)), std(drha(241:480, 4)), std(drha(481:720, 4)), std(drha(721:960, 4)), std(drha(961:1200, 4)), std(drha(1201:1440, 4)) );
# #         p = polyfit( mean(mytime), mean(mytemp), 1);
# #         fprintf(outfile,'hpop #s  p = #f (t) + #f \n', currktrstr, p[0], p[1] )
# #         if ktr == 8
# #             # xlabel('time (days)');
# #             # ylabel('error (m)');
# #             xlabel(gcf,'String','time (days)','FontName','Times New Roman','FontWeight','normal','FontSize', 9);
# #             ylabel(gcf,'String','error (m)','FontName','Times New Roman','FontWeight','normal','FontSize', 9);
# #             #        set(gca,'YLim',[1.0, 1000.0]);
# #             # set(gca,'YScale','log');
# #             set(gca,'XMinorTick','on');
# #             set(gca,'YMinorTick','on');
# #             hold on;
# #             drawnow;
# #         end

# #         if ktr == 888
# #             figure;
# #             plot(rhillhpoparr(:, 2), rhillhpoparr(:, 1),'r-');
# #             hold on;
# #             plot(rhillarr(:, 2), rhillarr(:, 1),'b-');

# #             figure;
# #             plot(rhillhpoparr(:, 2), rhillhpoparr(:, 1),'r-');
# #             hold on;
# #             plot(rhillarr(:, 2), rhillarr(:, 1),'b-');
# #             axis equal;

# #             figure;
# #             plot(rhillhpoparr(:, 2), rhillhpoparr(:, 3),'r-');
# #             hold on;
# #             plot(rhillarr(:, 2), rhillarr(:, 3),'b-');

# #             figure;
# #             plot(vhillhpoparr(:, 2), vhillhpoparr(:, 1),'r-');
# #             hold on;
# #             plot(vhillarr(:, 2), vhillarr(:, 1),'b-');

# #             figure;
# #             plot(vhillhpoparr(:, 2), vhillhpoparr(:, 1),'r-');
# #             hold on;
# #             plot(vhillarr(:, 2), vhillarr(:, 1),'b-');
# #             axis equal;

# #             figure;
# #             plot(vhillhpoparr(:, 2), vhillhpoparr(:, 3),'r-');
# #             hold on;
# #             plot(vhillarr(:, 2), vhillarr(:, 3),'b-');

# #         end

# # old end of loop...

# print('finished with diffs \n')
# # draw error bars for a segment and at a point in time
# #errorbar(.25, mean(dra(1:240, 4)), std(dra(1:240, 4)),'b-','LineWidth', 1.5); # magnitude

# #  mytemp = reshape( dra(:, 4), 1440/20, 20 );
# #  mytime = reshape( res1.t(1:numsteps), 1440/20, 20 );
# #  errorbar( mean(mytime), mean(mytemp), std(mytemp), 'b-' )
# #  legend( 'plot1', 'plot2','plot3','plot34', 'Location', 'NorthWest' )

# print('results in testouput/testhill.out  \n')
# print('results in testoutput/testhilldet.out  \n')
# print('should be done now. sort the outfile and plot in excel \n')
# outfile.close()
# #   load handel
# #   sound(y, Fs)

# ------------------------------------ test support question for hills
rtgteci = np.array([0.0, 9380.502297, 5093.197188])
vtgteci = np.array([- 6.110899, 3.222217, 1.749521])
p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = sc.rv2coe(rtgteci, vtgteci)
print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n'
      % (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg,
         nu * rad2deg, m * rad2deg, arglat * rad2deg, truelon * rad2deg,
         lonper * rad2deg))
hrokm = np.array([79.37875, 0.0, -93.843398])
hvokm = np.array([- 0.038118, -0.068157, -0.031665])
# do straight RSW conversion
# get transofmration using tgt eci position and velocity
rrsw, vrsw, transmat = sc.rv2rsw(rtgteci, vtgteci)
# now transform hills from RSW to eci
rhrsw = transmat.T @ hrokm
vhrsw = transmat.T @ hvokm
# form int vector
rinteci = rhrsw + rtgteci
vinteci = vhrsw + vtgteci
print('rinteci RSW only %20.13f  %20.13f  %20.13f  %20.13f  %20.13f  %20.13f \n'
      % (rinteci[0], rinteci[1], rinteci[2], vinteci[0], vinteci[1], vinteci[2]))
# now try EQCm conversion
rinteci, vinteci = sc.hilleqcm2eci(rtgteci, vtgteci, hrokm, hvokm)

print('rinteci EQCM hold pos %20.13f  %20.13f  %20.13f  %20.13f  %20.13f  %20.13f \n'
      % (rinteci[0], rinteci[1], rinteci[2], vinteci[0], vinteci[1], vinteci[2]))
#stk ans
rintans = np.array([1e-12, 9495.039900814565, 5048.602269265176])
vintans = np.array([- 6.088186588973, 3.203827617391, 1.703505361584])
print('STK answer            %20.13f  %20.13f  %20.13f  %20.13f  %20.13f  %20.13f \n'
      % (rintans[0], rintans[1], rintans[2],
         vintans[0], vintans[1], vintans[2]))
rdiff = rinteci - rintans
print('diff from STK answer  %20.13f  %20.13f  %20.13f \n'
      % (rdiff[0], rdiff[1], rdiff[2]))
# plotorb(rtgteci, vtgteci, 16678.1385,'red')
# now try circular orbit at same rtgt altitude
p = smu.mag(rtgteci)
r, v = sc.coe2rv(p, 0.0, incl, omega, argp, nu, arglat, truelon, lonper)
rinteciC, vinteciC = sc.hilleqcm2eci(r, v, hrokm, hvokm)

print('rinteci EQ rtn cir %20.13f  %20.13f  %20.13f  %20.13f  %20.13f  %20.13f \n'
      % (rinteciC[0], rinteciC[1], rinteciC[2],
         vinteciC[0], vinteciC[1], vinteciC[2]))
# hold('on')
# plotorb(r, v, p,'green')
rinthill, vinthill = sc.eci2hilleqcm(rtgteci, vtgteci, rinteci, vinteci)
print('rinthill EQCM     %20.13f  %20.13f  %20.13f  %20.13f  %20.13f  %20.13f \n'
      % (rinthill[0], rinthill[1], rinthill[2],
         vinthill[0], vinthill[1], vinthill[2]))
rinthill, vinthill = sc.eci2hilleqcm(rtgteci, vtgteci, rinteciC, vinteciC)
print('rinthill EQCM cir %20.13f  %20.13f  %20.13f  %20.13f  %20.13f  %20.13f \n'
      % (rinthill[0], rinthill[1], rinthill[2],
         vinthill[0], vinthill[1], vinthill[2]))
# try NTW frame - hold velocity direction constant
# not sure this one is correct
rinteciz, vinteciz = sc.hilleqcm2eci(rtgteci, vtgteci, hrokm, hvokm)

print('rinteci EQCM ntw     %20.13f  %20.13f  %20.13f  %20.13f  %20.13f  %20.13f \n'
      % (rinteciz[0], rinteciz[1], rinteciz[2],
         vinteciz[0], vinteciz[1], vinteciz[2]))
# test kozai back and forth for tle
rad = 180.0 / np.pi
x2o3 = 2.0 / 3.0
mu = 398600.8

re = 6378.135

xke = 60.0 / np.sqrt(re * re * re / mu)

tumin = 1.0 / xke

j2 = 0.001082616
xpdotp = 1440.0 / (2.0 * np.pi)

no_kozai = 15.5911407

ecco = 0.000577
inclo = 51.619 * deg2rad
eccsq = ecco * ecco
omeosq = 1.0 - eccsq
rteosq = np.sqrt(omeosq)
cosio = np.cos(inclo)
cosio2 = cosio * cosio
# twoline2rv:
no_kozai = no_kozai / xpdotp

a1 = (no_kozai * tumin) ** (- 2.0 / 3.0)

# initl (inside init)
ak = (xke / no_kozai) ** x2o3

d1 = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq)
del1 = d1 / (ak * ak)
adel = ak * (1.0 - del1 * del1 - del1 * (1.0 / 3.0 + 134.0 * del1 * del1 / 81.0))
del_ = d1 / (adel * adel)
no_unkozai = no_kozai / (1.0 + del_)

ao = (xke / no_unkozai) ** x2o3

print('nokoz  %11.6f   unkoz  %11.6f  a  %11.6f er a1  %11.6f er \n'
      % (no_kozai, no_unkozai, ao, a1))
print('nokoz  %11.6f   unkoz  %11.6f  a  %11.6f km a1  %11.6f km \n'
      % (no_kozai, no_unkozai, ao * re, a1 * re))
# now go in reverse

# assume for 1st guess that no_kozai = no_unkozai
ak = (xke / no_unkozai) ** x2o3

d1 = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq)

del1 = d1 / (ak * ak)
adel = ak * (1.0 - del1 * del1 - del1 * (1.0 / 3.0 + 134.0 * del1 * del1 / 81.0))
del_ = d1 / (adel * adel)
no_kozai = no_unkozai * (1.0 + del_)
print('nokoz  %11.6f   unkoz  %11.6f  a  %11.6f km a1  %11.6f km \n'
      % (no_kozai * 1440 / (2 * np.pi), no_unkozai * 1440 / (2 * np.pi),
         ao * re, a1 * re))
