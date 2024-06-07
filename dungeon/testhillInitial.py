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
print('------------------------------------------------ initial accuracy ----------------------------- \n' % ())
# now test the ability to convert eci - hills and back
constmath
constastro
casenumo = 1

casetest = 7

ang_step = 1e-08
# case where you read in the int and tgt ephemerides
ropt = 'm'

outfilehill = open(strcat('d:/STKFiles Educational Files/Hills/thillc',int2str(casenumo),'n',int2str(casetest),'.out'),'wt')

# -----------------------------------------------------------------------------------------------
# -------------------------------- do initial accuracy checks fwd, back, ------------------------
# -----------------------------------------------------------------------------------------------
# --- this one also sets the case for the rest of the program here...
for ktr in np.arange(1,10+1).reshape(-1):
    # set last one to be the case of interest
    if ktr < 10:
        casenum = ktr
    else:
        casenum = casenumo
    if 1 == casenum:
        rtgteci[1] = (6378.137 + 500.0)
        rtgteci[2] = 0.0
        rtgteci[3] = 0.0
        magrt = mag(rtgteci)
        vtgteci[1] = 0.0
        vtgteci[2] = np.sqrt(mu / magrt)
        vtgteci[3] = 0.0
        circ = 'y'
        rtgteci = np.transpose(rtgteci)
        vtgteci = np.transpose(vtgteci)
        #  previous doesn't work in STK92 - fixed in 9.3 for smaller inclination
        a = 6378.137 + 500.0
        ecc = 0.0
        p = a * (1.0 - ecc * ecc)
        rtgteci,vtgteci = coe2rv(p,0.0,0.001 / rad,0.0,0.0,0.0,0.0,0.0,0.0)
        circ = 'y'
        print('rtgt = [%20.13f %20.13f %20.13f]; \n vtgt = [%20.13f %20.13f %20.13f]; \n' % (rtgteci(1),rtgteci(2),rtgteci(3),vtgteci(1),vtgteci(2),vtgteci(3)))
    else:
        if 2 == casenum:
            a = 26500.0
            ecc = 0.0
            p = a * (1.0 - ecc * ecc)
            rtgteci,vtgteci = coe2rv(p,0.0,0.001 / rad,0.0,0.0,0.0,0.0,0.0,0.0)
            circ = 'y'
        else:
            if 3 == casenum:
                a = 26500.0
                ecc = 0.0
                p = a * (1.0 - ecc * ecc)
                rtgteci,vtgteci = coe2rv(p,0.0,55.0 / rad,0.0,0.0,0.0,0.0,0.0,0.0)
                circ = 'y'
            else:
                if 4 == casenum:
                    a = 42164.0
                    ecc = 0.0
                    p = a * (1.0 - ecc * ecc)
                    rtgteci,vtgteci = coe2rv(p,0.0,0.001 / rad,0.0,0.0,0.0,0.0,0.0,0.0)
                    circ = 'y'
                else:
                    if 5 == casenum:
                        rtgteci = np.transpose(np.array([- 605.7904308,- 5870.230407,3493.052004]))
                        vtgteci = np.transpose(np.array([- 1.568251615,- 3.702348353,- 6.479484915]))
                        circ = 'n'
                    else:
                        if 6 == casenum:
                            rtgteci = np.transpose(np.array([- 761.24075519,22699.40899449,13644.176032]))
                            vtgteci = np.transpose(np.array([- 2.630715586,1.396719096,- 2.491891]))
                            circ = 'n'
                        else:
                            if 7 == casenum:
                                rtgteci = np.transpose(np.array([- 40588.150362,- 11462.167028,27.147649]))
                                vtgteci = np.transpose(np.array([0.834787457,- 2.958305691,- 0.001173016]))
                                circ = 'n'
                            else:
                                if 8 == casenum:
                                    rtgteci = np.transpose(np.array([- 5551.898646,- 2563.049696,3257.756165]))
                                    vtgteci = np.transpose(np.array([2.149073,- 7.539457,- 2.185709]))
                                    circ = 'n'
                                else:
                                    if 9 == casenum:
                                        rtgteci = np.transpose(np.array([9668.14551571,6999.07240705,4041.43303674]))
                                        vtgteci = np.transpose(np.array([- 3.65292306,0.64966519,- 5.82123516]))
                                        circ = 'n'
                                    else:
                                        if 10 == casenum:
                                            a = 8164.7188
                                            ecc = 0.08
                                            p = a * (1.0 - ecc * ecc)
                                            rtgteci,vtgteci = coe2rv(p,ecc,32.861 / rad,0.0,0.0,0.0,0.0,0.0,0.0)
                                            circ = 'n'
    x = 10.0
    y = 10.0
    z = 10.0
    xd = 0.01
    yd = 0.01
    zd = 0.01
    hro1[1] = x / 1000.0
    hro1[2] = y / 1000.0
    hro1[3] = z / 1000.0
    hrokm = np.transpose(hro1)
    hvo1[1] = xd / 1000.0
    hvo1[2] = yd / 1000.0
    hvo1[3] = zd / 1000.0
    hvokm = np.transpose(hvo1)
    rintecix,vintecix = EQCM_to_ECI_RTN_sal(rtgteci,vtgteci,hrokm,hvokm)
    rhillx,vhillx = ECI_to_EQCM_RTN_sal(rtgteci,vtgteci,rintecix,vintecix,outfilehill)
    dr = rhillx - hrokm
    dv = vhillx - hvokm
    mdr = mag(dr * 1000000)
    mdv = mag(dv * 1000000)
    print('dr m         %20.13f %20.13f %20.13f %15.8f mm %20.13f %20.13f %20.13f %15.8f mm/s \n' % (dr(1) * 1000,dr(2) * 1000,dr(3) * 1000,mdr,dv(1) * 1000,dv(2) * 1000,dv(3) * 1000,mdv))
    p,a,ecc,incl,omega,argpx1,nux1,m,arglat,truelon,lonper = rv2coe(rtgteci,vtgteci)
    print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n' % (p,a,ecc,incl * rad,omega * rad,argpx1 * rad,nux1 * rad,m * rad,arglat * rad,truelon * rad,lonper * rad))
    p,a,ecc,incl,omega,argpx1,nux1,m,arglat,truelon,lonper = rv2coe(rintecix,vintecix)
    print('coes %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n' % (p,a,ecc,incl * rad,omega * rad,argpx1 * rad,nux1 * rad,m * rad,arglat * rad,truelon * rad,lonper * rad))


pause
# -----------------------------------------------------------------------------------------------
# ------------------- check various positions to determine if the veocity is correct ------------
# -----------------------------------------------------------------------------------------------
for ktr in np.arange(1,6+1).reshape(-1):
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
    hro1[1] = x / 1000.0
    hro1[2] = y / 1000.0
    hro1[3] = z / 1000.0
    hrokm = np.transpose(hro1)
    hvo1[1] = xd / 1000.0
    hvo1[2] = yd / 1000.0
    hvo1[3] = zd / 1000.0
    hvokm = np.transpose(hvo1)
    rintecisal,vintecisal = EQCM_to_ECI_RTN_sal(rtgteci,vtgteci,hrokm,hvokm)
    print('hillsin       %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n' % (hrokm(1) * 1000,hrokm(2) * 1000,hrokm(3) * 1000,hvokm(1) * 1000,hvokm(2) * 1000,hvokm(3) * 1000))
    mmagrti = mag(rtgteci)
    mmagvti = mag(vtgteci)
    print('rtgteci        %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f %15.8f %15.8f \n' % (rtgteci(1),rtgteci(2),rtgteci(3),vtgteci(1),vtgteci(2),vtgteci(3),mmagrti,mmagvti))
    mmagri = mag(rintecisal)
    mmagvi = mag(vintecisal)
    print('rinteci        %20.13f %20.13f %20.13f   %20.13f %20.13f %20.13f  %15.8f %15.8f \n' % (rintecisal(1),rintecisal(2),rintecisal(3),vintecisal(1),vintecisal(2),vintecisal(3),mmagri,mmagvi))

pause
# -----------------------------------------------------------------------------------------------
# -------------------------- perturb each one, or all of the vector components ------------------
# -----------------------------------------------------------------------------------------------
print('\n\n------------------------------- perturb each one ---------------------------------------- \n' % ())
# initialize
x = 0.0

y = 0.0

z = 0.0

xd = 0.0

yd = 0.0

zd = 0.0

fid = 2
for ktr in np.arange(1,6+1).reshape(-1):
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
    hro1[1] = x / 1000.0
    hro1[2] = y / 1000.0
    hro1[3] = z / 1000.0
    hrokm = np.transpose(hro1)
    hvo1[1] = xd / 1000.0
    hvo1[2] = yd / 1000.0
    hvo1[3] = zd / 1000.0
    hvokm = np.transpose(hvo1)
    # reset this!! so m for hills call!!
    hro[1] = x
    hro[2] = y
    hro[3] = z
    hvo[1] = xd
    hvo[2] = yd
    hvo[3] = zd
    rintecisal,vintecisal = EQCM_to_ECI_RTN_sal(rtgteci,vtgteci,hrokm,hvokm)
    rhill2,vhill2 = ECI_to_EQCM_RTN_sal(rtgteci,vtgteci,rintecisal,vintecisal,outfilehill)
    rintecisal1,vintecisal1 = EQCM_to_ECI_RTN_sal(rtgteci,vtgteci,rhill2,vhill2)
    rhillsal,vhillsal = ECI_to_EQCM_RTN_sal(rtgteci,vtgteci,rintecisal1,vintecisal1,outfilehill)
    rhillsal1,vhillsal1 = ECI_to_EQCM_RTN_sal(rtgteci,vtgteci,rintecisal,vintecisal,outfilehill)
    # ------------ set the vectors for future use in program
    rinteci[1] = rintecisal(1)
    rinteci[2] = rintecisal(2)
    rinteci[3] = rintecisal(3)
    vinteci[1] = vintecisal(1)
    vinteci[2] = vintecisal(2)
    vinteci[3] = vintecisal(3)
    # set initial values to compare and make sure STK doesn't have an error
    rintecio = rinteci
    vintecio = vinteci
    rtgtecio = rtgteci
    vtgtecio = vtgteci
    # --------------------- write out various transformation methods --- just at ktr = 8
    if (ktr <= 8):
        print('\nhills in ell   %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n' % (hro(1),hro(2),hro(3),hvo(1),hvo(2),hvo(3)))
        print('\nhills in       %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n' % (hrokm(1) * 1000,hrokm(2) * 1000,hrokm(3) * 1000,hvokm(1) * 1000,hvokm(2) * 1000,hvokm(3) * 1000))
        mmagrti = mag(rtgteci)
        mmagvti = mag(vtgteci)
        print('rtgteci        %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n' % (rtgteci(1),rtgteci(2),rtgteci(3),mmagrti,vtgteci(1),vtgteci(2),vtgteci(3),mmagvti))
        print('rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n' % (rhill2(1) * 1000,rhill2(2) * 1000,rhill2(3) * 1000,vhill2(1) * 1000,vhill2(2) * 1000,vhill2(3) * 1000))
        mmagrti = mag(rinteci)
        mmagvti = mag(vinteci)
        print('rinteci        %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n' % (rinteci(1),rinteci(2),rinteci(3),mmagrti,vinteci(1),vinteci(2),vinteci(3),mmagvti))
        mmagrti = mag(rintecisal)
        mmagvti = mag(vintecisal)
        print('rintecisal hil %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n' % (rintecisal(1),rintecisal(2),rintecisal(3),mmagrti,vintecisal(1),vintecisal(2),vintecisal(3),mmagvti))
        mmagrti = mag(rintecisal1)
        mmagvti = mag(vintecisal1)
        print('rintecisal1 eq %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n' % (rintecisal1(1),rintecisal1(2),rintecisal1(3),mmagrti,vintecisal1(1),vintecisal1(2),vintecisal1(3),mmagvti))
        print('hillssal hill  %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n' % (rhillsal(1) * 1000,rhillsal(2) * 1000,rhillsal(3) * 1000,vhillsal(1) * 1000,vhillsal(2) * 1000,vhillsal(3) * 1000))
        print('hillssal1 eq   %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n' % (rhillsal1(1) * 1000,rhillsal1(2) * 1000,rhillsal1(3) * 1000,vhillsal1(1) * 1000,vhillsal1(2) * 1000,vhillsal1(3) * 1000))
        print('==== sals codes ====\n' % ())
        #---            [rintecisal,vintecisal]   = EQCM_to_ECI_NTW_sal( rtgteci*1000,vtgteci*1000,hro',hvo' );      # do in m~!!!find int pos eqcm
#---            [rhill2, vhill2]          = ECI_to_EQCM_NTW_sal( rtgteci*1000,vtgteci*1000,rintecisal,vintecisal ); # find how much this would be for true hills
        rintecisal,vintecisal = EQCM_to_ECI_RTN_sal(rtgteci,vtgteci,np.transpose(hro),np.transpose(hvo))
        rhill2,vhill2 = ECI_to_EQCM_RTN_sal(rtgteci,vtgteci,rintecisal,vintecisal,outfilehill)
        mmagrti = mag(rintecisal)
        mmagvti = mag(vintecisal)
        print('rintecisal hil %20.13f %20.13f %20.13f %15.8f %20.13f %20.13f %20.13f %15.8f \n' % (rintecisal(1) * 0.001,rintecisal(2) * 0.001,rintecisal(3) * 0.001,mmagrti * 0.001,vintecisal(1) * 0.001,vintecisal(2) * 0.001,vintecisal(3) * 0.001,mmagvti * 0.001))
        print('rhills2 in     %20.13f %20.13f %20.13f %20.13f %20.13f %20.13f \n' % (rhill2(1),rhill2(2),rhill2(3),vhill2(1),vhill2(2),vhill2(3)))

fclose(outfilehill)
