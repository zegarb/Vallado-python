import numpy as np
import os
from space_constants import *
import spacemath_utils as smu
import space_conversions as sc
import orbit_utils as obu


def dolam(kt, outfile, nrev, kepmov, rtgto, vtgto, rinto, direc, energy):
    outfile.write('xx  \n')
    outfile.write('xx \n')
    for i in range((nrev) * 100 , 501):
        dt = (i+1) * 60.0
        # make target moving...
        if kepmov == 'y':
            rtgt1, vtgt1, errork = obu.kepler(rtgto, vtgto, dt)
        else:
            errork = 'ok'
            rtgt1 = rtgto
            vtgt1 = vtgto

        # if i == nrev * 100:
        #     # check min energy condition and min time for that
        #     cosdeltanu = np.dot(rinto, rtgt1) / (smu.mag(rinto) * smu.mag(rtgt1))
        #     chord = np.sqrt(smu.mag(rinto) ** 2 + smu.mag(rtgt1) ** 2 - 2.0 * smu.mag(rinto) * smu.mag(rtgt1) * cosdeltanu)
        #     s = (smu.mag(rinto) + smu.mag(rtgt1) + chord) * 0.5
        #     amin = s / 2.0
        #     betam = 2.0 * np.arcsin(np.sqrt((s - chord) / s))
        #     alpham = np.pi
        #     if (direc == 's'):
        #         # this one works - gets min for min a on each case
        #         ttran = np.sqrt(amin ** 3 / mu) * (2.0 * nrev * np.pi + alpham - np.sin(alpham) - betam + np.sin(betam))
        #         #    tpar  = (s^1.5 - sin(trangle)*(s-chord)^1.5) * sqrt(2)/(3.0*sqrt(mu));
        #         #    tpar  =  sqrt(s^3/(8.0*mu)) * (pi - betam + sin(betam));
        #         tpar = (s ** 1.5 - (s - chord) ** 1.5) * np.sqrt(2) / (3.0 * np.sqrt(mu))
        #     else:
        #         ttran = np.sqrt(amin ** 3 / mu) * (2.0 * nrev * np.pi + alpham - np.sin(alpham) + betam - np.sin(betam))
        #         #    tpar  = (s^1.5 + sin(trangle)*(s-chord)^1.5) * sqrt(2)/(3.0*sqrt(mu));
        #         #    tpar  =  sqrt(s^3/(8.0*mu)) * (pi + betam - sin(betam));
        #         # this one works - chk when rest is good if olong/short are +-
        #         tpar = (s ** 1.5 + (s - chord) ** 1.5) * np.sqrt(2) / (3.0 * np.sqrt(mu))
        #     beta = 2.0 * np.arcsin(np.sqrt((s - chord) / (s)))
        #     tmin = np.sqrt(amin ** 3 / mu) * ((2.0 * nrev + 1.0) * np.pi - beta + np.sin(beta))
        #     outfile.write('dnu %11.7f mins c %11.7f s %11.7f a %11.7f be %11.7f tranmin %11.7f tmin %11.7f  tpar %11.7f\n'
        #           % (np.arccos(cosdeltanu) * 180 / np.pi, chord, s, amin,
        #              betam, ttran, tmin, tpar))




        outfile.write(f'{kt = }, {direc = }, {energy = }, {nrev = }, {i = }, {dt = }\n')
        tmin, _, _, _ = obu.lambertmintof(rinto, rtgt1, direc, nrev, outfile)
        if dt < tmin:
            outfile.write(f'dt is less than minimum tof \n')
        else:
            tbi, _ = obu.lambgettbiu(rinto, rtgt1, nrev)
            vtrans1, vtrans2, errorl = \
                obu.lambertu(rinto, vtgto, rtgt1, direc, energy, nrev,
                                dt, tbi, outfile)
            outfile.write(f'Univ:\n{vtrans1 = }, \n{vtrans2 = }\n')
            if errorl == 'ok' and errork == 'ok':
                p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper \
                    = sc.rv2coe(rinto, vtrans1); # of trans orbit
                dv1 = smu.mag(vinto - vtrans1)
                dv2 = smu.mag(vtrans2 - vtgt1)
                outfile.write(' %11.5f %11.5f %11.5f %11.5f %11.5f \n'
                                % (a, ecc, dv1, dv2, dv1+dv2))
            else:
                outfile.write('  0  0 %s \n' % errorl)

            vtrans1, vtrans2, errorl = obu.lambertb(rinto, vinto, rtgt1, direc, energy, nrev, dt )
            outfile.write(f'Batt:\n{vtrans1 = }, \n{vtrans2 = }\n')
            if errorl == 'ok' and errork == 'ok':
                p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper \
                    = sc.rv2coe(rinto, vtrans1); # of trans orbit
                dv1 = smu.mag(vinto - vtrans1)
                dv2 = smu.mag(vtrans2 - vtgt1)
                outfile.write(' %11.5f %11.5f %11.5f %11.5f %11.5f \n\n'
                                % (a, ecc, dv1, dv2, dv1+dv2) )
            else:
                outfile.write('  0  0 %s \n' % errorl)

st = 'y'
outpath = os.path.join(os.path.dirname(__file__), 'testoutput', 'testlam.out')
outfile = open(outpath, 'wt')
for kt in range(1, 8):
    # ---- use for fixed r and r locations, and fig 7-12, 16, and 20 ---- }
    #      readln( infile, ro(1), ro(2), ro(3), rtgt(1), rtgt(2), rtgt(3), dt, dmy, direc );
    #      fprintf( kt:3,'fig 7-12 ro rtgt ', ro(4):11:7, rtgt(4):11:7 );
    kepmov = 'n'
    if kt == 1:
        rinto = np.array([1.04, 0.0, 0.0]) * re
        rtgto = np.array([-1.04, 0.2, 0.0]) * re
    if kt == 2:
        rinto = np.array([2.0, 0.0, 0.0]) * re
        rtgto = np.array([-2.0, 0.2, 0.0]) * re
    if kt == 3:
        rinto = np.array([4.0, 0.0, 0.0]) * re
        rtgto = np.array([-4.0, 0.2, 0.0]) * re
    if kt == 4:
        rinto = np.array([6.6107, 0.0, 0.0]) * re
        rtgto = np.array([-6.6107, 0.2, 0.0]) * re
    # make up a velocity vector (circular orbit) for dv calcs...
    if (kt > 0) and (kt < 5):
        vinto = np.array([0.0, math.sqrt(mu / smu.mag(rinto)), 0.0])
        vtgto = np.array([0.0, math.sqrt(mu / smu.mag(rtgto)), 0.0])
    if kt == 5:
        rinto = np.array([-6518.1083,-2403.8479,-22.1722])
        vinto = np.array([2.604057,-7.105717,-0.263218])
        rtgto = np.array([6697.4756, 1794.5831, 0.0])
        vtgto = np.array([-1.962372, 7.323674, 0.0])
        kepmov = 'n'
    if kt == 6:
        rinto = np.array([5328.7862, 4436.1273, 101.472])
        vinto = np.array([-4.864779, 5.816486, 0.240163])
        rtgto = np.array([6697.4756, 1794.5831, 0.0])
        vtgto = np.array([-1.962372, 7.323674, 0.0])
        kepmov = 'y'
    if kt == 7:
        trangle = 90.0 * deg2rad
        rinto = np.array([1.0, 0.0, 0.0]) * re
        vinto = np.array([2.604057,-7.105717,- 0.263218])
        rtgto = np.array([1.0 * math.cos(trangle), 1.0 * math.sin(trangle), 0.0]) * re
        vtgto = np.array([-1.962372, 7.323674, 0.0])
        kepmov = 'n'

    direcs = ['s', 'l']
    energys = ['l', 'h']
    for nrev in range(0, 3):
        for direc in direcs:
            if nrev == 0:
                dolam(kt, outfile, nrev, kepmov, rtgto, vtgto, rinto, direc, 'l')
                print(f'{kt = }, {direc = }, {nrev = }')
            else:
                for energy in energys:
                    dolam(kt, outfile, nrev, kepmov, rtgto, vtgto, rinto, direc, energy)
                    print(f'{kt = }, {direc = }, {nrev = }, {energy = }')

    # direc = 's'
    # nrev = 0
    # dolam(outfile, nrev, kepmov, rtgto, vtgto, rinto, direc)
    # direc = 'l'
    # nrev = 0
    # dolam(outfile, nrev, kepmov, rtgto, vtgto, rinto, direc)
    # direc = 's'
    # nrev = 1
    # # kt = 1: X falls below 0 in lambertu when nrev is set
    # # to 1 for this test case. Creates a negative square root in an
    # # atan function error. Can't figure out what is going wrong. - mjc
    # dolam(outfile, nrev, kepmov, rtgto, vtgto, rinto, direc)
    # direc = 'l'
    # nrev = 1
    # dolam(outfile, nrev, kepmov, rtgto, vtgto, rinto, direc)
    # direc = 's'
    # nrev = 2
    # dolam(outfile, nrev, kepmov, rtgto, vtgto, rinto, direc)
    # direc = 'l'
    # nrev = 2
    # dolam(outfile, nrev, kepmov, rtgto, vtgto, rinto, direc)

outfile.close()
print('done with fig 7-13 data ')