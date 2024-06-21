#
#
# dolam - code insertion for testlam program
#

import numpy as np
from space_constants import *
import orbit_utils as obu
import spacemath_utils as smu

muin = mu

# ---- do the short way multiple revolution cases ---- }
fid.write('xx  \n' % ())
print('xx 501    psinew         dt       x      a          e \n' % ())
for i in range((nrev) * 100, 500+1):
    dt = i * 60.0
    # make target moving...
    if kepmov == 'y':
        rtgt1, vtgt1, errork = obu.kepler(rtgto, vtgto, dt, 0)
    else:
        errork = '      ok'
        rtgt1 = rtgto
        vtgt1 = vtgto
    if i == nrev * 100:
        # check min energy condition and min time for that
        cosdeltanu = np.dot(rinto, rtgt1) / (smu.mag(rinto) * smu.mag(rtgt1))
        chord = np.sqrt(smu.mag(rinto) ** 2 + smu.mag(rtgt1) ** 2 - 2.0 * smu.mag(rinto) * smu.mag(rtgt1) * cosdeltanu)
        s = (smu.mag(rinto) + smu.mag(rtgt1) + chord) * 0.5
        amin = s / 2.0
        betam = 2.0 * np.arcsin(np.sqrt((s - chord) / s))
        alpham = np.pi
        if (direc == 's'):
            # this one works - gets min for min a on each case
            ttran = np.sqrt(amin ** 3 / mu) * (2.0 * nrev * np.pi + alpham - np.sin(alpham) - betam + np.sin(betam))
            #    tpar  = (s^1.5 - sin(trangle)*(s-chord)^1.5) * sqrt(2)/(3.0*sqrt(mu));
            #    tpar  =  sqrt(s^3/(8.0*mu)) * (pi - betam + sin(betam));
            tpar = (s ** 1.5 - (s - chord) ** 1.5) * np.sqrt(2) / (3.0 * np.sqrt(mu))
        else:
            ttran = np.sqrt(amin ** 3 / mu) * (2.0 * nrev * np.pi + alpham - np.sin(alpham) + betam - np.sin(betam))
            #    tpar  = (s^1.5 + sin(trangle)*(s-chord)^1.5) * sqrt(2)/(3.0*sqrt(mu));
            #    tpar  =  sqrt(s^3/(8.0*mu)) * (pi + betam - sin(betam));
            # this one works - chk when rest is good if olong/short are +-
            tpar = (s ** 1.5 + (s - chord) ** 1.5) * np.sqrt(2) / (3.0 * np.sqrt(mu))
        beta = 2.0 * np.arcsin(np.sqrt((s - chord) / (s)))
        tmin = np.sqrt(amin ** 3 / mu) * ((2.0 * nrev + 1.0) * np.pi - beta + np.sin(beta))
        print('dnu %11.7f mins c %11.7f s %11.7f a %11.7f be %11.7f tranmin %11.7f tmin %11.7f  tpar %11.7f\n' % (np.arccos(cosdeltanu) * 180 / np.pi, chord, s, amin, betam, ttran, tmin, tpar))
    vtrans1, vtrans2, errorl = obu.lambertu(rinto, rtgt1, direc, nrev, dt, fid)
    #          [vtrans1, vtrans2, errorl] = lambertb( rinto, rtgt1, direc, overrev, dt );
    #           if strcmp(errorl, '      ok') ==0 && strcmp(errork, '      ok') ==0
    #               [p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper ] = rv2coe (rinto, vtrans1, muin); # of trans orbit
    #               dv1 = smu.mag(vinto - vtrans1);
    #               dv2 = smu.mag(vtrans2 - vtgt1);
    #               fprintf( fid,' #11.5f #11.5f #11.5f #11.5f #11.5f \n', a, ecc, dv1, dv2, dv1+dv2 );
    #           else
    #               fprintf( fid,'  0  0 #s \n', errorl );
    #           end;
