import numpy as np
import math
from space_constants import *
from space_constants import sethelp
# ------------------------------------------------------------------------------
#
#                           procedure mincomb
#
#  this procedure calculates the delta v's and the change in inclination
#    necessary for the minimum change in velocity when traveling between two
#    non-coplanar orbits.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    rinit       - initial position magnitude     er
#    rfinal      - final position magnitude       er
#    einit       - ecc of first orbit
#    e2          - ecc of trans orbit
#    efinal      - ecc of final orbit
#    nuinit      - true anomaly of first orbit    0 or pi rad
#    nufinal     - true anomaly of final orbit    0 or pi rad
#    iinit       - incl of the first orbit        rad
#    ifinal      - incl of the second orbit       rad
#
#  outputs       :
#    deltai1     - amount of incl chg req at a    rad
#    deltava     - change in velocity at point a  er / tu
#    deltavb     - change in velocity at point b  er / tu
#    dttu        - time of flight for the trans   tu
#    numiter     - number of iterations
#
#  locals        :
#    sme1        - mech energy of first orbit     er2 / tu
#    sme2        - mech energy of transfer orbit  er2 / tu
#    sme3        - mech energy of final orbit     er2 / tu
#    vinit       - velocity of first orbit at a   er / tu
#    vtransa     - velocity of trans orbit at a   er / tu
#    vtransb     - velocity of trans orbit at b   er / tu
#    vfinal      - velocity of final orbit at b   er / tu
#    ainit       - semimajor axis of first orbit  er
#    atrans      - semimajor axis of trans orbit  er
#    afinal      - semimajor axis of final orbit  er
#    e2          - eccentricity of second orbit
#
#  coupling      :
#    power       - raise a base to a power
#    asin      - arc sine routine
#
#  references    :
#    vallado       2007, 355, alg 42, table 6-3
#function [deltai,deltai1,deltava,deltavb,dttu ] = mincomb(rinit,rfinal,einit,efinal,nuinit,nufinal,iinit,ifinal);
# -----------------------------------------------------------------------------

def mincomb(rinit: float, rfinal: float, einit: float, e2: float,
            efinal: float, nuinit: float, nufinal: float, iinit: float,
            ifinal: float):
    a1 = (rinit * (1.0 + einit * math.cos(nuinit))) / (1.0 - einit*einit)
    a2 = 0.5 * (rinit + rfinal)
    a3 = (rfinal * (1.0 + efinal * math.cos(nufinal))) / (1.0 - efinal*efinal)
    sme1 = -1.0 / (2.0 * a1)
    sme2 = -1.0 / (2.0 * a2)
    sme3 = -1.0 / (2.0 * a3)

    vinit = math.sqrt(2.0 * ((1.0 / rinit) + sme1))
    v1t = math.sqrt(2.0 * ((1.0 / rinit) + sme2))
    vfinal = math.sqrt(2.0 * ((1.0 / rfinal) + sme3))
    v3t = math.sqrt(2.0 * ((1.0 / rfinal) + sme2))

    tdi = ifinal - iinit

    temp = (1.0/tdi) * math.atan(((rfinal/rinit)**1.5 - math.cos(tdi))
                                 / math.sin(tdi))
    deltai = temp * tdi
    deltava = math.sqrt(v1t*v1t + vinit*vinit - 2.0 * v1t * vinit *
                        math.cos(deltai))

    temp = (1.0/tdi) * math.atan(math.sin(tdi) /
                                 ((rfinal / rinit)**1.5 + math.cos(tdi)))
    deltai1 = tdi * (1.0 - temp)
    deltavb = math.sqrt(v3t*v3t + vfinal*vfinal - 2.0 * v3t * vfinal *
                        math.cos(deltai1))

    dttu = math.pi * math.sqrt(a2**3)

    if sethelp.iauhelp == 'y':
        dvold = abs(v1t-vinit) + math.sqrt(v3t*v3t + vfinal*vfinal
                                          - 2.0*v3t*vfinal*math.cos(tdi))
        print(f's = {temp:11.7}')
        print(f'rinit: {rinit:14.7} er, {rinit*re:14.7} km')
        print(f'rfinal: {rfinal:14.7} er, {rfinal*re:14.7} km')
        print(f'deltai: {deltai*rad2deg:13.7}, deltai1: {deltai1*rad2deg:13.7} er/tu')
        print(f'deltava: {deltava:13.7}, deltavb: {deltavb:13.7}')
        print(f'deltava: {deltava*velkmps:13.7}, deltavb: {deltavb*velkmps:13.7} km/s')
        print(f'{1000 * (deltava+deltavb)*velkmps:13.7} m/s')
        print(f'dv old way: {1000*dvold*velkmps:13.7} m/s')
        print(f'dttu: {dttu*tumin:13.7} min')

    deltainew = deltai
    deltai1 = 100
    numiter = 0

    while abs(deltainew - deltai1) > smalle6:
        deltai1= deltainew
        deltava= math.sqrt(v1t*v1t + vinit*vinit - 2.0 * v1t * vinit * math.cos(deltai1))
        deltavb= math.sqrt(v3t*v3t + vfinal*vfinal - 2.0 * v3t * vfinal * math.cos(tdi-deltai1))
        deltainew= math.asin((deltava * vfinal * v3t *
                              math.sin(tdi - deltai1)) / (vinit * v1t * deltavb))
        numiter = numiter + 1
    print(f'iter id {deltai1*rad2deg:14.6}; {numiter:3},'
          f'{deltava+deltavb*1000*velkmps}')

    return deltai, deltai1, deltava, deltavb, dttu