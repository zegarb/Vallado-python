import math
import numpy as np
from numpy.polynomial import Polynomial
from pprint import pprint as pp
from space_constants import *
import spacemath_utils as smu
import spacetime_utils as stu
import orbit_utils as obu
from space_constants import sethelp as sh
from arclength_ellipse import arclength_ellipse
from elliptic12 import elliptic12
from inverselliptic2 import inverselliptic2

# ------------------------------------------------------------------------------
#
#                           function ap2kp
#
#  this function converts ap to kp using cubic splines
#
#  author        : david vallado                  719-573-2600   4 aug  2005
#
#  revisions
#
#  inputs          description                      range / units
#    apin        - geomagnetic planetary amplitude  gamma (10^-9 Tesla)
#
#  outputs       :
#    kpout       - geomagnetic planetary index
#
#  locals        :
#                -
#
#  coupling      :
#
#  references    :
#    vallado       2004, 899-901
#
# kpout = ap2kp(apin)
# ------------------------------------------------------------------------------

def ap2kp(apin: float):
    """this function converts ap to kp using cubic splines

    Parameters
    ----------
    apin : float
        geomagnetic planetary amplitude

    Returns
    -------
    kpout : float
        geomagnetic planetary index
    """

    ap = np.array([- 1e-05, - 0.001, 0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22,
                   27, 32, 39, 48, 56, 67, 80, 94, 111, 132, 154, 179, 207,
                   236, 300, 400, 900])
    kp = np.array([- 0.66666667, - 0.33333, 0.0, 0.33333, 0.66667, 1.0,
                   1.33333, 1.66667, 2.0, 2.33333, 2.66667, 3.0, 3.33333,
                   3.66667, 4.0, 4.33333,4.66667, 5.0, 5.33333, 5.66667, 6.0,
                   6.33333, 6.66667, 7.0, 7.33333, 7.66667, 8.0, 8.33333,
                   8.66667, 9.0, 9.33333])
    kpout = 0.0
    # find starting point in files based on input
    i = 0
    while ((i < 30) and (apin > ap[i])):
        i = i + 1

    if (i > 2):
        # -------- assign function points ---------
        p1 = kp[i - 2]
        p2 = kp[i - 1]
        p3 = kp[i]
        p4 = kp[i + 1]
        kac0, kac1, kac2, kac3 = smu.cubicspl(p1, p2, p3, p4)
        p1 = ap[i - 2]
        p2 = ap[i - 1]
        p3 = ap[i]
        p4 = ap[i + 1]
        ac0, ac1, ac2, ac3 = smu.cubicspl(p1, p2, p3, p4)

        # recover the original function values
        # use the normalized time first, but at an arbitrary interval
        # r1r, r1i, r2r, r2i, r3r, r3i = smu.cubic(kac3, kac2, kac1, kac0 - apin, 'R')
        # r = np.zeros(3)
        # r[0] = r1r
        # r[1] = r2r
        # r[2] = r3r
        r = Polynomial([ac0 - apin, ac1, ac2, ac3]).roots().real
        if ((r[0] >= 0.0) and (r[0] <= 1.001)):
            apval = r[0]
        elif ((r[1] >= 0.0) and (r[1] <= 1.001)):
            apval = r[1]
        elif ((r[2] >= 0.0) and (r[2] <= 1.001)):
            apval = r[2]
        else:
            apval = 0.0
            print('error in root %11.7f %3i %11.7f %11.7f %11.7f \n'
                   % (apin, i, r[0], r[1], r[2]))
        kpout = (kac3 * apval ** 3 + kac2 * apval ** 2 + kac1 * apval + kac0)


    return kpout

# ------------------------------------------------------------------------------
#
#                           function kp2ap
#
#  this function converts kp to ap using cubic splines
#
#  author        : david vallado                  719-573-2600   4 aug  2005
#
#  revisions
#
#  inputs          description                      range / units
#    kpin        - geomagnetic planetary index      gamma (10^-9 Tesla)
#
#  outputs       :
#    apout       - geomagnetic planetary amplitude
#
#  locals        :
#                -
#
#  coupling      :
#
#  references    :
#    vallado       2004, 899-901
#
# apout = kp2ap(kpin)
# ------------------------------------------------------------------------------

def kp2ap(kpin: float):
    """this function converts kp to ap using cubic splines

    Parameters
    ----------
    kpin : float
        geomagnetic planetary index

    Returns
    -------
    apout : float
        geomagnetic planetary amplitude
    """

    ap = np.array([- 1e-05, - 0.001, 0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22,
                   27, 32, 39, 48, 56, 67, 80, 94, 111, 132, 154, 179, 207,
                   236, 300, 400, 900])
    kp = np.array([- 0.66666667, - 0.33333, 0.0, 0.33333, 0.66667, 1.0,
                   1.33333, 1.66667, 2.0, 2.33333, 2.66667, 3.0, 3.33333,
                   3.66667, 4.0, 4.33333,4.66667, 5.0, 5.33333, 5.66667, 6.0,
                   6.33333, 6.66667, 7.0, 7.33333, 7.66667, 8.0, 8.33333,
                   8.66667, 9.0, 9.33333])
    apout = 0.0
    # find starting point in files based on input
    i = 0
    while ((i < 30) and (kpin > kp[i])):
        i = i + 1

    if (i > 2):
        # -------- assign function points ---------
        p1 = kp[i - 2]
        p2 = kp[i - 1]
        p3 = kp[i]
        p4 = kp[i + 1]
        kac0, kac1, kac2, kac3 = smu.cubicspl(p1, p2, p3, p4)
        p1 = ap[i - 2]
        p2 = ap[i - 1]
        p3 = ap[i]
        p4 = ap[i + 1]
        ac0, ac1, ac2, ac3 = smu.cubicspl(p1, p2, p3, p4)

        # recover the original function values
        # use the normalized time first, but at an arbitrary interval
        # r1r, r1i, r2r, r2i, r3r, r3i = smu.cubic(kac3, kac2, kac1, kac0 - apin, 'R')
        # r = np.zeros(3)
        # r[0] = r1r
        # r[1] = r2r
        # r[2] = r3r
        r = Polynomial([kac0 - kpin, kac1, kac2, kac3]).roots().real
        if ((r[0] >= 0.0) and (r[0] <= 1.001)):
            apval = r[0]
        elif ((r[1] >= 0.0) and (r[1] <= 1.001)):
            apval = r[1]
        elif ((r[2] >= 0.0) and (r[2] <= 1.001)):
            apval = r[2]
        else:
            apval = 0.0
            print('error in root %11.7f %3i %11.7f %11.7f %11.7f \n'
                   % (kpin, i, r[0], r[1], r[2]))
        apout = (ac3 * apval ** 3 + ac2 * apval ** 2 + ac1 * apval + ac0)


    return apout

# ----------------------------------------------------------------------------
#
#                           function rv2eq
#
#  this function transforms a position and velocity vector into the equinoctial
#  coordinate system.
#
#  author        : david vallado                  719-573-2600    7 jun 2002
#
#  revisions
#    vallado     - fix special orbit types (ee)                   5 sep 2002
#    vallado     - add constant file use                         29 jun 2003
#
#  inputs          description                               range / units
#    reci        - eci position vector                       km
#    veci        - eci velocity vector                       km/s
#
#  outputs       :
#    n           - mean motion                               rad
#    a           - semi major axis                           km
#    af          - component of ecc vector
#    ag          - component of ecc vector
#    chi         - component of node vector in eqw
#    psi         - component of node vector in eqw
#    meanlon     - mean longitude                            rad
#    truelon     - true longitude                            rad
#
#  locals        :
#    none        -
#
#  coupling      :
#    none        -
#
#  references    :
#    vallado       2013, 108
#    chobotov            30
#
# [a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr] = rv2eq (r, v)
# ----------------------------------------------------------------------------

def rv2eq(reci: np.ndarray, veci: np.ndarray):
    """this function transforms a position and velocity vector into the equinoctial
    coordinate system.

    Parameters
    ----------
    reci : ndarray
        eci position vector
    veci : ndarray
        eci velocity vector

    Returns
    -------
    n: float
        mean motion: rad
    a: float
        semi major axis: km
    af
        component of ecc vector
    ag
        component of ecc vector
    chi
        component of node vector in eqw
    psi
        component of node vector in eqw
    meanlon: float
        mean longitude: rad
    truelon: float
        true longitude: rad
    """

    # -------- convert to classical elements ----------------------
    _, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = \
        rv2coe(reci, veci)
    # -------- setup retrograde factor ----------------------------
    fr = 1.0
    # ---------- set this so it is -1 only for orbits near 180 deg !! ---------
    if np.abs(incl - np.pi) < 0.0001:
        fr = - 1.0

    coe = True
    if (coe == True):
        if (ecc < small):
            # ----------------  circular equatorial  ------------------
            if (incl < small) or (np.abs(incl - np.pi) < small):
                argp = 0.0
                omega = 0.0
                nu = truelon
                m = truelon
            else:
                # --------------  circular inclined  ------------------
                argp = 0.0
                nu = arglat
                m = arglat
        else:
            # ---------------  elliptical equatorial  -----------------
            if ((incl < small) or (np.abs(incl - np.pi) < small)):
                argp = lonper
                omega = 0.0
        n = np.sqrt(mu / (a * a * a))
        af = ecc * np.cos(fr * omega + argp)
        ag = ecc * np.sin(fr * omega + argp)
        if (fr > 0):
            chi = np.tan(incl * 0.5) * np.sin(omega)
            psi = np.tan(incl * 0.5) * np.cos(omega)
        else:
            chi = np.cot(incl * 0.5) * np.sin(omega)
            psi = np.cot(incl * 0.5) * np.cos(omega)
        meanlonM = fr * omega + argp + m
        meanlonM = np.fmod(meanlonM + twopi, twopi)
        meanlonNu = fr * omega + argp + nu
        meanlonNu = np.fmod(meanlonNu + twopi, twopi)
        eccanom, nu = smu.newtonm(ecc, m)
        #         print('rv2eq F #11.7f  L #11.7f \n', (fr*omega + argp + eccanom) * rad2deg, meanlonNu * rad2deg)
#         print('rv2eq F #11.7f  L #11.7f \n', (fr*omega + argp + eccanom) * rad2deg, (fr*omega + argp + nu) * rad2deg)
    else:
        magr = smu.mag(reci)
        magv = smu.mag(veci)
        a = 1.0 / (2.0 / magr - magv ** 2 / mu)
        n = np.sqrt(mu / (a * a * a))
        wvec = np.cross(reci, veci) / smu.mag(np.cross(reci, veci))
        chi = wvec[0] / (1.0 + fr * wvec[2])
        psi = - wvec[1] / (1.0 + fr * wvec[2])
        p0 = 1.0 / (1.0 + chi ** 2 + psi ** 2)
        fe = p0 * (1.0 - chi ** 2 + psi ** 2)
        fq = p0 * 2.0 * chi * psi
        fw = p0 * - 2.0 * fr * chi
        fvec = np.array([fe, fq, fw])
        ge = p0 * 2.0 * fr * chi * psi
        gq = p0 * fr * (1.0 + chi ** 2 - psi ** 2)
        gw = p0 * 2.0 * psi
        gvec = np.array([ge, gq, gw])
        we = wvec[0]
        wq = wvec[1]
        ww = wvec[2]
        we = p0 * 2.0 * chi
        wq = p0 * - 2.0 * psi
        ww = p0 * fr * (1.0 - chi ** 2 - psi ** 2)
        wvec1 = np.array([we, wq, ww])
        # alt formulation - same as other, EXCEPT for retrograde (fr-1) so
# do not use!!
#         fe = 1.0 - we**2/(1.0 + ww)
#         fq = -we*wq / (1.0 + ww)
#         fw = -we
#         fvec = [fe, fq, fw]
#         gvec = np.cross(wvec, fvec)
        evec = -reci / magr + np.cross(veci, np.cross(reci, veci)) / mu
        ag = np.dot(evec, gvec)
        af = np.dot(evec, fvec)
        X = np.dot(reci, fvec)
        Y = np.dot(reci, gvec)
        b = 1.0 / (1.0 + np.sqrt(1.0 - af ** 2 - ag ** 2))
        sinF = ag + (((1.0 - ag ** 2 * b) * Y - ag * af * b * X)
                     / (a * np.sqrt(1.0 - ag ** 2 - af ** 2)))
        cosF = af + (((1.0 - af ** 2 * b) * X - ag * af * b * Y)
                     / (a * np.sqrt(1.0 - ag ** 2 - af ** 2)))
        F = math.atan2(sinF, cosF)
        #        print('rv2eq fe #11.7f #11.7f #11.7f ge #11.7f  #11.7f #11.7f \n X #11.7f Y #11.7f b #11.7f sF #11.7f cF #11.7f \n', fe, fq, fw, ge, gq, gw, X, Y, b, sinF, cosF)
#        print('F = 316.20515  L = 13.61834  M = 288.88793 \n')
        sinZeta = ag / np.sqrt(af ** 2 + ag ** 2)
        cosZeta = af / np.sqrt(af ** 2 + ag ** 2)
        zeta = math.atan2(sinZeta, cosZeta)
        meanlonM = F + ag * np.cos(F) - af * np.sin(F)
        if meanlonM < 0.0:
            meanlonM = twopi + meanlonM
        Eccanom = F - zeta
        M = Eccanom - ecc * np.sin(Eccanom)
        if M < 0.0:
            M = 2.0 * np.pi + M
        #M = meanlonM - zeta  # same
        sinL = (((1.0 - af ** 2 * b) * np.sin(F) + ag * af * b * np.cos(F) - ag)
                / (1.0 - ag * np.sin(F) - af * np.cos(F)))
        cosL = (((1.0 - ag ** 2 * b) * np.cos(F) + ag * af * b * np.sin(F) - af)
                / (1.0 - ag * np.sin(F) - af * np.cos(F)))
        L = math.atan2(sinL, cosL)
        meanlonNu = L
        if meanlonNu < 0.0:
            meanlonNu = twopi + meanlonNu
        nu = L - zeta
        #       print('rv2eq F #11.7f  L #11.7f  #11.7f #11.7f #11.7f #11.7f  \n', F * rad2deg, L * rad2deg, X, Y, zeta * rad2deg, M * rad2deg)

    return a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr

# ------------------------------------------------------------------------------
#
#                           function eq2rv
#
#  this function finds the classical orbital elements given the equinoctial
#    elements.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#    vallado     - fix elliptical equatorial orbits case         19 oct 2002
#    vallado     - add constant file use                         29 jun 2003
#
#  inputs          description                    range / units
#    a           - semimajor axis                 km
#    af          - component of ecc vector
#    ag          - component of ecc vector
#    chi         - component of node vector in eqw
#    psi         - component of node vector in eqw
#    meanlonM    - mean longitude                 rad
#    fr          - retrograde factor              +1 or -1
#
#  outputs       :
#    reci        - position vector                km
#    veci        - velocity vector                km/s
#
#  locals        :
#    n           - mean motion                    rad
#    temp        - temporary variable
#    p           - semilatus rectum               km
#    ecc         - eccentricity
#    incl        - inclination                    0.0  to pi rad
#    omega       - longitude of ascending node    0.0  to 2pi rad
#    argp        - argument of perigee            0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#    m           - mean anomaly                   0.0  to 2pi rad
#    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
#    truelon     - true longitude            (ce) 0.0  to 2pi rad
#    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
#
#  coupling      :
#
#  references    :
#    vallado 2013:108
#
# [r, v] = eq2rv(a, af, ag, chi, psi, meanlonM, fr)
# ------------------------------------------------------------------------------

def eq2rv(a: float, af: float, ag: float, chi: float, psi: float,
          meanlonM: float, fr: float):
    """this function finds the classical orbital elements given the equinoctial
    elements.

    Parameters
    ----------
    a : float
        semimajor axis
    af : float
        component of ecc vector
    ag : float
        component of ecc vector
    chi : float
        component of node vector in eqw
    psi : float
        component of node vector in eqw
    meanlonM : float
        mean longitude
    fr : float
        retrograde factor (+1 for regular orbits, -1 for retrograde)

    Returns
    -------
    reci: ndarray
        eci position vector
    veci: ndarray
        eci velocity vector
    """
    # -------------------------  implementation   -----------------
    arglat = undefined
    lonper = undefined
    truelon = undefined
    coe = True

    if (coe == True):
        # ---- if n is input ----
        #a = (mu/n^2)^(1.0/3.0)
        ecc = np.sqrt(af ** 2 + ag ** 2)
        p = a * (1.0 - ecc * ecc)
        incl = (np.pi * ((1.0 - fr) * 0.5) + 2.0 * fr
                * np.arctan(np.sqrt(chi ** 2 + psi ** 2)))
        omega = math.atan2(chi, psi)
        argp = math.atan2(ag, af) - fr * omega
        if (ecc < small):
            # ----------------  circular equatorial  ------------------
            if (incl < small) or (np.abs(incl - np.pi) < small):
                #argp = 0.0
                #omega = 0.0
                truelon = omega
            else:
                # --------------  circular inclined  ------------------
                #argp = 0.0
                arglat = argp
        else:
            # ---------------  elliptical equatorial  -----------------
            if ((incl < small) or (np.abs(incl - np.pi) < small)):
                #argp = lonper
                lonper = omega
        m = meanlonM - fr * omega - argp
        m = np.fmod(m + twopi, twopi)
        e0, nu = smu.newtonm(ecc, m)
        # ----------  fix for elliptical equatorial orbits ------------
        if (ecc < small):
            # ----------------  circular equatorial  ------------------
            if (incl < small) or (np.abs(incl - np.pi) < small):
                argp = undefined
                omega = undefined
                truelon = nu
            else:
                # --------------  circular inclined  ------------------
                argp = undefined
                arglat = nu - fr * omega
            nu = undefined
        else:
            # ---------------  elliptical equatorial  -----------------
            if ((incl < small) or (np.abs(incl - np.pi) < small)):
                lonper = argp
                argp = undefined
                omega = undefined
        # -------- now convert back to position and velocity vectors
        reci, veci = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper)
    else:
        p0 = 1.0 / (1.0 + chi ** 2 + psi ** 2)
        fe = p0 * (1.0 - chi ** 2 + psi ** 2)
        fq = p0 * 2.0 * chi * psi
        fw = p0 * - 2.0 * fr * chi
        fvec = np.array([fe, fq, fw])
        ge = p0 * 2.0 * fr * chi * psi
        gq = p0 * fr * (1.0 + chi ** 2 - psi ** 2)
        gw = p0 * 2.0 * psi
        gvec = np.array([ge, gq, gw])
        we = p0 * 2.0 * chi
        wq = p0 * - 2.0 * psi
        ww = p0 * fr * (1.0 - chi ** 2 - psi ** 2)
        wvec = np.array([we, wq, ww])
        F0 = meanlonM
        numiter = 25
        ktr = 1
        F1 = (F0 - (F0 + ag * np.cos(F0) - af * np.sin(F0) - meanlonM)
              / (1.0 - ag * np.sin(F0) - af * np.cos(F0)))
        while ((np.abs(F1 - F0) > small) and (ktr <= numiter)):

            ktr = ktr + 1
            F0 = F1
            F1 = (F0 - (F0 + ag * np.cos(F0) - af * np.sin(F0) - meanlonM)
                  / (1.0 - ag * np.sin(F0) - af * np.cos(F0)))
            #            fprintf(1, 'iters #7i  #11.7f  #11.7f \n', ktr, F0, F1)

        F = F1
        F = np.fmod(F + twopi, twopi)
        n = np.sqrt(mu / (a * a * a))
        b = 1.0 / (1.0 + np.sqrt(1.0 - af ** 2 - ag ** 2))
        sinL = (((1.0 - af ** 2 * b) * np.sin(F) + ag * af * b * np.cos(F) - ag)
                / (1.0 - ag * np.sin(F) - af * np.cos(F)))
        cosL = (((1.0 - ag ** 2 * b) * np.cos(F) + ag * af * b * np.sin(F) - af)
                / (1.0 - ag * np.sin(F) - af * np.cos(F)))
        L = math.atan2(sinL, cosL)
        meanlonNu = L
        meanlonNu = np.fmod(meanlonNu + twopi, twopi)
        rr = (a * (1.0 - ag ** 2 - af ** 2)
              / (1.0 + ag * np.sin(L) + af * np.cos(L)))
        rr1 = a * (1.0 - ag * np.sin(F) - af * np.cos(F))
        # coordinates in equinoctial space
        X = a * ((1.0 - ag ** 2 * b) * np.cos(F) + af * ag * b * np.sin(F) - af)
        Y = a * ((1.0 - af ** 2 * b) * np.sin(F) + af * ag * b * np.cos(F) - ag)
        XD = - n * a * (ag + sinL) / (np.sqrt(1.0 - af ** 2 - ag ** 2))
        YD = n * a * (af + cosL) / (np.sqrt(1.0 - af ** 2 - ag ** 2))
        # alt forms all are the same now
        #XD = n*a^2/magr * (af*ag*b*cos(F) - (1.0 - ag^2*b)*sin(F))
        #YD = n*a^2/magr * ((1.0 - af^2*b)*cos(F) - af*ag*b*sin(F))
        #        fprintf(1, 'eq2rv fe #11.7f #11.7f #11.7f ge #11.7f  #11.7f #11.7f \n X #11.7f Y #11.7f b #11.7f sF #11.7f cF #11.7f \n', fe, fq, fw, ge, gq, gw, X, Y, b, sin(F), cos(F))
        #         r = a*(1.0-af^2-ag^2) / (1.0 + ag*sinL + af*cosL)
        #         r = a*(1.0 - ag*sin(F) - af*cos(F))
        #         r = magr
        #         X = r*cosL
        #         Y = r*sinL
        reci = X * fvec + Y * gvec
        veci = XD * fvec + YD * gvec
        #   fprintf(1, 'eq2rv F #11.7f  L #11.7f \n', F * rad2deg, L * rad2deg)

    # test for eqw axes
    #     r = [6524.834000000,  6862.875000000,  6448.296000000]
    #     v = [4.9013270000,    5.5337560000,   -1.9763410000]

    #     incl = 87.8691262/rad
    #     raan = 227.89826    /rad


    #     [outvec] = rot3(r, raan)
    #     [outvec1] = rot1(outvec, incl)
    #     [ans] = rot3(outvec1, -raan)
    #    ans = 1.0e+04 * 1.113447759019948   0.269748810750719  0.000000005640130 correct

    return reci, veci




# ----------------------------------------------------------------------------
#
#                           function covct2eq
#
#  this function transforms a six by six covariance matrix expressed in
#    cartesian vectors into one expressed in equinoctial elements.
#
#  author        : david vallado                  719-573-2600   14 jul 2002
#
#  revisions
#    vallado     - major update                                  26 aug 2015
#
#  inputs          description                    range / units
#    cartcov     - 6x6 cartesian covariance matrix
#    cartstate   - 6x1 cartesian orbit state      (x y z vx vy vz)
#    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
#    fr          - retrograde factor               +1, -1
#
#  outputs       :
#    eqcov       - 6x6 equinoctial covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    r           - matrix of partial derivatives
#    a           - semimajor axis                 km
#    ecc         - eccentricity
#    incl        - inclination                    0.0  to pi rad
#    omaga       - longitude of ascending node    0.0  to 2pi rad
#    argp        - argument of perigee            0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#    m           - mean anomaly                   0.0  to 2pi rad
#    e0          - eccentric anomaly              0.0  to 2pi rad
#    tau         - time from perigee passage
#    n           - mean motion                    rad
#    af          - component of ecc vector
#    ag          - component of ecc vector
#    chi         - component of node vector in eqw
#    psi         - component of node vector in eqw
#    meanlon     - mean longitude                 rad
#
#  coupling      :
#
#  references    :
#    Vallado and Alfano 2015
#
# [eqcov, tm] = covct2eq (cartcov, cartstate, anom, fr)
# ----------------------------------------------------------------------------

def covct2eq (cartcov, cartstate, anom, fr):

    # -------- parse the input vectors into cartesian and classical components
    rx = cartstate[0, 0] * 1000.0
    ry = cartstate[1, 0] * 1000.0
    rz = cartstate[2, 0] * 1000.0
    vx = cartstate[3, 0] * 1000.0
    vy = cartstate[4, 0] * 1000.0
    vz = cartstate[5, 0] * 1000.0
    reci = np.array([rx, ry, rz])
    veci = np.array([vx, vy, vz])
    magr = smu.mag(reci) # m
    magv = smu.mag(veci) # m

    a = 1.0 / (2.0/magr - magv**2/mum)
    n = np.sqrt(mum / a**3)

    hx = ry*vz - rz*vy
    hy = -rx*vz + rz*vx
    hz = rx*vy - ry*vx
    h_vec = np.array([hx, hy, hz]).T
    #h_vec = cross(reci, veci)  #same
    w_vec = h_vec/smu.mag(h_vec)
    chi = w_vec[0]/(1.0 + fr * w_vec[2])
    psi = -w_vec[1]/(1.0 + fr * w_vec[2])

    # components of equinoctial system
    p0 = 1.0/(1.0 + chi**2 + psi**2)
    fe = p0 * (1.0 - chi**2 + psi**2)
    fq = p0 * 2.0 * chi * psi
    fw = p0 * -2.0 * fr * chi
    f_vec = np.array([fe, fq, fw]).T
    ge = p0 * 2.0 * fr* chi * psi
    gq = p0 * fr * (1.0 + chi**2 - psi**2)
    gw = p0 * 2.0 * psi
    g_vec = np.array([ge, gq, gw]).T
    we = w_vec[0]
    wq = w_vec[1]
    ww = w_vec[2]
    w_vec = np.array([we, wq, ww]).T

    r_dot_v = np.dot(reci, veci)
    p1 = magv*magv - mum/magr
    ecc_x = (p1*rx - r_dot_v*vx)/mum
    ecc_y = (p1*ry - r_dot_v*vy)/mum
    ecc_z = (p1*rz - r_dot_v*vz)/mum
    ecc_vec = np.array([ecc_x, ecc_y, ecc_z]).T

    af = np.dot(ecc_vec, f_vec)
    ag = np.dot(ecc_vec, g_vec)

    X = np.dot(reci, f_vec)
    Y = np.dot(reci, g_vec)

    b = 1.0 / (1.0 + np.sqrt(1.0 - af**2 - ag**2))
    p0 = 1.0 / (a*np.sqrt(1.0 - af**2 - ag**2))
    sinF = ag + p0*((1.0 - ag**2*b)*Y - ag*af*b*X)
    cosF = af + p0*((1.0 - af**2*b)*X - ag*af*b*Y)
    F = math.atan2(sinF, cosF)
    if F < 0.0:
        F = F + 2.0 * np.pi

    XD = n*a**2/magr * (af*ag*b*np.cos(F) - (1.0 - ag**2*b)*np.sin(F))
    YD = n*a**2/magr * ((1.0 - af**2*b)*np.cos(F) - af*ag*b*np.sin(F))

    A = np.sqrt(mum*a)
    B = np.sqrt(1.0 - ag**2 - af**2)
    C = 1.0 + chi**2 + psi**2

    partXDaf = a*XD*YD / (A*B) - A/(magr**3)*(a*ag*X/(1 + B) + X*Y/B)
    partYDaf = -a*XD**2 / (A*B) - A/(magr**3)*(a*ag*Y/(1 + B) - X**2/B)
    partXDag = a*YD**2 / (A*B) + A/(magr**3)*(a*af*X/(1 + B) - Y**2/B)
    partYDag = -a*XD*YD / (A*B) + A/(magr**3)*(a*af*Y/(1 + B) + X*Y/B)

    print('na %11.7f %11.7f  \n'%(n, a))
    print('XY %11.7f %11.7f %11.7f %11.7f\n'%(X, Y, XD, YD))
    print('ABC %11.7f %11.7f %11.7f \n'%(A, B, C))
    print('part %11.7f %11.7f %11.7f %11.7f \n'%(partXDaf, partYDaf, partXDag,
                                                 partYDag))

    # ---------------- calculate matrix elements ------------------
    # ---- partials of a wrt (rx ry rz vx vy vz)
    if (anom == 'truea') or (anom == 'meana'):
        p0 = 2.0*a**2 / magr**3
        p1 = 2.0 / (n**2*a)
    elif (anom == 'truen') or (anom == 'meann'):
        p0 = -3.0*n*a/magr**3
        p1 = -3.0 / (n*a**2)
    tm = np.zeros((6, 6))

    tm[0, 0] = p0*rx
    tm[0, 1] = p0*ry
    tm[0, 2] = p0*rz
    tm[0, 3] = p1*vx
    tm[0, 4] = p1*vy
    tm[0, 5] = p1*vz

    # ---- partials of v wrt ag
    tm34 = partXDag*fe + partYDag*ge
    tm35 = partXDag*fq + partYDag*gq
    tm36 = partXDag*fw + partYDag*gw
    # ---- partials of af wrt (rx ry rz vx vy vz)
    p0 = 1.0 / mum
    tm[1, 0] = (-a*b*af*B*rx/(magr**3) - (ag*(chi*XD - psi*fr*YD)*we)/(A*B)
               + (B/A)*tm34)
    tm[1, 1] = (-a*b*af*B*ry/(magr**3) - (ag*(chi*XD - psi*fr*YD)*wq)/(A*B)
               + (B/A)*tm35)
    tm[1, 2] = (-a*b*af*B*rz/(magr**3) - (ag*(chi*XD - psi*fr*YD)*ww)/(A*B)
               + (B/A)*tm36)
    tm[1, 3] = (p0*((2.0*X*YD - XD*Y)*ge - Y*YD*fe)
               - (ag*(psi*fr*Y - chi*X)*we) / (A*B))
    tm[1, 4] = (p0*((2.0*X*YD - XD*Y)*gq - Y*YD*fq)
               - (ag*(psi*fr*Y - chi*X)*wq) / (A*B))
    tm[1, 5] = (p0*((2.0*X*YD - XD*Y)*gw - Y*YD*fw)
               - (ag*(psi*fr*Y - chi*X)*ww) / (A*B))

    # ---- partials of v wrt af
    tm24 = partXDaf*fe + partYDaf*ge
    tm25 = partXDaf*fq + partYDaf*gq
    tm26 = partXDaf*fw + partYDaf*gw
    # ---- partials of ag wrt (rx ry rz vx vy vz)
    p0 = 1.0 / mum
    tm[2, 0] = (-a*b*ag*B*rx/(magr**3) + (af*(chi*XD - psi*fr*YD)*we)/(A*B)
               - (B/A)*tm24)
    tm[2, 1] = (-a*b*ag*B*ry/(magr**3) + (af*(chi*XD - psi*fr*YD)*wq)/(A*B)
               - (B/A)*tm25)
    tm[2, 2] = (-a*b*ag*B*rz/(magr**3) + (af*(chi*XD - psi*fr*YD)*ww)/(A*B)
               - (B/A)*tm26)
    tm[2, 3] = (p0*((2.0*XD*Y - X*YD)*fe - X*XD*ge)
               + (af*(psi*fr*Y - chi*X)*we) / (A*B))
    tm[2, 4] = (p0*((2.0*XD*Y - X*YD)*fq - X*XD*gq)
               + (af*(psi*fr*Y - chi*X)*wq) / (A*B))
    tm[2, 5] = (p0*((2.0*XD*Y - X*YD)*fw - X*XD*gw)
               + (af*(psi*fr*Y - chi*X)*ww) / (A*B))

    # ---- partials of chi wrt (rx ry rz vx vy vz)
    tm[3, 0] = -C*YD*we / (2.0*A*B)
    tm[3, 1] = -C*YD*wq / (2.0*A*B)
    tm[3, 2] = -C*YD*ww / (2.0*A*B)
    tm[3, 3] = C*Y*we / (2.0*A*B)
    tm[3, 4] = C*Y*wq / (2.0*A*B)
    tm[3, 5] = C*Y*ww / (2.0*A*B)

    # ---- partials of psi wrt (rx ry rz vx vy vz)
    tm[4, 0] = -fr*C*XD*we / (2.0*A*B)
    tm[4, 1] = -fr*C*XD*wq / (2.0*A*B)
    tm[4, 2] = -fr*C*XD*ww / (2.0*A*B)
    tm[4, 3] = fr*C*X*we / (2.0*A*B)
    tm[4, 4] = fr*C*X*wq / (2.0*A*B)
    tm[4, 5] = fr*C*X*ww / (2.0*A*B)

    if (anom == 'truea') or (anom == 'truen'):
        # not ready yet
        #p0 = -sign(argp)/np.sqrt(1-np.cos(argp)*np.cos(argp))
        #p1 = 1.0 / mum
        #tm[5, 0] = p0*(p1*()/(n*ecc) - ()/n*()/(n**2*ecc) - tm(ecc/ry*()) + fr*-vz*nodey/n**2 + ...
        #
        #tm[5, 1] = p0*()
        #tm[5, 2] = p0*()
        #tm[5, 3] = p0*()
        #tm[5, 4] = p0*()
        #tm[5, 5] = p0*()
        tm[5, 0] = 0
        tm[5, 1] = 0
        tm[5, 2] = 0
        tm[5, 3] = 0
        tm[5, 4] = 0
        tm[5, 5] = 0
    elif (anom == 'meana') or (anom == 'meann'):
        # ---- partials of meanlon wrt (rx ry rz vx vy vz)
        tm[5, 0] = (-vx/A + (chi*XD - psi*fr*YD)*we/(A*B)
                   - (b*B/A)*(ag*tm34 + af*tm24))
        tm[5, 1] = (-vy/A + (chi*XD - psi*fr*YD)*wq/(A*B)
                   - (b*B/A)*(ag*tm35 + af*tm25))
        tm[5, 2] = (-vz/A + (chi*XD - psi*fr*YD)*ww/(A*B)
                   - (b*B/A)*(ag*tm36 + af*tm26))
        tm[5, 3] = (-2.0*rx/A + (af*tm[2, 3] - ag*tm[1, 3])/(1.0 + B)
                   + (fr*psi*Y - chi*X)*we/A)
        tm[5, 4] = (-2.0*ry/A + (af*tm[2, 4] - ag*tm[1, 4])/(1.0 + B)
                   + (fr*psi*Y - chi*X)*wq/A)
        tm[5, 5] = (-2.0*rz/A + (af*tm[2, 5] - ag*tm[1, 5])/(1.0 + B)
                   + (fr*psi*Y - chi*X)*ww/A)

    # ---------- calculate the output covariance matrix -----------
    eqcov = tm*cartcov*tm.T
    return eqcov, tm

# ----------------------------------------------------------------------------
#
#                           function covct2fl
#
#  this function transforms a six by six covariance matrix expressed in cartesian elements
#    into one expressed in flight parameters
#
#  author        : david vallado                  719-573-2600   21 jun 2002
#
#  revisions
#    vallado     - major update                                  26 aug 2015
#
#  inputs          description                           range / units
#    cartcov     - 6x6 cartesian covariance matrix
#    cartstate   - 6x1 cartesian orbit state             (x y z vx vy vz)
#    anom        - anomaly                               'latlon', 'radec'
#    ttt         - julian centuries of tt                centuries
#    jdut1       - julian date of ut1                    days from 4713 bc
#    lod         - excess length of day                  sec
#    xp          - polar motion coefficient              arc sec
#    yp          - polar motion coefficient              arc sec
#    terms       - number of terms for ast calculation   0, 2
#    ddpsi       - delta psi correction to gcrf          rad
#    ddeps       - delta eps correction to gcrf          rad
#
#  outputs       :
#    flcov       - 6x6 flight covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    r           - matrix of partial derivatives
#    x, y, z       - components of position vector  km
#    vx, vy, vz    - components of position vector  km/s
#    magr        - eci position vector magnitude  km
#    magv        - eci velocity vector magnitude  km/sec
#    d           - r dot v
#    h           - angular momentum vector
#    hx, hy, hz    - components of angular momentum vector
#    hcrossrx, y, z- components of h cross r vector
#    p1, p2       - denominator terms for the partials
#
#  coupling      :
#    ecef2eci    - convert eci vectors to ecef
#
#  references    :
#    Vallado and Alfano 2015
#
#   [flcov, tm] = covct2fl(cartcov, cartstate, anom, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ----------------------------------------------------------------------------

def covct2fl(cartcov, cartstate, anom, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps):

    # -------- parse the input vectors into cartesian components
    rx = cartstate[0, 0] * 1000  # keep all in m, m/s
    ry = cartstate[1, 0] * 1000  # this is eci always
    rz = cartstate[2, 0] * 1000
    vx = cartstate[3, 0] * 1000
    vy = cartstate[4, 0] * 1000
    vz = cartstate[5, 0] * 1000

    if (anom == 'latlon'):
        # -------- convert r to eci
        reci = np.array([rx, ry, rz])/1000
        veci = np.array([vx, vy, vz])/1000
        aeci = np.zeros(3)
        # now find ecef coordinates
        [recef, vecef, aecef] = eci2ecef(reci, veci, aeci, ttt, jdut1, lod, xp, yp,
                                       terms, ddpsi, ddeps)
        recef = recef * 1000  # m
        vecef = vecef * 1000
        rxf = recef[0]
        ryf = recef[1]
        rzf = recef[2]
# [rxf; ryf; rzf]
#    elif (anom, 'radec')

    # -------- calculate common quantities
    magr = np.sqrt(rx**2 + ry**2 + rz**2)   # in m
    magv = np.sqrt(vx**2 + vy**2 + vz**2)
    rdotv = rx*vx + ry*vy + rz*vz

    h = np.sqrt((rx*vy-ry*vx)**2 + (rz*vx-rx*vz)**2 + (ry*vz-rz*vy)**2)
#          hx = ry*vz - rz*vy  #
#          hy = rz*vx - rx*vz
#          hz = rx*vy - ry*vx
#          hcrossrx = (rz*hy - ry*hz)
#          hcrossry = (rx*hz - rz*hx)
#          hcrossrz = (ry*hx - rx*hy)

    tm = np.zeros((6, 6))
    # ---------------- calculate matrix elements ------------------
    if (anom == 'latlon'):
        # partial of lon wrt (x y z vx vy vz)
        p0 = 1.0 / (rxf**2 + ryf**2)
        tm[0, 0] = -p0*ryf
        tm[0, 1] = p0*rxf
        tm[0, 2] = 0.0
        tm[0, 3] = 0.0
        tm[0, 4] = 0.0
        tm[0, 5] = 0.0

        # partial of latgc wrt (x y z vx vy vz)
        p0 = 1.0 / (magr**2*np.sqrt(rxf**2 + ryf**2))
        tm[1, 0] = -p0*(rxf*rzf)
        tm[1, 1] = -p0*(ryf*rzf)
        tm[1, 2] = np.sqrt(rxf**2 + ryf**2) / magr**2
        tm[1, 3] = 0.0
        tm[1, 4] = 0.0
        tm[1, 5] = 0.0
    elif (anom, 'radec'):
        # partial of lon wrt (x y z vx vy vz)
        p0 = 1.0 / (rx**2 + ry**2)
        tm[0, 0] = -p0*ry
        tm[0, 1] = p0*rx
        tm[0, 2] = 0.0
        tm[0, 3] = 0.0
        tm[0, 4] = 0.0
        tm[0, 5] = 0.0

        # partial of latgc wrt (x y z vx vy vz)
        p0 = 1.0 / (magr**2*np.sqrt(rx**2 + ry**2))
        tm[1, 0] = -p0*(rx*rz)
        tm[1, 1] = -p0*(ry*rz)
        tm[1, 2] = np.sqrt(rx**2 + ry**2) / magr**2
        tm[1, 3] = 0.0
        tm[1, 4] = 0.0
        tm[1, 5] = 0.0

    # partial of fpa wrt (x y z vx vy vz)
    rdot = rdotv / magr  # (r dot v) / r
#        p1 = -1.0 / (magr*np.sqrt(magv**2 - rdot**2))
#         tm[2, 0] = p1*(rdot*rx/magr - vx)
#         tm[2, 1] = p1*(rdot*ry/magr - vy)
#         tm[2, 2] = p1*(rdot*rz/magr - vz)
#         tm[2, 3] = p1*(rdot*magr*vx/magv**2 - rx)
#         tm[2, 4] = p1*(rdot*magr*vy/magv**2 - ry)
#         tm[2, 5] = p1*(rdot*magr*vz/magv**2 - rz)
    # Sal from mathcad matches previous with - sign on previous
    p0 = 1.0 / (magr**2*h)
    p1 = 1.0 / (magv**2*h)
    tm[2, 0] = p0*(vx*(ry**2 + rz**2) - rx*(ry*vy + rz*vz))
    tm[2, 1] = p0*(vy*(rx**2 + rz**2) - ry*(rx*vx + rz*vz))
    tm[2, 2] = p0*(vz*(rx**2 + ry**2) - rz*(rx*vx + ry*vy))
    tm[2, 3] = p1*(rx*(vy**2 + vz**2) - vx*(ry*vy + rz*vz))
    tm[2, 4] = p1*(ry*(vx**2 + vz**2) - vy*(rx*vx + rz*vz))
    tm[2, 5] = p1*(rz*(vx**2 + vy**2) - vz*(rx*vx + ry*vy))

    # partial of az wrt (x y z vx vy vz)
#         p2 = 1.0 / ((magv**2 - rdot**2) * (rx**2 + ry**2))
#         tm[3, 0] = p2*(vy*(magr*vz - rz*rdot) - (rx*vy - ry*vx) * (rx*vz - rz*vx + rx*rz*rdot/magr) * (1.0 / magr))
#         tm[3, 1] = p2*(-vx*(magr*vz - rz*rdot) + (rx*vy - ry*vx) * (ry*vz - rz*vy + ry*rz*rdot/magr) * (1.0 / magr))
#         p2 = 1.0 / (magr**2 * (magv**2 - rdot**2))
#         tm[3, 2] = p2 * rdot * (rx*vy - ry*vx)
#         p2 = 1.0 / (magr * (magv**2 - rdot**2))
#         tm[3, 3] = -p2 * (ry*vz - rz*vy)
#         tm[3, 4] = p2 * (rx*vz - rz*vx)
#         tm[3, 5] = -p2 * (rx*vy - ry*vx)
    # sal from mathcad
    p2 = 1.0 / ((magv**2 - rdot**2)*(rx**2 + ry**2))
    k1 = np.sqrt(rx**2 + ry**2 + rz**2)*(rx*vy - ry*vx)
    k2 = ry*(ry*vz - rz*vy) + rx*(rx*vz - rz*vx)
    tm[3, 0] = p2 * (vy*(magr*vz - rz*rdot) -
                    (rx*vy - ry*vx)/magr * (rx*vz - rz*vx + rx*ry*rdot/magr))
    p2 = 1.0 / (magr*(k1**2 + k2**2))
    tm[3, 0] = p2 * (k1*magr*(rz*vx - 2*rx*vz) +
                    k2*(-ry*vx*rx + vy*rx**2 + vy*magr**2))
    tm[3, 1] = p2 * (k1*magr*(rz*vy - 2*ry*vz) +
                    k2*(rx*vy*ry - vx*ry**2 - vx*magr**2))
    p2 = k1 / (magr**2*(k1**2 + k2**2))
    tm[3, 2] = p2 * (k2*rz + (rx*vx + ry*vy)*magr**2)
    p2 = 1.0 / (k1**2 + k2**2)
    tm[3, 3] = p2 * (k1*rx*rz - k2*ry*magr)
    tm[3, 4] = p2 * (k1*ry*rz + k2*rx*magr)
    tm[3, 5] = -p2 * (k1*(rx**2 + ry**2))

    # partial of r wrt (x y z vx vy vz)
    p0 = 1.0 / magr
    tm[4, 0] = p0*rx
    tm[4, 1] = p0*ry
    tm[4, 2] = p0*rz
    tm[4, 3] = 0.0
    tm[4, 4] = 0.0
    tm[4, 5] = 0.0

    # partial of v wrt (x y z vx vy vz)
    p0 = 1.0 / magv
    tm[5, 0] = 0.0
    tm[5, 1] = 0.0
    tm[5, 2] = 0.0
    tm[5, 3] = p0*vx
    tm[5, 4] = p0*vy
    tm[5, 5] = p0*vz

    # ---------- calculate the output covariance matrix -----------
    flcov = tm*cartcov*tm.T
    return flcov, tm

# ----------------------------------------------------------------------------
#
#                           function covcl2eq
#
#  this function transforms a six by six covariance matrix expressed in
#    classical elements into one expressed in equinoctial elements.
#
#  author        : david vallado                  719-573-2600   14 jul 2002
#
#  revisions
#    vallado     - major update                                  26 aug 2015
#
#  inputs          description                    range / units
#    classcov    - 6x6 classical covariance matrix
#    classstate  - 6x1 classical orbit state      (a e i O w nu/m)
#    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
#    fr          - retrograde factor               +1, -1
#
#  outputs       :
#    eqcov       - 6x6 equinoctial covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    a           - semimajor axis                 km
#    ecc         - eccentricity
#    incl        - inclination                    0.0  to pi rad
#    omaga       - longitude of ascending node    0.0  to 2pi rad
#    argp        - argument of perigee            0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#    m           - mean anomaly                   0.0  to 2pi rad
#    mum         - gravitational paramater        m**3/s**2 NOTE Meters!
#
#  coupling      :
#
#  references    :
#    Vallado and Alfano 2015
#
# [eqcov, tm] = covcl2eq (classcov, classstate, anom, fr)
# ----------------------------------------------------------------------------

def covcl2eq (classcov, classstate, anom, fr):

    # --------- determine which set of variables is in use ---------
    # -------- parse the orbit state
    a = classstate[0]  # in m
    n = np.sqrt(mum/a**3)
    ecc = classstate[1]
    incl = classstate[2]
    omega = classstate[3]
    argp = classstate[4]
    if (anom == 'meana') or (anom == 'meann'):
        m = classstate[5]
    elif (anom == 'truea') or (anom == 'truen'):
        nu = classstate[5]

    tm = np.zeros((6, 6))
    # ---- partials of a/n wrt (a ecc incl node argp nu/M)
    if (anom == 'truea') or (anom == 'meana'):
        tm[0, 0] = 1.0
    elif (anom == 'truen') or (anom == 'meann'):
        #tm[0, 0] = 1.0
        tm[0, 0] = -(3.0*np.sqrt(mum/a**3)) / (2.0*a)  # if class = a, equin = n

    # ---------------- calculate matrix elements ------------------
    tm[0, 1] = 0.0
    tm[0, 2] = 0.0
    tm[0, 3] = 0.0
    tm[0, 4] = 0.0
    tm[0, 5] = 0.0

    # ---- partials of af wrt (a ecc incl node argp nu/M)
    tm[1, 0] = 0.0
    tm[1, 1] = np.cos(fr*omega + argp)
    tm[1, 2] = 0.0
    tm[1, 3] = -ecc*fr*np.sin(fr*omega + argp)
    tm[1, 4] = -ecc*np.sin(fr*omega + argp)
    tm[1, 5] = 0.0

    # ---- partials of ag wrt (a ecc incl node argp nu/M)
    tm[2, 0] = 0.0
    tm[2, 1] = np.sin(fr*omega + argp)
    tm[2, 2] = 0.0
    tm[2, 3] = ecc*fr*np.cos(fr*omega + argp)
    tm[2, 4] = ecc*np.cos(fr*omega + argp)
    tm[2, 5] = 0.0

    # ---- partials of chi wrt (a ecc incl node argp nu/M)
    tm[3, 0] = 0.0
    tm[3, 1] = 0.0
    tm[3, 2] = (np.sin(omega) * (0.5*np.tan(incl*0.5)*np.tan(incl*0.5) + 0.5)
                * fr * np.tan(incl*0.5)**(fr-1))
    tm[3, 3] = np.tan(incl*0.5)**fr * np.cos(omega)
    tm[3, 4] = 0.0
    tm[3, 5] = 0.0

    # ---- partials of psi wrt (a ecc incl node argp nu/M)
    tm[4, 0] = 0.0
    tm[4, 1] = 0.0
    tm[4, 2] = (np.cos(omega) * (0.5*np.tan(incl*0.5)*np.tan(incl*0.5) + 0.5)
                * fr * np.tan(incl*0.5)**(fr-1))
    tm[4, 3] = -np.tan(incl*0.5)**fr * np.sin(omega)
    tm[4, 4] = 0.0
    tm[4, 5] = 0.0

    # ---- partials of meanlonM/meanlonNu wrt (a ecc incl node argp nu/M)
    if (anom == 'truea') or (anom == 'truen'):
        #[e0, nu] = smu.newtonm (ecc, anomaly)
        tm[5, 0] = 0.0
        tm[5, 1] = 0.0
        tm[5, 2] = 0.0
        tm[5, 3] = fr
        tm[5, 4] = 1.0
        tm[5, 5] = 1.0
    elif (anom == 'meana') or (anom == 'meann'):
        tm[5, 0] = 0.0
        tm[5, 1] = 0.0
        tm[5, 2] = 0.0
        tm[5, 3] = fr
        tm[5, 4] = 1.0
        tm[5, 5] = 1.0

    # ---------- calculate the output covariance matrix -----------
    eqcov = tm*classcov*tm.T
    return eqcov, tm

# ----------------------------------------------------------------------------
#
#                           function coveq2ct
#
#  this function transforms a six by six covariance matrix expressed in
#    equinoctial elements into one expressed in cartesian elements.
#
#  author        : david vallado                  719-573-2600   24 jul 2003
#
#  revisions
#    vallado     - major update                                  26 aug 2015
#
#  inputs          description                    range / units
#    eqcov       - 6x6 equinoctial covariance matrix
#    eqstate     - 6x1 equinoctial orbit state    (a/n af ag chi psi lm/ln)
#    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
#    fr          - retrograde factor               +1, -1
#
#  outputs       :
#    cartcov     - 6x6 cartesian covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    n           - mean motion                    rad
#    af          - component of ecc vector
#    ag          - component of ecc vector
#    chi         - component of node vector in eqw
#    psi         - component of node vector in eqw
#    meanlon     - mean longitude                 rad
#    nu          - true anomaly                   0.0  to 2pi rad
#    m           - mean anomaly                   0.0  to 2pi rad
#    r           - matrix of partial derivatives
#    e0          - eccentric anomaly              0.0  to 2pi rad
#
#  coupling      :
#
#  references    :
#    Vallado and Alfano 2015
#
# [cartcov, tm] = coveq2ct(eqcov, eqstate, anom, fr)
# ----------------------------------------------------------------------------

def coveq2ct(eqcov, eqstate, anom, fr):

    # --------- determine which set of variables is in use ---------
    # -------- parse the orbit state
    if (anom == 'truea') or (anom == 'meana'):
        a = eqstate[0] * 1000.0  # in m
        n = np.sqrt(mum/a**3)  # mum
    elif (anom == 'truen') or (anom == 'meann'):
        n = eqstate[0]  # rad
        a = (mum / n**2)**(1/3)   # in m
    af = eqstate[1]
    ag = eqstate[2]
    chi = eqstate[3]
    psi = eqstate[4]
    if (anom == 'meana') or (anom == 'meann'):
        meanlonM = eqstate[5]
        omega = math.atan2(chi, psi)
        argp = math.atan2(ag, af) - fr*math.atan2(chi, psi)
        ecc = np.sqrt (af**2 + ag**2)
        m = meanlonM - fr*omega - argp
        [eccanom, nu] = smu.newtonm (ecc, m)
        meanlonNu = nu + fr*omega + argp
    elif (anom == 'truea') or (anom == 'truen'):
        meanlonNu = eqstate[5]
        omega = math.atan2(chi, psi)
        argp = math.atan2(ag, af) - fr*math.atan2(chi, psi)
        nu = meanlonNu - fr*omega - argp
        nu = np.fmod(nu + twopi, twopi)
        ecc = np.sqrt (af**2 + ag**2)
        [eccanom, m] = smu.newtonnu (ecc, nu)
        meanlonM = fr*omega + argp + m
        meanlonM = np.fmod(meanlonM, 2.0*np.pi)

    # needs to be mean longitude for eq2rv
    [reci, veci] = eq2rv(a/1000.0, af, ag, chi, psi, meanlonM, fr)
    rx = reci[0] * 1000.0  # in m
    ry = reci[1] * 1000.0
    rz = reci[2] * 1000.0
    vx = veci[0] * 1000.0
    vy = veci[1] * 1000.0
    vz = veci[2] * 1000.0

    magr = smu.mag(reci)*1000 # m
    magv = smu.mag(veci)*1000 # m/s

    A = n * a**2
    B = np.sqrt(1.0 - ag**2 - af**2)
    C = 1.0 + chi**2 + psi**2
    b = 1.0 / (1.0 + B)

    # same
    #G = n * a**2 * np.sqrt(1.0 - af**2 - ag**2)
    G = A*B

    # -----------  initial guess -------------
    F0 = meanlonM
    numiter = 25
    ktr = 1
    F1 = F0 - ((F0 + ag*np.cos(F0) - af*np.sin(F0) - meanlonM)
               / (1.0 - ag*np.sin(F0) - af*np.cos(F0)))
    while ((abs(F1-F0) > small) and (ktr <= numiter)):
        ktr = ktr + 1
        F0= F1
        F1 = F0 - ((F0 + ag*np.cos(F0) - af*np.sin(F0) - meanlonM)
                   / (1.0 - ag*np.sin(F0) - af*np.cos(F0)))

    F = F1

    X = a*((1.0 - ag**2 * b) * np.cos(F) + af*ag*b*np.sin(F) - af)
    Y = a*((1.0 - af**2 * b) * np.sin(F) + af*ag*b*np.cos(F) - ag)

    XD = n*a**2/magr * (af*ag*b*np.cos(F) - (1.0 - ag**2*b)*np.sin(F))
    YD = n*a**2/magr * ((1.0 - af**2*b)*np.cos(F) - af*ag*b*np.sin(F))

    # alt formulation, but needs L
    #XD = -n*a*(ag + sinL) / (np.sqrt(1.0 - af**2 - ag**2))
    #YD = n*a*(af + cosL) / (np.sqrt(1.0 - af**2 - ag**2))

    # alt forms all are the same now
    #         sinL = ((1-af**2*b)*np.sin(F) + ag*af*b*np.cos(F) - ag) / (1 - ag*np.sin(F) - af*np.cos(F))
    #         cosL = ((1-ag**2*b)*np.cos(F) + ag*af*b*np.sin(F) - af) / (1 - ag*np.sin(F) - af*np.cos(F))
    #         XD = -n*a*(ag + sinL) / (B)
    #         YD = n*a*(af + cosL) / (B)
    #         r = a*(1-af**2-ag**2) / (1 + ag*sinL + af*cosL)
    #         r = a*(1.0 - ag*np.sin(F) - af*np.cos(F))
    #         r = magr
    #         X = r*cosL
    #         Y = r*sinL

    # components of equinoctial system
    p0 = 1.0 / (1.0 + chi**2 + psi**2)
    fe = p0 * (1.0 - chi**2 + psi**2)  # 2nd one is minus???? no, + seems correct
    fq = p0 * 2.0 * chi * psi
    fw = p0 * -2.0 * fr*chi
    fvec = [fe, fq, fw]
    ge = p0 * 2.0 * fr*chi * psi
    gq = p0 * fr*(1.0 + chi**2 - psi**2)
    gw = p0 * 2.0 * psi
    gvec = [ge, gq, gw]
    we = p0 * 2.0 * chi
    wq = p0 * -2.0 * psi
    ww = p0 * fr*(1.0 - chi**2 - psi**2)
    wvec = [we, wq, ww]

    partXaf = ag*b*XD/n + a*Y*XD/G - a
    partYaf = ag*b*YD/n - a*X*XD/G
    partXag = -af*b*XD/n + a*Y*YD/G
    partYag = -af*b*YD/n - a*X*YD/G - a

    partXDaf = a*XD*YD / G - A/(magr**3) * (a*ag*X/(1 + B) + X*Y/B)
    partYDaf = -a*XD**2 / G - A/(magr**3) * (a*ag*Y/(1 + B) - X**2/B)
    partXDag = a*YD**2 / G + A/(magr**3) * (a*af*X/(1 + B) - Y**2/B)
    partYDag = -a*XD*YD / G + A/(magr**3) * (a*af*Y/(1 + B) + X*Y/B)

    print('na %11.7f %11.7f  \n'%(n, a))
    print('XY %11.7f %11.7f %11.7f %11.7f\n'%(X, Y, XD, YD))
    print('ABC %11.7f %11.7f %11.7f \n'%(A, B, C))
    print('part %11.7f %11.7f %11.7f %11.7f \n'
          %(partXDaf, partYDaf, partXDag, partYDag))

    dMdnu = (1.0 - ecc**2)**1.5 / ((1.0 + ecc*np.cos(nu))**2)  # dm/dv

    # ---------------- calculate matrix elements ------------------
    # ---- partials of (rx ry rz vx vy vz) wrt a
    if (anom == 'truea') or (anom == 'meana'):
        p0 = 1.0 / a
        p1 = -1.0 / (2.0*a)
    elif (anom == 'truen') or (anom == 'meann'):
        p0 = -2.0/(3.0*n)
        p1 = 1.0/(3.0*n)

    tm = np.zeros((6, 6))
    tm[0, 0] = p0 * rx
    tm[1, 0] = p0 * ry
    tm[2, 0] = p0 * rz
    tm[3, 0] = p1 * vx
    tm[4, 0] = p1 * vy
    if (anom == 'meana') or (anom == 'meann'):
        tm[5, 0] = p1 * vz
    elif (anom == 'truea') or (anom == 'truen'):
        #                 cosl = np.cos(meanlonNu)
        #                 sinl = np.sin(meanlonNu)
        #                 ch = a*(2.0*af*ag*cosl + 2.0*ag + (ag**2 - af**2 + 1.0)*sinl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 ck = a*(2.0*af*ag*sinl + 2.0*af + (af**2 - ag**2 + 1.0)*cosl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 dhkv = 1.0*(chi*cosl - psi*sinl)
        #
        #                 tm(1:3, 2) = dXdh * fe    + dYdh * ge
        #                 sensh = (ch*dhkv[0] - tm[0, 1])/tm[0, 5]
        #                 tm(:, 2) = tm(:, 2) + tm(:, 6)*sensh
        #
        #                 tm(1:3, 3) = dXdk * fe    + dYdk * ge
        #                 sensk = (ck * dhkVec[0] - tm[0, 2])/tm[0, 5]
        #                 tm(:, 3) = tm(:, 3) + tm(:, 6)*sensk
        #
        tm[5, 0] = p1 * vz
        #tm[5, 0] = dMdnu * tm[5, 0]
        #tm[5, 0] = 0.0

    # ---- partials of (rx ry rz vx vy vz) wrt af
    tm[0, 1] = partXaf * fe + partYaf * ge
    tm[1, 1] = partXaf * fq + partYaf * gq
    tm[2, 1] = partXaf * fw + partYaf * gw
    tm[3, 1] = partXDaf * fe + partYDaf * ge
    tm[4, 1] = partXDaf * fq + partYDaf * gq
    if (anom == 'meana') or (anom == 'meann'):
        tm[5, 1] = partXDaf * fw + partYDaf * gw
    elif (anom == 'truea') or (anom == 'truen'):
        #                 cosl = np.cos(meanlonNu)
        #                 sinl = np.sin(meanlonNu)
        #                 ch = a*(2.0*af*ag*cosl + 2.0*ag + (ag**2 - af**2 + 1.0)*sinl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 ck = a*(2.0*af*ag*sinl + 2.0*af + (af**2 - ag**2 + 1.0)*cosl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 dhkv = 1.0*(chi*cosl - psi*sinl)
        #
        #                 tm(1:3, 2) = dXdh * fq    + dYdh * gq
        #                 sensh = (ch*dhkv[0] - tm[0, 1])/tm[0, 5]
        #                 tm(:, 2) = tm(:, 2) + tm(:, 6)*sensh
        #
        #                 tm(1:3, 3) = dXdk * fq    + dYdk * gq
        #                 sensk = (ck * dhkVec[0] - tm[0, 2])/tm[0, 5]
        #                 tm(:, 3) = tm(:, 3) + tm(:, 6)*sensk
        #
        tm[5, 1] = partXDaf * fw + partYDaf * gw
        #tm[5, 1] = dMdnu * tm[5, 1]
        #tm[5, 1] = 0.0

    # ---- partials of (rx ry rz vx vy vz) wrt ag
    tm[0, 2] = partXag * fe + partYag * ge
    tm[1, 2] = partXag * fq + partYag * gq
    tm[2, 2] = partXag * fw + partYag * gw
    tm[3, 2] = partXDag * fe + partYDag * ge
    tm[4, 2] = partXDag * fq + partYDag * gq
    if (anom == 'meana') or (anom == 'meann'):
        tm[5, 2] = partXDag * fw + partYDag * gw
    elif (anom == 'truea') or (anom == 'truen'):
        #                 cosl = np.cos(meanlonNu)
        #                 sinl = np.sin(meanlonNu)
        #                 ch = a*(2.0*af*ag*cosl + 2.0*ag + (ag**2 - af**2 + 1.0)*sinl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 ck = a*(2.0*af*ag*sinl + 2.0*af + (af**2 - ag**2 + 1.0)*cosl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 dhkv = 1.0*(chi*cosl - psi*sinl)
        #
        #                 tm(1:3, 2) = dXdh * fw    + dYdh * gw
        #                 sensh = (ch*dhkv[0] - tm[0, 1])/tm[0, 5]
        #                 tm(:, 2) = tm(:, 2) + tm(:, 6)*sensh
        #
        #                 tm(1:3, 3) = dXdk * fw    + dYdk * gw
        #                 sensk = (ck * dhkVec[0] - tm[0, 2])/tm[0, 5]
        #                 tm(:, 3) = tm(:, 3) + tm(:, 6)*sensk
        #
        tm[5, 2] = partXDag * fw + partYDag * gw
        #tm[5, 2] = dMdnu * tm[5, 2]
        #tm[5, 2] = 0.0

    # ---- partials of (rx ry rz vx vy vz) wrt chi
    p0 = 2.0 * fr / C
    tm[0, 3] = p0 * (psi * (Y*fe - X*ge) - X*we)
    tm[1, 3] = p0 * (psi * (Y*fq - X*gq) - X*wq)
    tm[2, 3] = p0 * (psi * (Y*fw - X*gw) - X*ww)
    tm[3, 3] = p0 * (psi * (YD*fe - XD*ge) - XD*we)
    tm[4, 3] = p0 * (psi * (YD*fq - XD*gq) - XD*wq)
    if (anom == 'meana') or (anom == 'meann'):
        tm[5, 3] = p0 * (psi * (YD*fw - XD*gw) - XD*ww)
    elif (anom == 'truea') or (anom == 'truen'):
        #                 cosl = np.cos(meanlonNu)
        #                 sinl = np.sin(meanlonNu)
        #                 ch = a*(2.0*af*ag*cosl + 2.0*ag + (ag**2 - af**2 + 1.0)*sinl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 ck = a*(2.0*af*ag*sinl + 2.0*af + (af**2 - ag**2 + 1.0)*cosl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 dhkv = 1.0*(chi*cosl - psi*sinl)
        #
        #                 tm(4:6, 2) = dXdotdh * fe + dYdotdh * ge
        #                 sensh = (ch*dhkv[0] - tm[0, 1])/tm[0, 5]
        #                 tm(:, 2) = tm(:, 2) + tm(:, 6)*sensh
        #
        #                 tm(4:6, 3) = dXdotdk * fe + dYdotdk * ge
        #                 sensk = (ck * dhkVec[0] - tm[0, 2])/tm[0, 5]
        #                 tm(:, 3) = tm(:, 3) + tm(:, 6)*sensk
        #
        tm[5, 3] = p0 * (psi * (YD*fw - XD*gw) - XD*ww)
        #tm[5, 3] = dMdnu * tm[5, 3]
        #tm[5, 3] = 0.0

    # ---- partials of (rx ry rz vx vy vz) wrt psi
    p0 = 2.0 * fr / C
    tm[0, 4] = p0 * (chi * (X*ge - Y*fe) + Y*we)
    tm[1, 4] = p0 * (chi * (X*gq - Y*fq) + Y*wq)
    tm[2, 4] = p0 * (chi * (X*gw - Y*fw) + Y*ww)
    tm[3, 4] = p0 * (chi * (XD*ge - YD*fe) + YD*we)
    tm[4, 4] = p0 * (chi * (XD*gq - YD*fq) + YD*wq)
    if (anom == 'meana') or (anom == 'meann'):
        tm[5, 4] = p0 * (chi * (XD*gw - YD*fw) + YD*ww)
    elif (anom == 'truea') or (anom == 'truen'):
        #                 cosl = np.cos(meanlonNu)
        #                 sinl = np.sin(meanlonNu)
        #                 ch = a*(2.0*af*ag*cosl + 2.0*ag + (ag**2 - af**2 + 1.0)*sinl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 ck = a*(2.0*af*ag*sinl + 2.0*af + (af**2 - ag**2 + 1.0)*cosl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 dhkv = 1.0*(chi*cosl - psi*sinl)
        #
        #                 tm(4:6, 2) = dXdotdh * fq + dYdotdh * gq
        #                 sensh = (ch*dhkv[0] - tm[0, 1])/tm[0, 5]
        #                 tm(:, 2) = tm(:, 2) + tm(:, 6)*sensh
        #
        #                 tm(4:6, 3) = dXdotdk * fq + dYdotdk * gq
        #                 sensk = (ck * dhkVec[0] - tm[0, 2])/tm[0, 5]
        #                 tm(:, 3) = tm(:, 3) + tm(:, 6)*sensk
        #
        tm[5, 4] = p0 * (chi * (XD*gw - YD*fw) + YD*ww)
        #tm[5, 4] = dMdnu * tm[5, 4]
        #tm[5, 4] = 0.0

    # ---- partials of (rx ry rz vx vy vz) wrt meanlon
    p0 = 1.0 / n
    p1 = - n * a**3 / magr**3
    tm[0, 5] = p0 * vx
    tm[1, 5] = p0 * vy
    tm[2, 5] = p0 * vz
    tm[3, 5] = p1 * rx
    tm[4, 5] = p1 * ry
    if (anom == 'meana') or (anom == 'meann'):
        tm[5, 5] = p1 * rz
    elif (anom == 'truea') or (anom == 'truen'):
        #                 cosl = np.cos(meanlonNu)
        #                 sinl = np.sin(meanlonNu)
        #                 ch = a*(2.0*af*ag*cosl + 2.0*ag + (ag**2 - af**2 + 1.0)*sinl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 ck = a*(2.0*af*ag*sinl + 2.0*af + (af**2 - ag**2 + 1.0)*cosl) / (C*(1.0 + af*cosl + ag*sinl)**2)
        #                 dhkv = 1.0*(chi*cosl - psi*sinl)
        #
        #                 tm(4:6, 2) = dXdotdh * fw + dYdotdh * gw
        #                 sensh = (ch*dhkv[0] - tm[0, 1])/tm[0, 5]
        #                 tm(:, 2) = tm(:, 2) + tm(:, 6)*sensh
        #
        #                 tm(4:6, 3) = dXdotdk * fw + dYdotdk * gw
        #                 sensk = (ck * dhkVec[0] - tm[0, 2])/tm[0, 5]
        #                 tm(:, 3) = tm(:, 3) + tm(:, 6)*sensk
        #
        tm[5, 5] = p1 * rz
        #tm[5, 5] = dMdnu * tm[5, 3]
        #tm[5, 5] = 0.0

            # similar to ct2cl true...
    #             r_dot_v = np.sqrt(rx*vx + ry*vy + rz*vz)
    #             ecc_term = magv*magv - mum/magr
    #             ecc_x = (ecc_term*rx - r_dot_v*vx)/mum
    #             ecc_y = (ecc_term*ry - r_dot_v*vy)/mum
    #             ecc_z = (ecc_term*rz - r_dot_v*vz)/mum
    #             r_dot_e = np.sqrt(rx*ecc_x + ry*ecc_y + rz*ecc_z)
    #             nu_scale = -sign(r_dot_v)/np.sqrt(1-np.cos(nu)*np.cos(nu))
    #             magr3 = magr**3
    #             temp = ry*(vx*vy - mum*rx*ry/magr3) - rx*ecc_term + rz*(vx*vz - mum*rx*rz/magr3)
    #             temp = temp - rx*(vy*vy + vz*vz - mum/magr + mum*rx*rx/magr3) + vx*r_dot_v
    #             temp = -temp/(mum*magr*ecc) - rx*r_dot_e/(magr3*ecc) - tm[1, 0]*r_dot_e/(magr*ecc*ecc)
    #             tm[5, 0] = temp*nu_scale
    #             temp = rx*(vx*vy - mum*rx*ry/magr3) - ry*ecc_term + rz*(vy*vz - mum*ry*rz/magr3)
    #             temp = temp-ry*(vx*vx + vz*vz - mum/magr + mum*ry*ry/magr3) + vy*r_dot_v
    #             temp = -temp/(mum*magr*ecc) - ry*r_dot_e/(magr3*ecc) - tm[1, 1]*r_dot_e/(magr*ecc*ecc)
    #             tm[5, 1] = temp*nu_scale
    #             temp = rx*(vx*vz - mum*rx*rz/magr3) - rz*ecc_term + ry*(vy*vz - mum*ry*rz/magr3)
    #             temp = temp - rz*(vx*vx + vy*vy - mum/magr + mum*rz*rz/magr3) + vz*r_dot_v
    #             temp = -temp/(mum*magr*ecc) - rz*r_dot_e/(magr3*ecc) - tm[1, 2]*r_dot_e/(magr*ecc*ecc)
    #             tm[5, 2] = temp*nu_scale
    #             temp = ry*(rx*vy - 2*ry*vx) + rx*(ry*vy + rz*vz) + rz*(rx*vz - 2*rz*vx)
    #             temp = -temp/(mum*magr*ecc) - tm[1, 3]*r_dot_e/(magr*ecc*ecc)
    #             tm[5, 3] = temp*nu_scale
    #             temp = rx*(ry*vx - 2*rx*vy) + ry*(rx*vx + rz*vz) + rz*(ry*vz - 2*rz*vy)
    #             temp = -temp/(mum*magr*ecc) - tm[1, 4]*r_dot_e/(magr*ecc*ecc)
    #             tm[5, 4] = temp*nu_scale
    #             temp = rz*(rx*vx + ry*vy) + rx*(rz*vx - 2*rx*vz) + ry*(rz*vy - 2*ry*vz)
    #             temp = -temp/(mum*magr*ecc) - tm[1, 5]*r_dot_e/(magr*ecc*ecc)
    #             tm[5, 5] = temp*nu_scale

    # ---------- calculate the output covariance matrix -----------
    cartcov = tm * eqcov * tm.T
    return cartcov, tm


# ----------------------------------------------------------------------------
#
#                           function coveq2cl
#
#  this function transforms a six by six covariance matrix expressed in
#    equinoctial elements into one expressed in classical orbital elements.
#
#  author        : david vallado                  719-573-2600   18 jun 2002
#
#  revisions
#    vallado     - major update                                  26 aug 2015
#
#  inputs          description                    range / units
#    eqcov       - 6x6 equinoctial covariance matrix
#    eqstate     - 6x1 equinoctial orbit state    (a/n af ag chi psi lm/ln)
#    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
#    fr          - retrograde factor               +1, -1
#
#  outputs       :
#    classcov    - 6x6 classical covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    n           - mean motion                    rad
#    af          - component of ecc vector
#    ag          - component of ecc vector
#    chi         - component of node vector in eqw
#    psi         - component of node vector in eqw
#    meanlon     - mean longitude                 rad
#    mum         - gravitational parameter        m**3/s**2 NOTE Meters!
#
#  coupling      :
#
#  references    :
#    Vallado and Alfano 2015
#
# [classcov, tm] = coveq2cl eqcov, eqstate, anom, fr)
# ----------------------------------------------------------------------------

def coveq2cl(eqcov, eqstate, anom, fr):

    # -------- parse the orbit state
    # --------- determine which set of variables is in use ---------
    if (anom == 'truea') or (anom == 'meana'):
        a = eqstate[0]  # in m
        n = np.sqrt(mum/a**3)
    elif (anom == 'truen') or (anom == 'meann'):
        n = eqstate[0]
        a = (mum / n**2)**(1/3)
    af = eqstate[1]
    ag = eqstate[2]
    chi = eqstate[3]
    psi = eqstate[4]
    if (anom == 'meana') or (anom == 'meann'):
        meanlonM = eqstate[5]
    elif (anom == 'truea') or (anom == 'truen'):
        meanlonNu = eqstate[5]


    tm = np.zeros((6, 6))
    # ---------------- calculate matrix elements ------------------
    # ---- partials of a wrt (a/n af ag chi psi l)
    if (anom == 'truea') or (anom == 'meana'):
        tm[0, 0] = 1.0
    elif (anom == 'truen') or (anom == 'meann'):
        tm[0, 0] = -2.0/(3*n) * (mum/n**2)**(1/3)  # if class = a, equin = n

    tm[0, 1] = 0.0
    tm[0, 2] = 0.0
    tm[0, 3] = 0.0
    tm[0, 4] = 0.0
    tm[0, 5] = 0.0

    # ---- partials of ecc wrt (a/n af ag chi psi l)
    p0 = 1.0 / np.sqrt(af**2 + ag**2)
    tm[1, 0] = 0.0
    tm[1, 1] = p0*af
    tm[1, 2] = p0*ag
    tm[1, 3] = 0.0
    tm[1, 4] = 0.0
    tm[1, 5] = 0.0

    # ---- partials of incl wrt (a/n af ag chi psi l)
    p1 = 2.0*fr / ((1.0 + chi**2 + psi**2) * np.sqrt(chi**2 + psi**2))
    tm[2, 0] = 0.0
    tm[2, 1] = 0.0
    tm[2, 2] = 0.0
    tm[2, 3] = p1*chi
    tm[2, 4] = p1*psi
    tm[2, 5] = 0

    # ---- partials of omega wrt (a/n af ag chi psi l)
    p2 = 1.0 / (chi**2 + psi**2)
    tm[3, 0] = 0
    tm[3, 1] = 0
    tm[3, 2] = 0
    tm[3, 3] = p2*psi
    tm[3, 4] = -p2*chi
    tm[3, 5] = 0

    # ---- partials of argp wrt (a/n af ag chi psi l)
    p3 = 1.0 / (af**2 + ag**2)
    tm[4, 0] = 0.0
    tm[4, 1] = -p3*ag
    tm[4, 2] = p3*af
    tm[4, 3] = -fr*p2*psi
    tm[4, 4] = fr*p2*chi
    tm[4, 5] = 0.0

    # ---- partials of anomaly wrt (a/n af ag chi psi l)
    p4 = 1.0 / (af**2 + ag**2)
    if (anom == 'truea') or (anom == 'truen'):
        tm[5, 0] = 0.0
        tm[5, 1] = p4*ag  # p3 on these. same
        tm[5, 2] = -p4*af
        tm[5, 3] = 0.0
        tm[5, 4] = 0.0
        tm[5, 5] = 1.0
    elif (anom == 'meana') or (anom == 'meann'):
        tm[5, 0] = 0.0
        tm[5, 1] = p3*ag
        tm[5, 2] = -p3*af
        tm[5, 3] = 0.0
        tm[5, 4] = 0.0
        tm[5, 5] = 1.0

    # ---------- calculate the output covariance matrix -----------
    classcov = tm@eqcov@tm.T
    return classcov, tm

# ----------------------------------------------------------------------------
#
#                           function covfl2ct
#
#  this function transforms a six by six covariance matrix expressed in
#    flight elements into one expressed in cartesian elements.
#
#  author        : david vallado                  719-573-2600   27 may 2003
#
#  revisions
#    vallado     - major update                                  26 aug 2015
#
#  inputs          description                    range / units
#    flcov       - 6x6 flight covariance matrix
#    flstate     - 6x1 flight orbit state         (r v latgc lon fpa az)
#    anom        - anomaly                        'latlon', 'radec'
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    terms       - number of terms for ast calculation 0, 2
#    ddpsi       - delta psi correction to gcrf          rad
#    ddeps       - delta eps correction to gcrf          rad
#
#  outputs       :
#    cartcov     - 6x6 cartesian covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    tm           - matrix of partial derivatives
#    magr        - eci position vector magnitude  km
#    magv        - eci velocity vector magnitude  km/sec
#    latgc       - geocentric latitude            rad
#    lon         - longitude                      rad
#    fpa         - sat flight path angle          rad
#    az          - sat flight path az             rad
#    fpav        - sat flight path anglefrom vert rad
#    xe, ye, ze    - ecef position vector componentskm
#
#  coupling      :
#    ecef2eci    - convert eci vectors to ecef
#
#  references    :
#    Vallado and Alfano 2015
#
# [cartcov, tm] = covfl2ct(flcov, flstate, anom, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ----------------------------------------------------------------------------


def covfl2ct(flcov, flstate, anom, ttt, jdut1, lod, xp, yp, terms, ddpsi,
             ddeps):

    small = 0.00000001

    # -------- parse the input vectors into components
    lon = flstate[0] # these will come in as either lon/lat or rtasc/decl depending on anom1
    latgc = flstate[1]
    fpa = flstate[2]
    az = flstate[3]
    magr = flstate[4]  # already converted to m in setcov
    magv = flstate[5]

    cfpa = np.cos(fpa)
    sfpa = np.sin(fpa)
    caz = np.cos(az)
    saz = np.sin(az)

    # --------- determine which set of variables is in use ---------
    # need to get eci vector
    # seems like it's easier to simply convert the lat-lon to rtasc
    # decl
    recef = np.zeros(3)
    vecef = np.zeros(3)
    reci = np.zeros(3)
    veci = np.zeros(3)
    if (anom, 'latlon'):
        craf = np.cos(lon)  # earth fixed needed for the lon lat partials only
        sraf = np.sin(lon)
        cdf = np.cos(latgc)
        sdf = np.sin(latgc)
        recef[0] = magr*0.001*np.cos(latgc)*np.cos(lon)  # in km
        recef[1] = magr*0.001*np.cos(latgc)*np.sin(lon)
        recef[2] = magr*0.001*np.sin(latgc)
        # -------- convert r to eci
        # this vel is wrong but not needed except for special case ahead
        vecef[0] = magv*0.001*(-np.cos(lon)*np.sin(latgc)*caz*cfpa
                               - np.sin(lon)*saz*cfpa
                               + np.cos(lon)*np.cos(latgc)*sfpa) # m/s
        vecef[1] = magv*0.001*(-np.sin(lon)*np.sin(latgc)*caz*cfpa
                               + np.cos(lon)*saz*cfpa
                               + np.sin(lon)*np.cos(latgc)*sfpa)
        vecef[2] = magv*0.001*(np.sin(lon)*sfpa + np.cos(latgc)*caz*cfpa)
        aecef = np.zeros(3)
        [reci, veci, a] = ecef2eci(recef.T, vecef.T, aecef, ttt,
                                 jdut1, lod, xp, yp, terms, ddpsi, ddeps)
        reci = reci*1000.0  # in m
        veci = veci*1000.0  # in m/s

        temp = np.sqrt(reci[0]*reci[0] + reci[1]*reci[1])
        if (temp < small):
            rtasc = math.atan2(veci[1] , veci[0])
            temp
        else:
            rtasc = math.atan2(reci[1] , reci[0])
        #decl = math.atan2(reci[2] , np.sqrt(reci[0]**2 + reci[1]**2))
        decl = math.asin(reci[2]/magr)
        cra = np.cos(rtasc)
        sra = np.sin(rtasc)
        cd = np.cos(decl)
        sd = np.sin(decl)
    elif (anom, 'radec'):
        rtasc = lon  # these come in as rtasc decl in this case
        decl = latgc
        reci[0] = magr*np.cos(decl)*np.cos(rtasc)
        reci[1] = magr*np.cos(decl)*np.sin(rtasc)
        reci[2] = magr*np.sin(decl)
        veci[0] = magv*(-np.cos(rtasc)*np.sin(decl)*caz*cfpa
                        - np.sin(rtasc)*saz*cfpa
                        + np.cos(rtasc)*np.cos(decl)*sfpa) # m/s
        veci[1] = magv*(-np.sin(rtasc)*np.sin(decl)*caz*cfpa
                        + np.cos(rtasc)*saz*cfpa
                        + np.sin(rtasc)*np.cos(decl)*sfpa)
        veci[2] = magv*(np.sin(decl)*sfpa + np.cos(decl)*caz*cfpa)
        cra = np.cos(rtasc)
        sra = np.sin(rtasc)
        cd = np.cos(decl)
        sd = np.sin(decl)

    tm = np.zeros((6, 6))
    # ---------------- calculate matrix elements ------------------
    # ---- partials of rx wrt (lon latgc fpa az r v)
    if (anom, 'radec'):
        tm[0, 0] = -magr*cd*sra
        tm[0, 1] = -magr*sd*cra
    else:  # latlon
        tm[0, 0] = -magr*cdf*sraf
        tm[0, 1] = -magr*sdf*craf
    tm[0, 2] = 0.0
    tm[0, 3] = 0.0
    tm[0, 4] = cd*cra
    tm[0, 5] = 0.0

    # ---- partials of ry wrt (lon latgc fpa az r v)
    if (anom, 'radec'):
        tm[1, 0] = magr*cd*cra
        tm[1, 1] = -magr*sd*sra
    else:   # latlon
        tm[1, 0] = magr*cdf*craf
        tm[1, 1] = -magr*sdf*sraf
    tm[1, 2] = 0.0
    tm[1, 3] = 0.0
    tm[1, 4] = cd*sra
    tm[1, 5] = 0.0

    # ---- partials of rz wrt (lon latgc fpa az r v)
    if (anom, 'radec'):
        tm[2, 0] = 0.0
        tm[2, 1] = magr*cd
    else:   # latlon
         tm[2, 0] = 0.0
         tm[2, 1] = magr*cdf
    tm[2, 2] = 0.0
    tm[2, 3] = 0.0
    tm[2, 4] = sd
    tm[2, 5] = 0.0

    # ---- partials of vx wrt (lon latgc fpa az r v)
    if (anom, 'radec'):
        tm[3, 0] = -magv*(-sra*caz*sd*cfpa + cra*saz*cfpa + cd*sra*sfpa)
        #  tm[3, 0] = -vy
        tm[3, 1] = -cra*magv*(sd*sfpa + cd*caz*cfpa)
        #  tm[3, 1] = -vz*cra
    else:   # latlon
        tm[3, 0] = -magv*(-sraf*caz*sdf*cfpa + craf*saz*cfpa + cdf*sraf*sfpa)
        #  tm[3, 0] = -vy
        tm[3, 1] = -craf*magv*(sdf*sfpa + cdf*caz*cfpa)
        #  tm[3, 1] = -vz*cra
    tm[3, 2] = magv*(cra*caz*sd*sfpa + sra*saz*sfpa + cd*cra*cfpa)
    tm[3, 3] = magv*(cra*saz*sd*cfpa - sra*caz*cfpa)
    tm[3, 4] = 0.0
    tm[3, 5] = -cra*caz*sd*cfpa - sra*saz*cfpa + cd*cra*sfpa

    # ---- partials of vy wrt (lon latgc fpa az r v)
    if (anom, 'radec'):
        tm[4, 0] = magv*(-cra*caz*sd*cfpa - sra*saz*cfpa + cd*cra*sfpa)
        #  tm[4, 0] = vx
        tm[4, 1] = -sra*magv*(sd*sfpa + cd*caz*cfpa)
        #   tm[4, 1] = -vz*sra
    else:   # latlon
        tm[4, 0] = magv*(-craf*caz*sdf*cfpa - sraf*saz*cfpa + cdf*craf*sfpa)
        #  tm[4, 0] = vx
        tm[4, 1] = -sraf*magv*(sdf*sfpa + cdf*caz*cfpa)
        #   tm[4, 1] = -vz*sra
    tm[4, 2] = magv*(sra*caz*sd*sfpa - cra*saz*sfpa + cd*sra*cfpa)
    tm[4, 3] = magv*(sra*saz*sd*cfpa + cra*caz*cfpa)
    tm[4, 4] = 0.0
    tm[4, 5] = -sra*caz*sd*cfpa + cra*saz*cfpa + cd*sra*sfpa

    # ---- partials of vz wrt (lon latgc fpa az r v)
    if (anom, 'radec'):
        tm[5, 0] = 0.0
        tm[5, 1] = magv*(cd*sfpa - sd*caz*cfpa)
    else:   # latlon
        tm[5, 0] = 0.0
        tm[5, 1] = magv*(cdf*sfpa - sdf*caz*cfpa)
    tm[5, 2] = magv*(sd*cfpa - cd*caz*sfpa)
    tm[5, 3] = -magv*cd*saz*cfpa
    tm[5, 4] = 0.0
    tm[5, 5] = sd*sfpa + cd*caz*cfpa

    # ---------- calculate the output covariance matrix -----------
    cartcov = tm*flcov*tm.T
    return cartcov, tm

# ------------------------------------------------------------------------------
#
#                           function rv2tradc
#
#  this function converts geocentric equatorial (eci) position and velocity
#    vectors into range, topcentric right acension, declination, and rates.
#    notice the value of small as it can affect the rate term calculations.
#    the solution uses the velocity vector to find the singular cases. also,
#    the right acension and declination rate terms are not observable unless
#    the acceleration vector is available.
#
#  author        : david vallado                  719-573-2600   19 jul 2004
#
#  revisions
#
#  inputs          description                    range / units
#    reci        - eci position vector            km
#    veci        - eci velocity vector            km/s
#    latgd       - geodetic latitude              -pi/2 to pi/2 rad
#    lon         - longitude of site              -2pi to 2pi rad
#    alt         - altitude                       km
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    ddpsi       - delta psi correction to gcrf          rad
#    ddeps       - delta eps correction to gcrf          rad
#
#  outputs       :
#    rho         - satellite range from site      km
#    rtasc       - topocentric right ascension    0.0 to 2pi rad
#    tdecl       - topocentric declination        -pi/2 to pi/2 rad
#    drho        - range rate                     km/s
#    daz         - xxazimuth rate                 rad / s
#    del         - xxelevation rate               rad / s
#
#  locals        :
#    rhoveci     - eci range vector from site     km
#    drhoveci    - eci velocity vector from site  km / s
#    rhoeci      - eci range vector from site     km
#    drhoeci     - sez velocity vector from site  km
#    wcrossr     - cross product result           km / s
#    earthrate   - eci earth's rotation rate vec  rad / s
#    tempvec     - temporary vector
#    temp        - temporary real*8 value
#    temp1       - temporary real*8 value
#    i           - index
#
#  coupling      :
#    mag         - magnitude of a vector
#    rot3        - rotation about the 3rd axis
#    rot2        - rotation about the 2nd axis
#
#  references    :
#    vallado       2001, 250-255, alg 27
#
# [rho, trtasc, tdecl, drho, dtrtasc, dtdecl] = rv2tradc (reci, veci, latgd, lon, alt, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ------------------------------------------------------------------------------

def rv2tradc(reci: np.ndarray, veci: np.ndarray, latgd: float, lon: float,
             alt: float, ttt: float, jdut1: float, lod: float, xp: float,
             yp: float, terms: int, ddpsi: float, ddeps: float) :
    """this function converts geocentric equatorial (eci) position and velocity
    vectors into range, topcentric right acension, declination, and rates.
    notice the value of small as it can affect the rate term calculations.
    the solution uses the velocity vector to find the singular cases. also,
    the right acension and declination rate terms are not observable unless
    the acceleration vector is available.

    Parameters
    ----------
    reci : ndarray
        eci position vector: km
    veci : ndarray
        eci velocity vector: km/s
    latgd : float
        geodetic latitude of site: -pi/2 to pi/2 rads
    lon : float
        longitude of site: -2pi to 2pi rads
    alt : float
        altitude of site: km
    ttt : float
        julian centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days from 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    terms : int
        # of terms for ast calculation: 0 or 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta eps correction to gcrf: rad

    Returns
    -------

    rho: float
        satellite range from site: km
    trtasc: float
        topocentric right ascension: 0 to 2pi rads
    tdecl: float
        topocentric declination: -pi/2 to pi/2 rads
    drho: float:
        range rate: km/s
    dtrtasc: float

    dtdecl
    """
    # --------------------- implementation ------------------------
# ----------------- get site vector in ecef -------------------
    rsecef, vsecef = obu.site(latgd, lon, alt)

# -------------------- convert ecef to eci --------------------
    a = np.zeros(3)
    rseci, vseci, _ = ecef2eci(rsecef, vsecef, a, ttt,
                                jdut1, lod, xp, yp, 2, ddpsi, ddeps)


    # ------- find eci range vector from site to satellite -------
    rhoeci = reci - rseci
    drhoeci = veci - vseci
    rho = smu.mag(rhoeci)
    # ------------- calculate azimuth and elevation ---------------
    temp = np.sqrt(rhoeci[0]**2 + rhoeci[1]**2)
    if (temp < small):
        trtasc = math.atan2(drhoeci[1], drhoeci[0])
    else:
        trtasc = math.atan2(rhoeci[1], rhoeci[0])

    if ((temp < small)):
        tdecl = np.sign(rhoeci[2]) * halfpi
    else:
        magrhoeci = smu.mag(rhoeci)
        tdecl = math.asin(rhoeci[2] / magrhoeci)

    if (trtasc < 0.0):
        trtasc = trtasc + 2.0 * np.pi

    # ------ calculate range, azimuth and elevation rates ---------
    temp1 = -rhoeci[1]**2 - rhoeci[0]**2
    drho = np.dot(rhoeci, drhoeci) / rho
    if (np.abs(temp1) > small):
        dtrtasc = (drhoeci[0] * rhoeci[1] - drhoeci[1] * rhoeci[0]) / temp1
    else:
        dtrtasc = 0.0

    if (abs(temp) > small):
        dtdecl = (drhoeci[2] - drho * math.sin(tdecl)) / temp
    else:
        dtdecl = 0.0

    return rho, trtasc, tdecl, drho, dtrtasc, dtdecl



# ------------------------------------------------------------------------------
#
#                           function rv2rsw
#
#  this function converts position and velocity vectors into radial, tangential (in-
#  track), and normal (cross-track) coordinates. note that there are numerous
#  nomenclatures for these systems. this is the rsw system of vallado. the reverse
#  values are found using the transmat transpose.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#                -
#  inputs          description                    range / units
#    reci        - eci position vector            km
#    veci        - eci velocity vector            km/s
#
#  outputs       :
#    rrsw        - rsw position vector            km
#    vrsw        - rsw velocity vector            km/s
#    transmat    - transformation matrix used to rotate vectors
#
#  locals        :
#    temp        - temporary position vector
#
#  coupling      :
#    mag         - magnitude of a vector
#    matvecmult  - multiply a square matrix by a vector (single column)
#
#  references    :
#    vallado       2007, 172
#
# [rrsw, vrsw, transmat] = rv2rsw(reci, veci)
# ------------------------------------------------------------------------------


def rv2rsw(reci: np.ndarray, veci: np.ndarray):
    """this function converts position and velocity vectors into radial, tangential (in-
    track), and normal (cross-track) coordinates. note that there are numerous
    nomenclatures for these systems. this is the rsw system of vallado. the reverse
    values are found using the transmat transpose.

    Parameters
    ----------
    reci : ndarray
        eci position vector: km
    veci : ndarray
        eci velocity vector ECI: km/s

    Returns
    -------
    rrsw: ndarray
        rsw position vector: km
    vrsw: ndarray
        rsw velocity vector: km/s
    transmat: ndarray
        transformation matrix used to rotate vectors
    """
    # each of the components must be unit vectors
    # radial component
    rvec = smu.unit(reci)

    # cross-track component
    wvec = np.cross(reci, veci)
    wvec = smu.unit(wvec)

    # along-track component
    svec = np.cross(wvec, rvec)
    svec = smu.unit(svec)

    # assemble transformation matrix from eci to rsw frame (individual
    # components arranged in row vectors)
    transmat = np.zeros((3, 3))
    transmat[0, 0] = rvec[0]
    transmat[0, 1] = rvec[1]
    transmat[0, 2] = rvec[2]
    transmat[1, 0] = svec[0]
    transmat[1, 1] = svec[1]
    transmat[1, 2] = svec[2]
    transmat[2, 0] = wvec[0]
    transmat[2, 1] = wvec[1]
    transmat[2, 2] = wvec[2]

    rrsw = transmat @ reci
    vrsw = transmat @ veci
    # rrsw = smu.matvecmult(transmat, reci, 3)
    # vrsw = smu.matvecmult(transmat, veci, 3)
    return rrsw, vrsw, transmat



# ------------------------------------------------------------------------------
#
#                           function rv2ntw
#
#  this function converts position and velocity vectors into normal (in-radial),
#    tangential (velocity), and normal (cross-track) coordinates. note that sometimes
#    the first vector is called along-radial. the tangential direction is
#    always aligned with the velocity vector. this is the ntw system of
#    vallado.
#
#  author        : david vallado                  719-573-2600    5 jul 2002
#
#  revisions
#                -
#  inputs          description                    range / units
#    r           - position vector                km
#    v           - velocity vector                km/s
#
#  outputs       :
#    rntw        - ntw position vector            km
#    vntw        - ntw velocity vector            km/s
#
#  locals        :
#    temp        - temporary position vector
#
#  coupling      :
#    mag         - magnitude of a vector
#    matvecmult  - multiply a square matrix by a vector (single column)
#
#  references    :
#    vallado       2007, 172
#
# [rntw, vntw, transmat] = rv2ntw(r, v)
# ------------------------------------------------------------------------------


def rv2ntw(r, v):
    # compute satellite velocity vector magnitude
    vmag = smu.mag(v)
    # in order to work correctly each of the components must be
#  unit vectors !
# in-velocity component
    tvec = np.divide(v, vmag)
    # cross-track component
    wvec = np.cross(r, v)
    wvec = smu.unit(wvec)
    # along-radial component
    nvec = np.cross(tvec, wvec)
    nvec = smu.unit(nvec)
    # assemble transformation matrix from to ntw frame (individual
#  components arranged in row vectors)
    transmat = np.zeros((3, 3))
    transmat[0, 0] = nvec[0]
    transmat[0, 1] = nvec[1]
    transmat[0, 2] = nvec[2]
    transmat[1, 0] = tvec[0]
    transmat[1, 1] = tvec[1]
    transmat[1, 2] = tvec[2]
    transmat[2, 0] = wvec[0]
    transmat[2, 1] = wvec[1]
    transmat[2, 2] = wvec[2]
    rntw = smu.matvecmult(transmat, r, 3)
    vntw = smu.matvecmult(transmat, v, 3)
    return rntw, vntw, transmat



#
# ----------------------------------------------------------------------------
#
#                           function rv2flt
#
#  this function transforms a position and velocity vector into the flight
#    elements - latgc, lon, fpa, az, position and velocity magnitude.
#
#  author        : david vallado                  719-573-2600   17 jun 2002
#
#  revisions
#    vallado     - add terms for ast calculation                 30 sep 2002
#    vallado     - chg magr var names                            23 may 2003
#
#  inputs          description                    range / units
#    reci        - eci position vector            km
#    veci        - eci velocity vector            km/s
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    terms       - number of terms for ast calculation 0, 2
#    ddpsi, ddeps - corrections for fk5 to gcrf    rad
#
#  outputs       :
#    lon         - longitude                      rad
#    latgc       - geocentric latitude            rad
#    rtasc       - right ascension                rad
#    decl        - declination                    rad
#    fpa         - sat flight path angle          rad
#    az          - sat flight path az             rad
#    magr        - eci position vector magnitude  km
#    magv        - eci velocity vector magnitude  km/sec
#
#  locals        :
#    fpav        - sat flight path anglefrom vert rad
#
#  coupling      :
#    none        -
#
#  references    :
#    vallado       2001, xx
#
# [lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt(r, v, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ----------------------------------------------------------------------------

def rv2flt(reci: np.ndarray, veci: np.ndarray, ttt: float, jdut1: float,
           lod: float, xp: float, yp: float, terms: int, ddpsi: float,
           ddeps: float):
    """this function transforms a position and velocity vector into the flight
    elements - latgc, lon, fpa, az, position and velocity magnitude.

    Parameters
    ----------
    reci : ndarray
        eci position vector: km
    veci : np.ndarray
        eci velocity vector: km/s
    ttt : float
        julian centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days from 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    terms : int
        number of terms for ast calculation: 0, 2
    ddpsi : float
        corrections for fk5 to gcrf: rad
    ddeps : float
        corrections for fk5 to gcrf: rad

    Returns
    -------
    lon: float
        longitude: rad
    latgc: float
        geocentric latitude: rad
    rtasc: float
        right ascension: rad
    decl: float
        declination: rad
    fpa: float
        satellite flight path angle: rad
    az: float
        satellite flight path azimuth: rad
    magr: float
        eci position vector magnitude:  km
    magv: float
        eci velocity vector magnitude: km/sec
    """
    small = 1e-08
    magr = smu.mag(reci)
    magv = smu.mag(veci)
    # -------- convert r to ecef for lat/lon calculation
    avec = np.zeros(3)
    recef, vecef, _ = eci2ecef(reci, veci, avec, ttt,
                                 jdut1, lod, xp, yp, terms, ddpsi, ddeps)
    # ----------------- find longitude value  ----------------- uses ecef
    temp = math.sqrt(recef[0] ** 2 + recef[1] ** 2)
    if (temp < small):
        lon = math.atan2(vecef[1], vecef[0])
    else:
        lon = math.atan2(recef[1], recef[0])

    #latgc = math.atan2(recef[2] , sqrt(recef[0]**2 + recef[1]**2))
    latgc = math.asin(recef[2] / magr)
    # ------------- calculate rtasc and decl ------------------ uses eci
    temp = math.sqrt(reci[0] ** 2 + reci[1] ** 2)
    if (temp < small):
        rtasc = math.atan2(veci[1], veci[0])
    else:
        rtasc = math.atan2(reci[1], reci[0])

    #decl = math.atan2(reci[2] , sqrt(reci[0]**2 + reci[1]**2))
    decl = math.asin(reci[2] / magr)
    h = np.cross(reci, veci)
    hmag = smu.mag(h)
    rdotv = np.dot(reci, veci)
    fpav = math.atan2(hmag, rdotv)
    fpa = np.pi * 0.5 - fpav
    hcrossr = np.cross(h, reci)
    az = math.atan2(reci[0] * hcrossr[1] - reci[1] * hcrossr[0],
                    hcrossr[2] * magr)
    return lon, latgc, rtasc, decl, fpa, az, magr, magv



#
# ----------------------------------------------------------------------------
#
#                           function rv2adbar
#
#  this function transforms a position and velocity vector into the adbarv
#    elements - rtasc, decl, fpav, azimuth, position and velocity magnitude.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    reci        - eci position vector            km
#    veci        - eci velocity vector            km/s
#
#  outputs       :
#    magr        - eci position vector magnitude  km
#    magv        - eci velocity vector magnitude  km/sec
#    rtasc       - right ascension of sateillite  rad
#    decl        - declination of satellite       rad
#    fpav        - sat flight path angle from vertrad
#    az          - sat flight path azimuth        rad
#
#  locals        :
#    none        -
#
#  coupling      :
#    none        -
#
#  references    :
#    vallado       2001, xx
#    chobotov            70
#
# [rmag, vmag, rtasc, decl, fpav, az] = rv2adbar (r, v)
# ----------------------------------------------------------------------------

def rv2adbar(reci: np.ndarray, veci: np.ndarray):
    """this function transforms a position and velocity vector into the adbarv
    elements - rtasc, decl, fpav, azimuth, position and velocity magnitude.

    Parameters
    ----------
    reci : ndarray
        eci position vector: km
    veci : ndarray
        eci velocity vector: km/s

    Returns
    -------
    magr: float
        eci position vector magnitude: km
    magv: float
        eci velocity vector magnitude: km/sec
    rtasc: float
        right ascension of sateillite: rad
    decl: float
        declination of satellite: rad
    fpav: float
        sat flight path angle from vert: rad
    az: float
        sat flight path azimuth: rad
    """

    small = 1e-08
    magr = smu.mag(reci)
    magv = smu.mag(veci)
    # ---------------- calculate rtasc and decl -------------------
    temp = math.sqrt(reci[0] * reci[0] + reci[1] * reci[1])
    if (temp < small):
        temp1 = math.sqrt(veci[0] * veci[0] + veci[1] * veci[1])
        if (np.abs(temp1) > small):
            rtasc = math.atan2(veci[1], veci[0])
        else:
            rtasc = 0.0
    else:
        rtasc = math.atan2(reci[1], reci[0])

    decl = math.asin(reci[2] / magr)
    h = np.cross(reci, veci)
    hmag = smu.mag(h)
    rdotv = np.dot(reci, veci)
    fpav = math.atan2(hmag, rdotv)
    hcrossr = np.cross(h, reci)
    az = math.atan2(reci[0] * hcrossr[1] - reci[1] * hcrossr[0], hcrossr[2] * magr)
    return magr, magv, rtasc, decl, fpav, az


# ----------------------------------------------------------------------------
#
#                           function printcov
#
#  this function prints a covariance matrix
#
#  author        : david vallado                  719-573-2600   23 may 2003
#
#  revisions
#
#  inputs          description                    range / units
#    covin       - 6x6 input covariance matrix
#    covtype     - type of covariance             'cl', 'ct', 'fl', 'sp', 'eq',
#    cu          - covariance units (deg or rad)  't' or 'm'
#    anom        - anomaly                        'mean' or 'true' or 'tau'
#
#  outputs       :
#
#  locals        :
#
#  references    :
#    none
#
#  printcov(covin, covtype, cu, anom)
# ----------------------------------------------------------------------------

def printcov(covin: np.ndarray, covtype: str, cu: str, anom: str):
    """this function prints a covariance matrix.

    Parameters
    ----------
    covin : ndarray
        6x6 input covariance matrix
    covtype : str
        type of covariance: 'cl', 'ct', 'fl', 'sp', 'eq'
    cu : str
        covariance units: 't' or 'm'
    anom : str
        anomaly: 'mean' or 'true' or 'tau'
    """
    print("in printcov, anom is ", anom)
    if (anom == 'truea') or (anom == 'meana'):
        semi = 'a m  '
    elif (anom == 'truen') or (anom == 'meann'):
        semi = 'n rad'

    if covtype == 'ct':
        print('cartesian covariance \n')
        print('        x  m            y m             z  m           '
              'xdot  m/s       ydot  m/s       zdot  m/s  \n')

    if covtype == 'cl':
        print('classical covariance \n')
        if (cu == 'm'):
            print('          %s          ecc           incl rad      '
                  'raan rad         argp rad        ' % (semi))
            if (anom == 'meana') or (anom == 'meann'):
                print(' m rad \n')
            elif (anom == 'truea') or (anom == 'truen'):
                print(' nu rad \n')
        else:
            print('          %s           ecc           incl deg      '
                  'raan deg         argp deg        ' % (semi))
            if (anom == 'meana') or (anom == 'meann'):
                print(' m deg \n')
            elif (anom == 'truea') or (anom == 'truen'):
                print(' nu deg \n')

    if covtype == 'eq':
        print('equinoctial covariance \n')
        #            if (cu == 'm')
        if (anom == 'meana') or (anom == 'meann'):
            print('         %5s           af              ag           '
                  'chi             psi         meanlonM rad\n' % (semi))
        elif (anom == 'truea') or (anom == 'truen'):
            print('         %5s           af              ag           '
                  'chi             psi         meanlonNu rad\n' % (semi))

    if covtype == 'fl':
        print('flight covariance \n')
        if (cu == 'm'):
            print('       lon  rad      latgc rad        fpa rad         '
                  'az rad           r  m           v  m/s  \n')
        else:
            print('       lon  deg      latgc deg        fpa deg         '
                  'az deg           r  m           v  m/s  \n')

    if covtype == 'sp':
        print('spherical covariance \n')
        if (cu == 'm'):
            print('      rtasc rad       decl rad        fpa rad         '
                  'az rad           r  m           v  m/s  \n')
        else:
            print('      rtasc deg       decl deg        fpa deg         '
                  'az deg           r  m           v  m/s  \n')

    print("covin:")
    print(covin)
    #        print('#16e#16e#16e#16e#16e#16e\n', covin')
    print('covin transpose\n', (covin.T))

# ----------------------------------------------------------------------------
#
#                           function setcov
#
#  this function sets the intilal states for the covariance calculations.
#
#  author        : david vallado                  719-573-2600   21 jul 2003
#
#  revisions
#
#  inputs          description                    range / units
#    reci        - position vector                km
#    veci        - velocity vector                km/s
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    terms       - number of terms for ast calculation 0, 2
#    anom, anom1 - anomaly values                 meann, meana, turen,
#    ddpsi       - delta psi correction to gcrf          rad
#    ddeps       - delta eps correction to gcrf          rad

#    truea; latlon, radec
#
#  outputs       :
#    cartstate   - 6x1 cartesian state            m, m/s
#    classstate  - 6x1 classical orbital elements
#    flstate     - 6x1 flight orbital elements
#    eqstate     - 6x1 equinoctial orbital elements
#    fr          - retrograde factor for orbits with inlc > 90 deg  1, -1
#
#  locals        :
#
#  coupling      :
#    none
#
#  references    :
#    none
#
#  [cartstate, classstate, flstate, eqstate] = setcov(reci, veci, ...
#                                           year, mon, day, hr, min, sec, dut1, dat, ...
#                                           ttt, jdut1, lod, xp, yp, terms, printopt, anom)
# ----------------------------------------------------------------------------

def setcov(reci: np.ndarray, veci: np.ndarray, ttt: float,
           jdut1: float, lod: float, xp: float, yp: float, terms: int,
           printopt: str, anom: str, anom1: str, ddpsi: float, ddeps: float):

    fr = 1.0
    cartstate = np.array([[reci[0]], [reci[1]], [reci[2]],
                          [veci[0]], [veci[1]], [veci[2]]])
    print("in setcov, reci aand veci are:")
    print(reci)
    print(veci)

    # -------- convert to a classical orbit state
    p, a, ecc, incl, omega, argp, nu, m, _, _, _ = rv2coe(reci, veci)
    classstate = np.zeros(6)
    classstate[0] = a * 1000
    classstate[1] = ecc
    classstate[2] = incl
    classstate[3] = omega
    classstate[4] = argp
    if (anom == 'meana') or (anom == 'meann'):
        classstate[5] = m
    elif (anom == 'truea') or (anom == 'truen'):
        classstate[5] = nu

    # -------- convert to a flight orbit state
    lon, latgc, rtasc, decl, fpa, az, magr, magv = \
        rv2flt(reci, veci, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
    flstate = np.zeros(6)
    if anom1 == 'radec':
        flstate[0] = rtasc
        flstate[1] = decl
    elif anom1 == 'latlon':
        flstate[0] = lon
        flstate[1] = latgc

    flstate[2] = fpa
    flstate[3] = az
    flstate[4] = magr * 1000

    flstate[5] = magv * 1000
    # test position and velocity going back
    # avec = np.array([[0.0], [0.0], [0.0]])
    # recef, vecef, aecef = eci2ecef(reci, veci, avec, ttt, jdut1,
    #                              lod, xp, yp, terms, ddpsi, ddeps)

    # sinlon, coslon, sinlat, coslat,
    sinaz, cosaz, sinfpa, cosfpa, sinrta, cosrta, sindec, cosdec = \
        smu.getsincos(az, fpa, rtasc, decl)

    # vx = magv * (-coslon * sinlat * cosaz * cosfpa - sinlon * sinaz * cosfpa
    #              + coslon * coslat * sinfpa)
    # vy = magv * (-sinlon * sinlat * cosaz * cosfpa + coslon * sinaz * cosfpa
    #              + sinlon * coslat * sinfpa)
    # vz = magv * (sinlat * sinfpa + coslat * cosaz * cosfpa)

    # correct:
    ve1 = magv * (-cosrta * sindec * cosaz * cosfpa - sinrta * sinaz * cosfpa
                  + cosrta * cosdec * sinfpa)
    ve2 = magv * (-sinrta * sindec * cosaz * cosfpa + cosrta * sinaz * cosfpa
                  + sinrta * cosdec * sinfpa)
    ve3 = magv * (sindec * sinfpa + cosdec * cosaz * cosfpa)

    # -------- convert to an equinoctial orbit state
    a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr = rv2eq(reci, veci)
    eqstate = np.zeros(6)
    if (anom == 'meana') or (anom == 'truea'):
        eqstate[0] = a
    elif (anom == 'meann') or (anom == 'truen'):
        eqstate[0] = n

    eqstate[1] = af
    eqstate[2] = ag
    eqstate[3] = chi
    eqstate[4] = psi
    if (anom == 'meana') or (anom == 'meann'):
        eqstate[5] = meanlonM
    elif (anom == 'truea') or (anom == 'truen'):
        eqstate[5] = meanlonNu

    if printopt == 'y':
        # --------------------- write out input data --------------------------
        # test velocity prints
        # vtecef = np.array([vx, vy, vz])
        vteci = np.array([ve1, ve2, ve3])

        if smu.mag(vteci - np.transpose(veci)) > 0.01:
            print('ERROR in test of vel in setcov %11.7f \n'
                  % (smu.mag(vteci - veci.T)))
        print('input data \n' % ())
        print(' re %8.6f km  \n' % (re))
        print(' mu %14.8f km3/s2  \n' % (mu))
        # print('year %5i ' % (year))
        # print('mon %4i ' % (mon))
        # print('day %3i ' % (day))
        # print('hr %3i:%2i:%8.6f\n' % (hr, min, sec))
        # print('dut1 %8.6f s' % (dut1))
        # print(' dat %3i s' % (dat))
        print(' xp %8.6f "' % (xp))
        print(' yp %8.6f "' % (yp))
        print(' lod %8.6f s\n' % (lod))
        print('r :')
        print(reci)
        print('v :')
        print(veci)
        print('          p km       a km      ecc      incl deg    ')
        print(' raan deg     argp deg      nu deg      m deg \n')
        print('coes %11.4f %11.4f %11.7f %11.5f %11.5f'
              % (p, a, ecc, incl * rad2deg, omega * rad2deg))
        print('%11.5f %11.5f %11.5f\n'
              % (argp * rad2deg, nu * rad2deg, m * rad2deg))
        print('          a (km)       af           ag')
        print('           chi        psi            meanlonM     '
              'meanLonNu   fr\n')
        print('eq   %14.7f %14.7f %14.7f %15.7f %14.7f %14.7f %14.7f %2.0f\n'
              % (a, af, ag, chi, psi, meanlonM * rad2deg, meanlonNu * rad2deg,
                 fr))
        print('       lon deg       latgc deg     rtasc deg      '
              'decl deg      fpa deg       ')
        print(' az deg       magr km      magv km/s\n')
        print('flt  %14.7f%14.7f%14.7f%14.7f%14.7f%15.7f%14.7f%14.7f\n'
              % (lon * rad2deg, latgc * rad2deg, rtasc * rad2deg, decl
                 * rad2deg, fpa * rad2deg, az * rad2deg, magr, magv))

    return cartstate, classstate, flstate, eqstate, fr

# ----------------------------------------------------------------------------
#
#                           function covcl2ct
#
#  this function transforms a six by six covariance matrix expressed in classical elements
#    into one expressed in cartesian elements
#
#  author        : david vallado
#
#  revisions
#    vallado     - simplify code using pqw-eci transformation    12 may 2017
#    vallado     - major update                                  26 aug 2015
#
#  inputs          description                    range / units
#    classcov    - 6x6 classical covariance matrix
#    classstate  - 6x1 classical orbit state      (a e i O w nu/m)
#    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
#
#  outputs       :
#    cartcov     - 6x6 cartesian covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    r           - matrix of partial derivatives
#    a           - semimajor axis                 km
#    ecc         - eccentricity
#    incl        - inclination                    0.0  to pi rad
#    omaga       - longitude of ascending node    0.0  to 2pi rad
#    argp        - argument of perigee            0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#    m           - mean anomaly                   0.0  to 2pi rad
#    p1, p2, p3, p4 - denominator terms for the partials
#    e0          - eccentric anomaly              0.0  to 2pi rad
#    true1, true2- temp true anomaly              0.0  to 2pi rad
#
#  coupling      :
#    newtonm     - newton iteration for m and ecc to nu
#    newtonnu    - newton iteration for nu and ecc to m
#
#  references    :
#    Vallado and Alfano 2015
#
#   [cartcov, tm] = covcl2ct(classcov, classstate, anom)
# ----------------------------------------------------------------------------

def covcl2ctnew(classcov, classstate, anom):
    # -------- define gravitational constant

    # --------- determine which set of variables is in use ---------
    # ---- parse the input vector into the classical elements -----
    a = classstate[0]

    n = np.sqrt(mum / a ** 3)
    ecc = classstate[1]
    incl = classstate[2]
    raan = classstate[3]
    argp = classstate[4]
    # -------- if mean anomaly is used, convert to true anomaly
    # -------- eccentric anomaly (e) is needed for both
    if (anom == 'meana') or (anom == 'meann'):
        mean = classstate[5]
        e, nu = smu.newtonm(ecc, mean)
    elif (anom == 'truea') or (anom == 'truen'):
        # note that mean is not used in the partials, but nu is!
        nu = classstate[5]
        e, mean = smu.newtonnu(ecc, nu)

    p = a * (1 - ecc ** 2) / 1000

    # r, v = coe2rv(p, ecc, incl, raan, argp, nu, 0.0, 0.0, 0.0)
    # rx = r[0] * 1000
    # ry = r[1] * 1000
    # rz = r[2] * 1000
    # vx = v[0] * 1000
    # vy = v[1] * 1000
    # vz = v[2] * 1000

    # assign trig values for efficiency
    # sin_inc = np.sin(incl)
    # cos_inc = np.cos(incl)
    # sin_raan = np.sin(raan)
    # cos_raan = np.cos(raan)
    # sin_w = np.sin(argp)
    # cos_w = np.cos(argp)
    # sin_nu = np.sin(nu)
    # cos_nu = np.cos(nu)
    sin_inc, cos_inc, sin_raan, cos_raan, sin_w, cos_w, sin_nu, cos_nu = \
        smu.getsincos(incl, raan, argp, nu)

    # assign elements of PQW to ECI transformation (pg 168)
    p11 = cos_raan * cos_w - sin_raan * sin_w * cos_inc
    p12 = -cos_raan * sin_w - sin_raan * cos_w * cos_inc
    p13 = sin_raan * sin_inc
    p21 = sin_raan * cos_w + cos_raan * sin_w * cos_inc
    p22 = -sin_raan * sin_w + cos_raan * cos_w * cos_inc
    p23 = -cos_raan * sin_inc
    p31 = sin_w * sin_inc
    p32 = cos_w * sin_inc
    # p33 = cos_inc
    # assign constants for efficiency
    p0 = math.sqrt(mum / (a * (1.0 - ecc * ecc)))
    p1 = (1.0 - ecc * ecc) / (1.0 + ecc * cos_nu)
    p2 = 1.0 / (2.0 * a) * p0
    p3 = ((2.0 * a * ecc + a * cos_nu + a * cos_nu * ecc * ecc)
          / ((1.0 + ecc * cos_nu) ** 2))
    p4 = ecc * mum / (a * (1.0 - ecc * ecc) ** 2 * p0)
    p5 = a * p1
    p6 = a * (1.0 - ecc * ecc) / ((1.0 + ecc * cos_nu) ** 2)
    dMdnu = (1.0 - ecc * ecc) ** 1.5 / ((1.0 + ecc * cos_nu) ** 2)
    dMde = (- sin_nu
            * ((ecc * cos_nu + 1) * (ecc + cos_nu) / np.sqrt((ecc + cos_nu) ** 2)
               + 1.0 - 2.0 * ecc ** 2 - ecc ** 3 * np.cos(nu))
            / ((ecc * np.cos(nu) + 1.0) ** 2 * np.sqrt(1 - ecc ** 2)))

    # does the sign() work for all cases? reduced form from Fabian
    # dMde1 = -sin(nu)*((ecc*cos(nu) + 1)*sign(ecc+cos(nu)) + 1.0 - 2.0*ecc**2 - ecc**3*cos(nu)) / ((ecc*cos(nu) + 1.0)**2 * sqrt(1-ecc**2))  # dm/de
    # dMdnu = (1.0 - ecc**2)**1.5 / ((1.0 + ecc*cos(nu))**2)  # dm/dv
    # dMde = -sin(nu)*((ecc*cos(nu) + 1)*(ecc+cos(nu))/sqrt((ecc + cos(nu))**2) + 1.0 - 2.0*ecc**2 - ecc**3*cos(nu)) / ((ecc*cos(nu) + 1.0)**2 * sqrt(1-ecc**2))  # dm/de

    # ---------------- calculate matrix elements ------------------
    # ---- partials of (a ecc incl node argp nu) wrt rx
    # same
    tm = np.zeros((6, 6))
    tm[0, 0] = p1 * (p11 * cos_nu + p12 * sin_nu)
    #tm[0, 0] = rx/a # alternate approach if vectors available
    tm[0, 1] = - p3 * (p11 * cos_nu + p12 * sin_nu)
    tm[0, 2] = p5 * p13 * (sin_w * cos_nu + cos_w * sin_nu)
    tm[0, 4] = - p5 * (p21 * cos_nu + p22 * sin_nu)
    tm[0, 5] = p5 * (p12 * cos_nu - p11 * sin_nu)
    # true anomaly same
    # p10 = a * (ecc**2 - 1.0) / (ecc*cos_nu + 1.0)**2
    # tm[0, 5] = p10 * (ecc*cos(raan)*sin(argp) + cos(raan)*cos(argp)*sin(nu) + ...
    #           cos(raan)*sin(argp)*cos_nu + ecc*cos(incl)*sin(raan)*cos(argp) + ...
    #           cos(incl)*sin(raan)*cos(argp)*cos(nu) - cos(incl)*sin(raan)*sin(argp)*sin(nu))
    tm[0, 5] = p6 * (- p11 * sin_nu + p12 * (ecc + cos_nu))
    #atm = tm[0, 5]
    if (anom == 'meana') or (anom == 'meann'):
        #tm[0, 5] = tm[0, 5] / dMdnu + tm[0, 1] / dMde
        #tm[0, 1] = tm[0, 5] * dMde - dMde*atm/dMdnu
        tm[0, 5] = tm[0, 5] / dMdnu
        tm[0, 1] = tm[0, 1] - tm[0, 5] * dMde

    # ---- partials of (a ecc incl node argp nu) wrt ry
    tm[1, 0] = p1 * (p21 * cos_nu + p22 * sin_nu)
    #tm[1, 0] = ry/a
    tm[1, 1] = - p3 * (p21 * cos_nu + p22 * sin_nu)
    tm[1, 2] = p5 * p23 * (sin_w * cos_nu + cos_w * sin_nu)
    tm[1, 4] = p5 * (p11 * cos_nu + p12 * sin_nu)
    tm[1, 5] = p5 * (p22 * cos_nu - p21 * sin_nu)
    # true anomaly, same
    #  p10 = a * (ecc**2 - 1.0) / (ecc*cos(nu) + 1.0)**2
    #  tm[1, 5] = p10 * (ecc*sin(raan)*sin(argp) + sin(raan)*cos(argp)*sin(nu) + ...
    #          sin(raan)*sin(argp)*cos(nu) - ecc*cos(incl)*cos(raan)*cos(argp) - ...
    #           cos(incl)*cos(raan)*cos(argp)*cos(nu)+ cos(incl)*cos(raan)*sin(argp)*sin(nu))
    tm[1, 5] = p6 * (- p21 * sin_nu + p22 * (ecc + cos_nu))
    #atm = tm[1, 5]
    if ((anom == 'meana') or (anom == 'meann')):
        #tm[1, 5] = tm[1, 5] / dMdnu + tm[1, 1] / dMde
        #tm[1, 1] = tm[1, 5] * dMde - dMde*atm/dMdnu
        tm[1, 5] = tm[1, 5] / dMdnu
        tm[1, 1] = tm[1, 1] - tm[1, 5] * dMde

    # ---- partials of (a ecc incl node argp nu) wrt rz
    tm[2, 0] = p1 * (p31 * cos_nu + p32 * sin_nu)
    #tm[2, 0] = rz/a
    tm[2, 1] = - p3 * sin_inc * (cos_w * sin_nu + sin_w * cos_nu)
    tm[2, 2] = p5 * cos_inc * (cos_w * sin_nu + sin_w * cos_nu)
    tm[2, 4] = 0.0
    tm[2, 5] = p5 * sin_inc * (cos_w * cos_nu - sin_w * sin_nu)
    #  p10 = -a * (ecc**2 - 1.0) / (ecc*cos(nu) + 1.0)**2
    #  tm[2, 5] = p10 * sin(incl)*(cos(argp+nu)+ecc*cos(argp))
    tm[2, 5] = p6 * (- p31 * sin_nu + p32 * (ecc + cos_nu))
    #atm = tm[2, 5]
    if (anom == 'meana') or (anom == 'meann'):
        #tm[2, 5] = tm[2, 5] / dMdnu + tm[2, 1] / dMde
        #tm[2, 1] = tm[2, 5] * dMde - dMde*atm/dMdnu
        tm[2, 5] = tm[2, 5] / dMdnu
        tm[2, 1] = tm[2, 1] - tm[2, 5] * dMde

    # ---- partials of (a ecc incl node argp nu) wrt vx
    tm[3, 0] = p2 * (p11 * sin_nu - p12 * (ecc + cos_nu))
    #tm[3, 0] = -vx/(2.0*a)
    tm[3, 1] = - p4 * (p11 * sin_nu - p12 * (ecc + cos_nu)) + p12 * p0
    tm[3, 2] = - p0 * sin_raan * (p31 * sin_nu - p32 * (ecc + cos_nu))
    tm[3, 3] = p0 * (p21 * sin_nu - p22 * (ecc + cos_nu))
    tm[3, 4] = - p0 * (p12 * sin_nu + p11 * (ecc + cos_nu))
    # same
    # p10 = sqrt(-mum/(a*(ecc**2-1.0)))
    # tm[3, 5] = p10 * (cos(raan)*sin(argp)*sin(nu)-cos(raan)*cos(argp)*cos(nu)+cos(incl)*sin(raan)*cos(argp)*sin(nu)+cos(incl)*sin(raan)*sin(argp)*cos(nu))
    tm[3, 5] = - p0 * (p11 * cos_nu + p12 * sin_nu)
    #atm = tm[3, 5]
    if (anom == 'meana') or (anom == 'meann'):
        #tm[3, 5] = tm[3, 5] / dMdnu + tm[3, 1] / dMde
        #tm[3, 1] = tm[3, 5] * dMde - dMde*atm/dMdnu
        tm[3, 5] = tm[3, 5] / dMdnu
        tm[3, 2] = tm[3, 2] - tm[3, 5] * dMde

    # ---- partials of (a ecc incl node argp nu) wrt vy
    tm[4, 0] = p2 * (p21 * sin_nu - p22 * (ecc + cos_nu))
    #tm[4, 0] = -vy/(2.0*a)
    tm[4, 1] = - p4 * (p21 * sin_nu - p22 * (ecc + cos_nu)) + p22 * p0
    tm[4, 2] = p0 * cos_raan * (p31 * sin_nu - p32 * (ecc + cos_nu))
    tm[4, 3] = p0 * (- p11 * sin_nu + p12 * (ecc + cos_nu))
    tm[4, 4] = - p0 * (p22 * sin_nu + p21 * (ecc + cos_nu))
    #  p10 = sqrt(-mum/(a*(ecc**2-1.0)))
    #  tm[4, 5] = -p10 * (sin(raan)*cos(argp)*cos(nu)-sin(raan)*sin(argp)*sin(nu)+cos(incl)*cos(raan)*cos(argp)*sin(nu)+cos(incl)*cos(raan)*sin(argp)*cos(nu))
    tm[4, 5] = - p0 * (p21 * cos_nu + p22 * sin_nu)
    #atm = tm[4, 5]
    if (anom == 'meana') or (anom == 'meann'):
        #tm[4, 5] = tm[4, 5] / dMdnu + tm[4, 1] / dMde
    #tm[4, 1] = tm[4, 5] * dMde - dMde*atm/dMdnu
        tm[4, 5] = tm[4, 5] / dMdnu
        tm[4, 1] = tm[4, 1] - tm[4, 5] * dMde

    # ---- partials of (a ecc incl node argp nu) wrt vz
    # same
    tm[5, 0] = p2 * (p31 * sin_nu - p32 * (ecc + cos_nu))
    # tm[5, 0] = -vz/(2.0*a)
    tm[5, 1] = - p4 * (p31 * sin_nu - p32 * (ecc + cos_nu)) + p32 * p0
    tm[5, 2] = p0 * cos_inc * (cos_w * cos_nu - sin_w * sin_nu + ecc * cos_w)
    tm[5, 3] = 0.0
    tm[5, 4] = - p0 * (p32 * sin_nu + p31 * (ecc + cos_nu))
    # same
    # p10 = sqrt(-mum/(a*(ecc**2-1.0)))
    # atm[5, 5] = p10 * (-sin(incl)*sin(argp+nu))
    tm[5, 5] = - p0 * (p31 * cos_nu + p32 * sin_nu)
    #atm = tm[5, 5]
    if (anom == 'meana') or (anom == 'meann'):
        #tm[5, 5] = tm[5, 5] / dMdnu + tm[5, 1] / dMde
        #tm[5, 1] = tm[5, 5] * dMde - dMde*atm/dMdnu
        tm[5, 5] = tm[5, 5] / dMdnu
        tm[5, 1] = tm[5, 1] - tm[5, 5] * dMde

    # ---------- calculate the output covariance matrix -----------
    cartcov = tm @ classcov @ tm.T
    return cartcov, tm

# ----------------------------------------------------------------------------
#
#                           function covct2cl
#
#  this function transforms a six by six covariance matrix expressed in cartesian elements
#    into one expressed in classical elements
#
#  author        : david vallado                  719-573-2600   21 jun 2002
#
#  revisions
#    vallado     - major update                                  26 aug 2015
#
#  inputs          description                    range / units
#    cartcov     - 6x6 cartesian covariance matrix
#    cartstate   - 6x1 cartesian orbit state      (rx ry rz vx vy vz)
#    anom        - anomaly                        'meana', 'truea', 'meann', 'truen'
#
#  outputs       :
#    classcov    - 6x6 classical covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    r           - matrix of partial derivatives
#    rj2000      - position vector                km
#    x, y, z       - components of position vector  km
#    vj2000      - velocity vector                km/s
#    vx, vy, vz    - components of position vector  km/s
#    p           - semilatus rectum               km
#    a           - semimajor axis                 km
#    ecc         - eccentricity
#    incl        - inclination                    0.0  to pi rad
#    omaga       - longitude of ascending node    0.0  to 2pi rad
#    argp        - argument of perigee            0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#    m           - mean anomaly                   0.0  to 2pi rad
#    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
#    truelon     - true longitude            (ce) 0.0  to 2pi rad
#    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
#    magr        - magnitude of position vector   km
#    magv        - magnitude of velocity vector   km/s
#
#  coupling      :
#    rv2coe      - position and velocity vectors to classical elements
#
#  references    :
#    Vallado and Alfano 2015
#
#   [classcov, tm] = covct2cl(cartcov, cartstate, anom)
# ----------------------------------------------------------------------------

def covct2clnew(cartcov, cartstate, anom):
    # -------- parse the input vectors into cartesian and classical components
    rx = cartstate[0, 0] * 1000.0
    ry = cartstate[1, 0] * 1000.0
    rz = cartstate[2, 0] * 1000.0
    vx = cartstate[3, 0] * 1000.0
    vy = cartstate[4, 0] * 1000.0
    vz = cartstate[5, 0] * 1000.0

    reci = np.array([rx / 1000, ry / 1000, rz / 1000])
    veci = np.array([vx / 1000, vy / 1000, vz / 1000])

    # -------- convert to a classical orbit state for ease of computation
    p, a, ecc, _, _, _, nu, _, _, _, _ = rv2coe(reci, veci)
    p = p * 1000.0
    a = a * 1000.0
    n = math.sqrt(mum / a ** 3)
    # -------- calculate common quantities
    # sqrt1me2 = np.sqrt(1.0 - ecc * ecc)
    magr = math.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
    magr3 = magr ** 3
    magv = math.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    # ----------  form pqw position and velocity vectors ----------
    r_dot_v = np.dot(reci, veci) * 1000 * 1000
    ecc_term = magv * magv - mum / magr
    ecc_x = (ecc_term * rx - r_dot_v * vx) / mum
    ecc_y = (ecc_term * ry - r_dot_v * vy) / mum
    ecc_z = (ecc_term * rz - r_dot_v * vz) / mum
    ecc_vec = np.array([ecc_x, ecc_y, ecc_z]).T
    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    h_vec = np.array([hx, hy, hz]).T
    h = smu.mag(h_vec)
    h_squared = h * h
    nx = -hy
    ny = hx
    nz = 0.0
    node_vec = np.array([nx, ny, nz]).T
    node = smu.mag(node_vec)
    n_squared = node * node
    n_dot_e = np.dot(node_vec, ecc_vec)
    # sign_anode = np.sign(ny)
    # cos_anode = nx / node
    # omega = sign_anode * np.arccos(cos_anode)
    sign_w = np.sign((magv ** 2 - mum / magr) * rz - r_dot_v * vz)
    cos_w = n_dot_e / (ecc * node)
    # argp = sign_w * np.arccos(cos_w)
    w_scale = - sign_w / np.sqrt(1 - cos_w * cos_w)
    r_dot_e = np.dot(reci, ecc_vec) * 1000
    cos_nu = r_dot_e / (magr * ecc)
    sign_nu = np.sign(r_dot_v)
    nu = sign_nu * np.arccos(cos_nu)
    nu_scale = - sign_nu / np.sqrt(1 - cos_nu * cos_nu)
    _, ax, eccx, inclx, nodex, argpx, nux, mx, _, _, _ = rv2coe(reci, veci)
    print(' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f \n'
          % (ax, eccx, inclx * rad2deg, nodex * rad2deg, argpx * rad2deg,
             nux * rad2deg, mx * rad2deg))
    # ---------------- calculate matrix elements ------------------
    # ---- partials of a wrt (rx ry rz vx vy vz)
    p0 = 2.0 * a ** 2 / magr ** 3
    p1 = 2.0 / (n ** 2 * a)
    tm = np.zeros((6, 6))
    tm[0, 0] = p0 * rx
    tm[0, 1] = p0 * ry
    tm[0, 2] = p0 * rz
    tm[0, 3] = p1 * vx
    tm[0, 4] = p1 * vy
    tm[0, 5] = p1 * vz
    # ---- partials of ecc wrt (rx ry rz vx vy vz)
    p0 = 1.0 / (mum * ecc)
    tm[1, 0] = -p0 * (((vx * vy - mum * rx * ry / magr3) * ecc_y)
                      + ((vx * vz - mum * rx * rz / magr3) * ecc_z)
                      - (vy * vy + vz * vz - mum / magr + mum * rx * rx / magr3)
                      * ecc_x)
    tm[1, 1] = -p0 * (((vx * vy - mum * rx * ry / magr3) * ecc_x)
                     + ((vy * vz - mum * ry * rz / magr3) * ecc_z)
                     - (vx * vx + vz * vz - mum / magr + mum * ry * ry / magr3)
                     * ecc_y)
    tm[1, 2] = -p0 * (((vx * vz - mum * rx * rz / magr3) * ecc_x)
                     + ((vy * vz - mum * ry * rz / magr3) * ecc_y)
                     - (vy * vy + vx * vx - mum / magr + mum * rz * rz / magr3)
                     * ecc_z)
    tm[1, 3] = -p0 * ((rx * vy - 2 * ry * vx) * ecc_y
                     + (ry * vy + rz * vz) * ecc_x
                     + (rx * vz - 2 * rz * vx) * ecc_z)
    tm[1, 4] = -p0 * ((ry * vx - 2 * rx * vy) * ecc_x
                     + (rx * vx + rz * vz) * ecc_y
                     + (ry * vz - 2 * rz * vy) * ecc_z)
    tm[1, 5] = -p0 * ((rx * vx + ry * vy) * ecc_z
                      + (rz * vx - 2 * rx * vz) * ecc_x
                      + (rz * vy - 2 * ry * vz) * ecc_y)
    # ---- partials of incl wrt (rx ry rz vx vy vz)
    p3 = 1.0 / node
    tm[2, 0] = -p3 * (vy - hz * (vy * hz - vz * hy) / h_squared)
    tm[2, 1] = p3 * (vx - hz * (vx * hz - vz * hx) / h_squared)
    tm[2, 2] = -p3 * (hz * (vy * hx - vx * hy) / h_squared)
    tm[2, 3] = p3 * (ry - hz * (ry * hz - rz * hy) / h_squared)
    tm[2, 4] = -p3 * (rx - hz * (rx * hz - rz * hx) / h_squared)
    tm[2, 5] = p3 * (hz * (ry * hx - rx * hy) / h_squared)

    # ---- partials of node wrt (rx ry rz vx vy vz)
    p4 = 1.0 / n_squared
    tm[3, 0] = -p4 * vz * ny
    tm[3, 1] = p4 * vz * nx
    tm[3, 2] = p4 * (vx * ny - vy * nx)
    tm[3, 3] = p4 * rz * ny
    tm[3, 4] = -p4 * rz * nx
    tm[3, 5] = p4 * (ry * nx - rx * ny)
    # ---- partials of argp wrt (rx ry rz vx vy vz)
    # p5 = 1.0 / (node * a * a)
    temp = - hy * (vy * vy + vz * vz - mum / magr + mum * rx * rx / magr3)
    temp = temp - hx * (vx * vy - mum * rx * ry / magr3) + vz * mum * ecc_x
    temp = (temp / (mum * node * ecc)
            + vz * hy * n_dot_e / (node * node * node * ecc)
            - tm[1, 0] * n_dot_e / (node * ecc * ecc))
    tm[4, 0] = temp * w_scale
    temp = hx * (vx * vx + vz * vz - mum / magr + mum * ry * ry / magr3)
    temp = temp + hy * (vx * vy - mum * rx * ry / magr3) + vz * mum * ecc_y
    temp = (temp / (mum * node * ecc)
            - vz * hx * n_dot_e / (node * node * node * ecc)
            - tm[1, 1] * n_dot_e / (node * ecc * ecc))
    tm[4, 1] = temp * w_scale
    temp = (- hy * (vx * vz - mum * rx * rz / magr3)
            + hx * (vy * vz - mum * ry * rz / magr3)
            + vx * mum * ecc_x
            + vy * mum * ecc_y)
    temp = (- temp / (mum * node * ecc)
            + (vy * hx - vx * hy) * n_dot_e / (node * node * node * ecc)
            - tm[1, 2] * n_dot_e / (node * ecc * ecc))
    tm[4, 2] = temp * w_scale
    temp = ((rx * vy - 2 * ry * vx) * hx - hy * (ry * vy + rz * vz)
            + rz * mum * ecc_x)
    temp = (- temp / (mum * node * ecc)
            - rz * hy * n_dot_e / (node * node * node * ecc)
            - tm[1, 3] * n_dot_e / (node * ecc * ecc))
    tm[4, 3] = temp * w_scale
    temp = (- (ry * vx - 2 * rx * vy) * hy + hx * (rx * vx + rz * vz)
            + rz * mum * ecc_y)
    temp = (- temp / (mum * node * ecc)
            + rz * hx * n_dot_e / (node * node * node * ecc)
            - tm[1, 4] * n_dot_e / (node * ecc * ecc))
    tm[4, 4] = temp * w_scale
    temp = (- (rz * vx - 2 * rx * vz) * hy + hx * (rz * vy - 2 * ry * vz)
            - rx * mum * ecc_x - ry * mum * ecc_y)
    temp = (- temp / (mum * node * ecc)
            + (rx * hy - ry * hx) * n_dot_e / (node * node * node * ecc)
            - tm[1, 5] * n_dot_e / (node * ecc * ecc))
    tm[4, 5] = temp * w_scale
    # ---- partials of nu/M wrt (rx ry rz vx vy vz)
    temp = (ry * (vx * vy - mum * rx * ry / magr3)
            - rx * ecc_term + rz * (vx * vz - mum * rx * rz / magr3))
    temp = (temp
            - rx * (vy * vy + vz * vz - mum / magr + mum * rx * rx / magr3)
            + vx * r_dot_v)
    temp = (- temp / (mum * magr * ecc)
            - rx * r_dot_e / (magr3 * ecc)
            - tm[1, 0] * r_dot_e / (magr * ecc * ecc))
    tm[5, 0] = temp * nu_scale
    temp = (rx * (vx * vy - mum * rx * ry / magr3)
            - ry * ecc_term + rz * (vy * vz - mum * ry * rz / magr3))
    temp = (temp
            - ry * (vx * vx + vz * vz - mum / magr + mum * ry * ry / magr3)
            + vy * r_dot_v)
    temp = (- temp / (mum * magr * ecc)
            - ry * r_dot_e / (magr3 * ecc)
            - tm[1, 1] * r_dot_e / (magr * ecc * ecc))
    tm[5, 1] = temp * nu_scale
    temp = (rx * (vx * vz - mum * rx * rz / magr3)
            - rz * ecc_term + ry * (vy * vz - mum * ry * rz / magr3))
    temp = (temp
            - rz * (vx * vx + vy * vy - mum / magr + mum * rz * rz / magr3)
            + vz * r_dot_v)
    temp = (- temp / (mum * magr * ecc)
            - rz * r_dot_e / (magr3 * ecc)
            - tm[1, 2] * r_dot_e / (magr * ecc * ecc))
    tm[5, 2] = temp * nu_scale
    temp = (ry * (rx * vy - 2 * ry * vx)
            + rx * (ry * vy + rz * vz)
            + rz * (rx * vz - 2 * rz * vx))
    temp = - temp / (mum * magr * ecc) - tm[1, 3] * r_dot_e / (magr * ecc * ecc)
    tm[5, 3] = temp * nu_scale
    temp = (rx * (ry * vx - 2 * rx * vy)
            + ry * (rx * vx + rz * vz) + rz * (ry * vz - 2 * rz * vy))
    temp = - temp / (mum * magr * ecc) - tm[1, 4] * r_dot_e / (magr * ecc * ecc)
    tm[5, 4] = temp * nu_scale
    temp = (rz * (rx * vx + ry * vy)
            + rx * (rz * vx - 2 * rx * vz) + ry * (rz * vy - 2 * ry * vz))
    temp = - temp / (mum * magr * ecc) - tm[1, 5] * r_dot_e / (magr * ecc * ecc)
    tm[5, 5] = temp * nu_scale
    #       # same answers as above
#       # ---- partials of (rx ry rz vx vy vz) wrt true anomaly
#       p8 = magr**2*magv**2 - mum*magr - r_dot_v**2
#       p9 = 1.0/(p8**2 + r_dot_v**2 * h**2)
#       tm[5, 0] = p9 * (p8 * (h*vx + r_dot_v*(vy*hz - vz*hy)/h) - r_dot_v*h*(2*rx*magv**2 - mum*rx/magr - 2*r_dot_v*vx))
#       tm[5, 1] = p9 * (p8 * (h*vy + r_dot_v*(vz*hx - vx*hz)/h) - r_dot_v*h*(2*ry*magv**2 - mum*ry/magr - 2*r_dot_v*vy))
#       tm[5, 2] = p9 * (p8 * (h*vz + r_dot_v*(vx*hy - vy*hx)/h) - r_dot_v*h*(2*rz*magv**2 - mum*rz/magr - 2*r_dot_v*vz))
#       tm[5, 3] = p9 * (p8 * (h*rx + r_dot_v*(rz*hy - ry*hz)/h) - r_dot_v*h*(2*vx*magr**2 - 2*r_dot_v*rx))
#       tm[5, 4] = p9 * (p8 * (h*ry + r_dot_v*(rx*hz - rz*hx)/h) - r_dot_v*h*(2*vy*magr**2 - 2*r_dot_v*ry))
#       tm[5, 5] = p9 * (p8 * (h*rz + r_dot_v*(ry*hx - rx*hy)/h) - r_dot_v*h*(2*vz*magr**2 - 2*r_dot_v*rz))

    if (anom == 'meana') or (anom == 'meann'):
        # ---- partials of (rx ry rz vx vy vz) wrt mean anomaly
        # then update for mean anomaly
        ecc = smu.mag(ecc_vec)
        dMdnu = (1.0 - ecc ** 2) ** 1.5 / ((1.0 + ecc * math.cos(nu)) ** 2)
        dMde = (-math.sin(nu) * ((ecc * math.cos(nu) + 1)
                               * (ecc + math.cos(nu))
                               / math.sqrt((ecc + math.cos(nu)) ** 2)
                               + 1.0 - 2.0 * ecc ** 2 - ecc ** 3 * math.cos(nu))
                / ((ecc * math.cos(nu) + 1.0) ** 2 * math.sqrt(1 - ecc ** 2)))
        # p6 = -sin(nu)*(sign((ecc*cos(nu) + 1)) + 1.0 - 2.0*ecc**2 - ecc**3*cos(nu)) / ((ecc*cos(nu) + 1.0)**2 * sqrt(1-ecc**2))  # dm/de
        tm[5, 0] = tm[5, 0] * dMdnu + tm[1, 0] * dMde
        tm[5, 1] = tm[5, 1] * dMdnu + tm[1, 1] * dMde
        tm[5, 2] = tm[5, 2] * dMdnu + tm[1, 2] * dMde
        tm[5, 3] = tm[5, 3] * dMdnu + tm[1, 3] * dMde
        tm[5, 4] = tm[5, 4] * dMdnu + tm[1, 4] * dMde
        tm[5, 5] = tm[5, 5] * dMdnu + tm[1, 5] * dMde

    # ---------- calculate the output covariance matrix -----------
    classcov = tm @ cartcov @ tm.T
    return classcov, tm

#
# ----------------------------------------------------------------------------
#
#                           function covct2ntw
#
#  this function transforms a six by six covariance matrix expressed in cartesian
#    into one expressed in orbit plane, ntw frame
#
#  author        : david vallado                  719-573-2600   20 may 2003
#
#  revisions
#    vallado     - fix indices                                   16 jul 2003
#    vallado     - send out tm                                   21 jul 2003
#
#  inputs          description                    range / units
#    cartcov     - 6x6 cartesian covariance matrix
#    cartstate   - 6x1 cartesian orbit state      (x y z vx vy vz)
#
#  outputs       :
#    covntw      - 6x6 orbit plane ntw covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    r           - position vector                m
#    v           - velocity vector                m/s
#    temv        - temporary vector
#
#  coupling      :
#    none
#
#  references    :
#    none
#
#  [covopntw, tm] = covct2ntw(cartcov, cartstate)
# ----------------------------------------------------------------------------

def covct2ntw(cartcov, cartstate):
    x = cartstate[0, 0]
    y = cartstate[1, 0]
    z = cartstate[2, 0]
    vx = cartstate[3, 0]
    vy = cartstate[4, 0]
    vz = cartstate[5, 0]
    r = np.array([x, y, z])
    v = np.array([vx, vy, vz])
    tv = smu.unit(v)
    temv = np.cross(r, v)
    wv = smu.unit(temv)
    nv = np.cross(tv, wv)
    tm = np.zeros((6, 6))

    tm[0, 0] = nv[0]
    tm[0, 1] = nv[1]
    tm[0, 2] = nv[2]
    tm[1, 0] = tv[0]
    tm[1, 1] = tv[1]
    tm[1, 2] = tv[2]
    tm[2, 0] = wv[0]
    tm[2, 1] = wv[1]
    tm[2, 2] = wv[2]
    tm[3, 3] = nv[0]
    tm[3, 4] = nv[1]
    tm[3, 5] = nv[2]
    tm[4, 3] = tv[0]
    tm[4, 4] = tv[1]
    tm[4, 5] = tv[2]
    tm[5, 3] = wv[0]
    tm[5, 4] = wv[1]
    tm[5, 5] = wv[2]
    covntw = tm @ cartcov @ tm.T
    return covntw, tm


#
# ----------------------------------------------------------------------------
#
#                           function covct2o2
#
#  this function transforms a six by six covariance matrix expressed in cartesian
#    into one expressed in orbit plane, rsw frame
#
#  author        : david vallado                  719-573-2600   20 may 2003
#
#  revisions
#    vallado     - fix indices                                   16 jul 2003
#    vallado     - send out tm                                   21 jul 2003
#
#  inputs          description                    range / units
#    cartcov     - 6x6 cartesian covariance matrix
#    cartstate   - 6x1 cartesian orbit state      (x y z vx vy vz)
#
#  outputs       :
#    covoprsw    - 6x6 orbit plane rsw covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    r           - position vector                m
#    v           - velocity vector                m/s
#    temv        - temporary vector
#
#  coupling      :
#    none
#
#  references    :
#    none
#
#  [covopntw, tm] = covct2o2(cartcov, cartstate)
# ----------------------------------------------------------------------------

def covct2o2(cartcov, cartstate):
    x = cartstate[0, 0]
    y = cartstate[1, 0]
    z = cartstate[2, 0]
    vx = cartstate[3, 0]
    vy = cartstate[4, 0]
    vz = cartstate[5, 0]
    r = np.array([x, y, z])
    v = np.array([vx, vy, vz])
    tv = smu.unit(v)
    temv = np.cross(r, v)
    wv = smu.unit(temv)
    nv = np.cross(tv, wv)
    tm = np.zeros((6, 6))

    tm[0, 0] = nv[0]
    tm[0, 1] = nv[1]
    tm[0, 2] = nv[2]
    tm[1, 0] = tv[0]
    tm[1, 1] = tv[1]
    tm[1, 2] = tv[2]
    tm[2, 0] = wv[0]
    tm[2, 1] = wv[1]
    tm[2, 2] = wv[2]
    tm[3, 3] = nv[0]
    tm[3, 4] = nv[1]
    tm[3, 5] = nv[2]
    tm[4, 3] = tv[0]
    tm[4, 4] = tv[1]
    tm[4, 5] = tv[2]
    tm[5, 3] = wv[0]
    tm[5, 4] = wv[1]
    tm[5, 5] = wv[2]
    covopntw = tm @ cartcov @ tm.T
    return covopntw, tm


#
# ----------------------------------------------------------------------------
#
#                           function covct2rsw
#
#  this function transforms a six by six covariance matrix expressed in cartesian
#    into one expressed in orbit plane, rsw frame
#
#  author        : david vallado                  719-573-2600   20 may 2003
#
#  revisions
#    vallado     - fix indices                                   16 jul 2003
#    vallado     - send out tm                                   21 jul 2003
#
#  inputs          description                    range / units
#    cartcov     - 6x6 cartesian covariance matrix
#    cartstate   - 6x1 cartesian orbit state      (x y z vx vy vz)
#
#  outputs       :
#    covoprsw    - 6x6 orbit plane rsw covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    r           - position vector                m
#    v           - velocity vector                m/s
#    temv        - temporary vector
#
#  coupling      :
#    none
#
#  references    :
#    none
#
#  [covoprsw, tm] = covct2rsw(cartcov, cartstate)
# ----------------------------------------------------------------------------

def covct2rsw(cartcov, cartstate):
    print("cartstate:")
    print(cartstate)
    x = cartstate[0, 0]
    y = cartstate[1, 0]
    z = cartstate[2, 0]
    vx = cartstate[3, 0]
    vy = cartstate[4, 0]
    vz = cartstate[5, 0]
    r = np.array([x, y, z])
    v = np.array([vx, vy, vz])
    rv = smu.unit(r)
    temv = np.cross(r, v)
    wv = smu.unit(temv)
    sv = np.cross(wv, rv)
    tm = np.zeros((6, 6))

    tm[0, 0] = rv[0]
    tm[0, 1] = rv[1]
    tm[0, 2] = rv[2]
    tm[1, 0] = sv[0]
    tm[1, 1] = sv[1]
    tm[1, 2] = sv[2]
    tm[2, 0] = wv[0]
    tm[2, 1] = wv[1]
    tm[2, 2] = wv[2]
    tm[3, 3] = rv[0]
    tm[3, 4] = rv[1]
    tm[3, 5] = rv[2]
    tm[4, 3] = sv[0]
    tm[4, 4] = sv[1]
    tm[4, 5] = sv[2]
    tm[5, 3] = wv[0]
    tm[5, 4] = wv[1]
    tm[5, 5] = wv[2]
    covoprsw = tm @ cartcov @ tm.T
    return covoprsw, tm


#
# ----------------------------------------------------------------------------
#
#                           function covo22ct
#
#  this function transforms a six by six covariance matrix expressed in the
#    orbit plane (ntw) into one expressed in cartesian
#
#  author        : david vallado                  719-573-2600   17 jul 2003
#
#  revisions
#    vallado     - send out tm                                   25 jul 2003
#
#  inputs          description                    range / units
#    covopntw    - 6x6 orbit plane ntw covariance matrix
#    cartstate   - 6x1 cartesian orbit state      (x y z vx vy vz)
#
#  outputs       :
#    cartcov     - 6x6 cartesian covariance matrix
#    tm          - transformation matrix
#
#  locals        :
#    r           - position vector                m
#    v           - velocity vector                m/s
#    temv        - temporary vector
#
#  coupling      :
#    none
#
#  references    :
#    none
#
#  [cartcov, tm] = covo22ct(covopntw, cartstate)
# ----------------------------------------------------------------------------

def covo22ct(covopntw, cartstate):
    x = cartstate[0, 0]
    y = cartstate[1, 0]
    z = cartstate[2, 0]
    vx = cartstate[3, 0]
    vy = cartstate[4, 0]
    vz = cartstate[5, 0]
    r = np.array([x, y, z])
    v = np.array([vx, vy, vz])
    tv = smu.unit(v)
    temv = np.cross(r, v)
    wv = smu.unit(temv)
    nv = np.cross(tv, wv)
    tm = np.zeros((6, 6))

    tm[0, 0] = nv[0]
    tm[0, 1] = nv[1]
    tm[0, 2] = nv[2]
    tm[1, 0] = tv[0]
    tm[1, 1] = tv[1]
    tm[1, 2] = tv[2]
    tm[2, 0] = wv[0]
    tm[2, 1] = wv[1]
    tm[2, 2] = wv[2]
    tm[3, 3] = nv[0]
    tm[3, 4] = nv[1]
    tm[3, 5] = nv[2]
    tm[4, 3] = tv[0]
    tm[4, 4] = tv[1]
    tm[4, 5] = tv[2]
    tm[5, 3] = wv[0]
    tm[5, 4] = wv[1]
    tm[5, 5] = wv[2]
    tm = tm.T
    cartcov = tm @ covopntw @ tm.T
    return cartcov, tm


# ----------------------------------------------------------------------------
#
#                           function adbar2rv.m
#
#  this function transforms the adbarv elements (rtasc, decl, fpa, azimuth,
#    position and velocity magnitude) into eci position and velocity vectors.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    magr        - eci position vector magnitude  km
#    magv        - eci velocity vector magnitude  km/sec
#    rtasc       - right ascension of sateillite  rad
#    decl        - declination of satellite       rad
#    fpav        - sat flight path angle from vertrad
#    az          - sat flight path azimuth        rad
#
#  outputs       :
#    reci        - eci position vector            km
#    veci        - eci velocity vector            km/s
#
#  locals        :
#    none        -
#
#  coupling      :
#    none        -
#
#  references    :
#    vallado       2001, xx
#    chobotov            70
#
# [r, v] = adbar2rv (rmag, vmag, rtasc, decl, fpav, az)
# ----------------------------------------------------------------------------

def adbar2rv(magr: float, magv: float, rtasc: float, decl: float, fpav:float,
             az: float):
    """this function transforms the adbarv elements (rtasc, decl, fpa, azimuth,
    position and velocity magnitude) into eci position and velocity vectors.

    Parameters
    ----------
    magr : float
        eci position vector magnitude: km
    magv : float
        eci velocity vector magnitude: km/s
    rtasc : float
        right ascension of satellite: rad
    decl : float
        declination of satellite: rad
    fpav : float
        satellite flight path from vert: rad
    az : float
        satellite flight path azimuth: rad

    Returns
    -------
    reci: ndarray
        eci position vector: km
    veci: ndarray
        eci velocity vector: km
    """
    # -------- form position vector
    reci = np.zeros(3)
    reci[0] = magr * math.cos(decl) * math.cos(rtasc)
    reci[1] = magr * math.cos(decl) * math.sin(rtasc)
    reci[2] = magr * math.sin(decl)
    # -------- form velocity vector
    veci = np.zeros(3)
    veci[0] = magv * (math.cos(rtasc) * (- math.cos(az) * math.sin(fpav) * math.sin(decl)
                                    + math.cos(fpav) * math.cos(decl))
                   - math.sin(az) * math.sin(fpav) * math.sin(rtasc))
    veci[1] = magv * (math.sin(rtasc) * (- math.cos(az) * math.sin(fpav) * math.sin(decl)
                                    + math.cos(fpav) * math.cos(decl))
                   + math.sin(az) * math.sin(fpav) * math.cos(rtasc))
    veci[2] = magv * (math.cos(az) * math.cos(decl) * math.sin(fpav)
                   + math.cos(fpav) * math.sin(decl))
    return reci, veci



# azl2radc
#
# this function finds the rtasc decl values given the az-el
#
#
#

def azel2radec(az: float, el: float, lat: float, lst: float):
    """This function finds the right ascension and declination values
    given the azimuth and elevation

    Parameters
    ----------
    az : float
        azimuth: rad
    el : float
        elevation: rad
    lat : float
        latitude: rad
    lst : float
        local mean sidreal time

    Returns
    -------
    rtasc : float
        right ascension: rad
    decl : float
        declination: rad
    """
    decl = math.asin(math.sin(el) * math.sin(lat)
                     + math.cos(el) * math.cos(lat) * math.cos(az))
    slha1 = (-(math.sin(az) * math.cos(el) * math.cos(lat))
             / (math.cos(decl) * math.cos(lat)))
    clha1 = ((math.sin(el) - math.sin(lat) * math.sin(decl))
             / (math.cos(decl) * math.cos(lat)))
    lha1 = math.atan2(slha1, clha1)
    #    print(' lha1 #13.7f \n', lha1 * rad2deg)

    # alt approach
    slha2 = - (math.sin(az) * math.cos(el)) / (math.cos(decl))
    clha2 = ((math.cos(lat) * math.sin(el) - math.sin(lat) * math.cos(el)
              * math.cos(az)) / (math.cos(decl)))
    lha2 = math.atan2(slha2, clha2)
    #    print(' lha2 #13.7f \n', lha2 * rad2deg)

    rtasc = lst - lha1
    return rtasc, decl

# ------------------------------------------------------------------------------
#
#                           function raz2rvs
#
#  this function converts range, azimuth, and elevation values with slant
#    range and velocity vectors for a satellite from a radar site in the
#    topocentric horizon (sez) system.
#
#  author        : david vallado                  719-573-2600   10 jun 2002
#
#  revisions
#    vallado     - del unnecessary code                          26 aug 2002
#
#  inputs          description                    range / units
#    rho         - satellite range from site      km
#    az          - azimuth                        0.0 to 2pi rad
#    el          - elevation                      -pi/2 to pi/2 rad
#    drho        - range rate                     km / s
#    daz         - azimuth rate                   rad / s
#    del_        - elevation rate                 rad / s
#
#  outputs       :
#    rhovec      - sez satellite range vector     km
#    drhovec     - sez satellite velocity vector  km / s
#
#  locals        :
#    sinel       - variable for sin(el)
#    cosel       - variable for cos(el)
#    sinaz       - variable for sin(az)
#    cosaz       - variable for cos(az)
#    temp        -
#    temp1       -
#
#  coupling      :
#    mag         - magnitude of a vector
#
#  references    :
#    vallado       2001, 250-251, eq 4-4, eq 4-5
#
# [rhosez, drhosez] = raz2rvs (rho, az, el, drho, daz, del)
# ------------------------------------------------------------------------------

def raz2rvs(rho: float, az: float, el: float, drho: float, daz: float,
            del_: float):
    """this function converts range, azimuth, and elevation values with slant
    range and velocity vectors for a satellite from a radar site in the
    topocentric horizon (sez) system.

    Parameters
    ----------
    rho : float
        satellite range from site: km
    az : float
        azimuth: 0 to 2pi rads
    el : float
        elevation: -pi/2 to pi/2 rads
    drho : float
        range rate: km/s
    daz : float
        azimuth rate: km/s
    del_ : float
        elevation rate: km/s

    Returns
    -------
    rhosez : ndarray
        sez satellite range vector: km
    drhosez : ndarray
        sez satellite velocity vector: km/s
    """
    # ----------------------- initialize values -------------------
    sinel = math.sin(el)
    cosel = math.cos(el)
    sinaz = math.sin(az)
    cosaz = math.cos(az)
    # ------------------- form sez range vector -------------------
    rhosez = np.zeros(3)
    rhosez[0] = -rho * cosel * cosaz
    rhosez[1] = rho * cosel * sinaz
    rhosez[2] = rho * sinel
    # ----------------- form sez velocity vector ------------------
    drhosez = np.zeros(3)
    drhosez[0] = (-drho * cosel * cosaz + rhosez[2] * cosaz * del_
                  + rhosez[1] * daz)
    drhosez[1] = (drho * cosel * sinaz - rhosez[2] * sinaz * del_
                  - rhosez[0] * daz)
    drhosez[2] = drho * sinel + rho * cosel * del_
    return rhosez, drhosez


#
# ------------------------------------------------------------------------------
#
#                           function razel2rv
#
#  this function converts range, azimuth, and elevation and their rates to
#    the geocentric equatorial (eci) position and velocity vectors.
#
#  author        : david vallado                  719-573-2600   30 may 2002
#
#  revisions
#    vallado     - add terms for ast calculation                 30 sep 2002
#
#  inputs          description                    range / units
#    rho         - satellite range from site      km
#    az          - azimuth                        0.0 to 2pi rad
#    el          - elevation                      -pi/2 to pi/2 rad
#    drho        - range rate                     km/s
#    daz         - azimuth rate                   rad / s
#    del_        - elevation rate                 rad / s
#    latgd       - geodetic latitude              -pi/2 to pi/2 rad
#    lon         - longitude of site              -2pi to 2pi rad
#    alt         - altitude                       km
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    terms       - number of terms for ast calculation 0, 2
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    reci        - eci position vector            km
#    veci        - eci velocity vector            km/s
#
#  locals        :
#    rs          - ecef site position vector      km
#    rhoecef     - ecef range vector from site    km
#    drhoecef    - ecef velocity vector from site km/s
#    rhosez      - sez range vector from site     km
#    drhosez     - sez velocity vector from site  km
#    tempvec     - temporary vector
#
#  coupling      :
#    raz2rvs     - find r and v from site in topocentric horizon (sez) system
#
#  references    :
#    vallado       2001, 250-255, alg 27
#
# [reci, veci] = razel2rv (rho, az, el, drho, daz, del, latgd, lon, alt, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ------------------------------------------------------------------------------

def razel2rv(rho: float, az: float, el: float, drho: float, daz: float,
             del_: float, latgd: float, lon: float, alt: float, ttt: float,
             jdut1: float, lod: float, xp: float, yp: float, terms: int,
             ddpsi: float, ddeps: float):
    """this function converts range, azimuth, and elevation and their rates to
    the geocentric equatorial (eci) position and velocity vectors.

    Parameters
    ----------
    rho : float
        satellite range from site: km
    az : float
        azimuth: 0 to 2pi rad
    el : float
        elevation: -pi/2 to pi/2 rad
    drho : float
        range rate: km/s
    daz : float
        azimuth rate: rad/s
    del_ : float
        elevation rate: rad/s
    latgd : float
        geodetic latitude of site: -pi/2 to pi2 rad
    lon : float
        longitute of site: -2pi to 2pi rad
    alt : float
        altitude of site: km
    ttt : float
        julian centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days from 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    terms : int
        # of terms for ast calculation: 0,2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta eps correction to gcrf: rad

    Returns
    -------
    reci: ndarray
        eci position vector of satellite: km
    veci: ndarray
        eci velocity vecotr of satellite: km/s
    """

    # -------------------------  implementation   -----------------
    # -----------  find sez range and velocity vectors ------------
    rhosez, drhosez = raz2rvs(rho, az, el, drho, daz, del_)
    # -----------  perform sez to ijk (ecef) transformation -------
    # tempvec = smu.rot2(rhosez, latgd - halfpi)
    # rhoecef = smu.rot3(tempvec, -lon)
    # tempvec = smu.rot2(drhosez, latgd - halfpi)
    # drhoecef = smu.rot3(tempvec, -lon)

    # alternate sez to ecef transformation; faster?
    sinlat = math.sin(latgd)
    coslat = math.cos(latgd)
    sinlon = math.sin(lon)
    coslon = math.cos(lon)
    sez2ecef = np.array([[sinlat*coslon, -sinlon, coslat*coslon],
                         [sinlat*sinlon, coslon, coslat*sinlon],
                         [-coslat, 0, sinlat]])
    rhoecef = sez2ecef @ rhosez
    drhoecef = sez2ecef @ drhosez


    # ----------  find ecef range and velocity vectors -------------
    rs, _ = obu.site(latgd, lon, alt)
    recef = rhoecef + rs
    vecef = drhoecef
    # -------- convert ecef to eci
    recef = recef
    print("in razel2rv")
    print("recef")
    print(recef)
    vecef = vecef
    print("vecef")
    print(vecef)
    a = np.zeros((3, 1))
    reci, veci, _ = ecef2eci(recef, vecef, a, ttt, jdut1,
                              lod, xp, yp, terms, ddpsi, ddeps)
    return reci, veci

# ------------------------------------------------------------------------------
#
#                           function coe2rv
#
#  this function finds the position and velocity vectors in geocentric
#    equatorial (ijk) system given the classical orbit elements.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#    vallado     - add constant file use                         29 jun 2003
#
#  inputs          description                    range / units
#    p           - semilatus rectum               km
#    ecc         - eccentricity
#    incl        - inclination                    0.0  to pi rad
#    omega       - longitude of ascending node    0.0  to 2pi rad
#    argp        - argument of perigee            0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
#    truelon     - true longitude            (ce) 0.0  to 2pi rad
#    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
#    mu          - grav parameter (optional)      km3/s2
#                  (default=earth)
#
#  outputs       :
#    r           - ijk position vector            km
#    v           - ijk velocity vector            km / s
#
#  locals        :
#    temp        - temporary real*8 value
#    rpqw        - pqw position vector            km
#    vpqw        - pqw velocity vector            km / s
#    sinnu       - sine of nu
#    cosnu       - cosine of nu
#    tempvec     - pqw velocity vector
#
#  coupling      :
#    mag         - magnitude of a vector
#    rot3        - rotation about the 3rd axis
#    rot1        - rotation about the 1st axis
#
#  references    :
#    vallado       2007, 126, alg 10, ex 2-5
#
# [r, v] = coe2rv (p, ecc, incl, omega, argp, nu, arglat, truelon, lonper)
# ------------------------------------------------------------------------------


def coe2rv(p: float, ecc: float, incl: float, omega: float, argp: float,
           nu: float, arglat: float, truelon: float, lonper: float, mu:float=mu):
    """This function finds the position and velocity vectors in geocentric
    equatorial (ijk) system given the classical orbit elements.

    Parameters
    ----------
    p : float
        semilatus rectum: km
    ecc : float
        eccentricity
    incl : float
        inclination: 0.0  to pi rad
    omega : float
        longitude of ascending node: 0.0  to 2pi rad
    argp : float
        argument of pedigree: 0.0  to 2pi rad
    nu : float
        true anomaly: 0.0  to 2pi rad
    arglat : float
        argument of latitude: 0.0  to 2pi rad
    truelon : float
        true longitude: 0.0  to 2pi rad
    lonper : float
        longitude of periapsis: 0.0  to 2pi rad
    mu: float
        graviational parameter: km3/s2, default = earth mu

    Returns
    -------
    r : ndarray
        ijk position vector: km
    v : ndarray
        ijk velocity vector: km/s
    """

    # -------------------------------------------------------------
    #       determine what type of orbit is involved and set up the
    #       set up angles for the special cases.
    # -------------------------------------------------------------
    if (ecc < small):
        # ----------------  circular equatorial  ------------------
        if (incl < small) or (abs(incl - math.pi) < small):
            argp = 0.0
            omega = 0.0
            nu = truelon
        else:
            # --------------  circular inclined  ------------------
            argp = 0.0
            nu = arglat
    else:
            # ---------------  elliptical equatorial  -----------------
        if ((incl < small) or (abs(incl - math.pi) < small)):
            argp = lonper
            omega = 0.0
    if sh.show:
        print("here, argp =%.3f, omega =%.3f, nu =%.3f" % (argp, omega, nu))
    # ----------  form pqw position and velocity vectors ----------
    cosnu = math.cos(nu)
    sinnu = math.sin(nu)
    temp = p / (1.0  + ecc * cosnu)
    rpqw = np.zeros((3))
    rpqw[0] = temp * cosnu
    rpqw[1] = temp * sinnu
    rpqw[2] = 0.0
    if (abs(p) < 0.0001):
        p = 0.0001
    vpqw = np.zeros(3)
    vpqw[0] = - sinnu * np.sqrt(mu) / np.sqrt(p)
    vpqw[1] = (ecc + cosnu) * np.sqrt(mu) / np.sqrt(p)
    vpqw[2] = 0.0

    # ----------------  perform transformation to ijk  ------------
    tempvec = smu.rot3(rpqw, -argp)
    tempvec = smu.rot1(tempvec, -incl)
    r = smu.rot3(tempvec, -omega)

    tempvec = smu.rot3(vpqw, -argp)
    tempvec = smu.rot1(tempvec, -incl)
    v = smu.rot3(tempvec, -omega)

    r = r.T #r = r' (transpose)
    v = v.T #v = v' (transpose)
    return r, v

# ------------------------------------------------------------------------------
#
#                           function rv2coe
#
#  this function finds the classical orbital elements given the geocentric
#    equatorial position and velocity vectors.
#
#  author        : david vallado                  719-573-2600   21 jun 2002
#
#  revisions
#    vallado     - fix special cases                              5 sep 2002
#    vallado     - delete extra check in inclination code        16 oct 2002
#    vallado     - add constant file use                         29 jun 2003
#    vallado     - add mu                                         2 apr 2007
#
#  inputs          description                          range / units
#    r           - ijk position vector                  km
#    v           - ijk velocity vector                  km / s
#    mu          - gravitational parameter (optional)   km3 / s2
#
#  outputs       :
#    p           - semilatus rectum                     km
#    a           - semimajor axis                       km
#    ecc         - eccentricity
#    incl        - inclination                          0.0  to pi rad
#    raan        - longitude of ascending node          0.0  to 2pi rad
#    argp        - argument of perigee                  0.0  to 2pi rad
#    nu          - true anomaly                         0.0  to 2pi rad
#    m           - mean anomaly                         0.0  to 2pi rad
#    arglat      - argument of latitude      (ci)       0.0  to 2pi rad
#    truelon     - true longitude            (ce)       0.0  to 2pi rad
#    lonper      - longitude of periapsis    (ee)       0.0  to 2pi rad
#
#  locals        :
#    hbar        - angular momentum h vector            km2 / s
#    ebar        - eccentricity     e vector
#    nbar        - line of nodes    n vector
#    c1          - v**2 - u/r
#    rdotv       - r dot v
#    hk          - hk unit vector
#    sme         - specfic mechanical energy            km2 / s2
#    i           - index
#    e           - eccentric, parabolic,
#                  hyperbolic anomaly                   rad
#    temp        - temporary variable
#    typeorbit   - type of orbit                        ee, ei, ce, ci
#
#  coupling      :
#    mag         - magnitude of a vector
#    angl        - find the angl between two vectors
#    newtonnu    - find the mean anomaly
#
#  references    :
#    vallado       2007, 121, alg 9, ex 2-5
#
# [p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r, v)
# ------------------------------------------------------------------------------


def rv2coe(r: np.ndarray, v: np.ndarray, mu=mu):
    """this function finds the classical orbital elements given the geocentric
    equatorial position and velocity vectors.

    Parameters
    ----------
    r : ndarray
        ijk position vector
    v : ndarray
        ijk velocity vector
    mu
        gravitational parameter: km3 / s2 (default set to earth)

    Returns
    -------
    p
        semilatus rectum: km
    a
        semimajor axis: km
    ecc
        eccentricity
    incl
        inclination: 0.0  to pi rad
    raan
        longitude of ascending node: 0.0  to 2pi rad
    argp
        argument of perigee: 0.0  to 2pi rad
    nu
        true anomaly: 0.0  to 2pi rad
    m
        mean anomaly: 0.0  to 2pi rad
    arglat
        argument of latitude: Circular Inclined: 0.0  to 2pi rad, else None
    truelon
        true longitude: Circular Equatorial: 0.0  to 2pi rad, else None
    lonper
        longitude of periapsis: Elliptical equatorial: 0.0  to 2pi rad, else None
    """

    m = None
    small = 1.0e-12
    # -------------------------  implementation   -----------------
    magr = smu.mag(r)
    magv = smu.mag(v)
    # ------------------  find h n and e vectors   ----------------
    hbar = np.cross(r, v)
    magh = smu.mag(hbar)
    if (magh >= small):
        nbar = np.zeros((3))
        nbar[0] = -hbar[1]
        nbar[1] = hbar[0]
        nbar[2] = 0.0
        magn = smu.mag(nbar)
        c1 = (magv * magv) - (mu / magr)
        rdotv = np.dot(r, v)
        ebar = np.zeros(3)
        for i in range(3):
            ebar[i] = (c1 * r[i] - rdotv * v[i]) / mu
        ecc = smu.mag(ebar)

        # ------------  find a e and semi-latus rectum   ----------
        sme = (magv * magv * 0.5) - (mu / magr)
        if (abs(sme) > small and ecc != 1.0):
            a = -mu / (2.0 * sme)
            p = a * (1 - ecc * ecc)
        else:
            a = infinite
            p = magh * magh / mu

        # -----------------  find inclination   -------------------

        hk = hbar[2] / magh
        incl = math.acos(hk)

        # --------  determine type of orbit for later use  --------
        # ------ elliptical, parabolic, hyperbolic inclined -------
        typeorbit = 'ei'
        if (ecc < small):
            # ----------------  circular equatorial ---------------
            if  (incl < small) or (abs(incl-math.pi)<small):
                typeorbit = 'ce'
            else:
                # --------------  circular inclined ---------------
                typeorbit = 'ci'
        else:
            # - elliptical, parabolic, hyperbolic equatorial --
            if  (incl < small) or (abs(incl-math.pi)<small):
                typeorbit = 'ee'

        if sh.show:
            print("Orbit type in rv2coe is ", typeorbit)

        # ----------  find right ascension of ascending node ------------
        if (magn > small):
            temp = nbar[0] / magn
            if (abs(temp > 1.0)):
                temp = np.sign(temp)

            raan = math.acos(temp)
            if (nbar[1] < 0.0):
                raan = twopi - raan

        else:
            raan = None


        # ---------------- find argument of perigee ---------------
        if (typeorbit == 'ei'):
            argp = smu.angl(nbar, ebar)
            if (ebar[2] < 0.0):
                argp = twopi - argp
        else:
            argp = None

        # ------------  find true anomaly at epoch    -------------
        if typeorbit.startswith('e'):
            nu = smu.angl(ebar, r)
            if (rdotv < 0.0):
                nu = twopi - nu

        # "True anomaly is not defined for circular orbits because
        # they have no periapsis. We can overcome this limitation by selecting a
        # direction in the orbit to replace periapsis as the location for the
        # initial measurement. Computer-software routines must account for this
        # special case." (pg 18)
        else:
            nu = None

        # ----  find argument of latitude - circular inclined -----
        # -- find in general cases too
        if (typeorbit == 'ci') or (typeorbit == 'ei'):
            arglat = smu.angl(nbar, r)
            if (r[2] < 0.0):
                arglat = twopi - arglat
            m = arglat
        else:
            arglat = None

        # -- find longitude of perigee - elliptical equatorial ----
        if  (ecc > small) and (typeorbit == 'ee'):
            temp = ebar[0] / ecc
            if (abs(temp) > 1.0):
                temp = np.sign(temp)

            lonper = math.acos(temp)
            if (ebar[1] < 0.0):
                lonper = twopi - lonper

            if (incl > halfpi):
                lonper = twopi - lonper

        else:
            lonper = None

        # -------- find true longitude - circular equatorial ------
        if  (magr > small) and (typeorbit == 'ce'):
            temp = r[0] / magr
            if (abs(temp) > 1.0):
                temp = np.sign(temp)
            truelon = math.acos(temp)
            if (r[1] < 0.0):
                truelon = twopi - truelon
            if (incl > halfpi):
                truelon = twopi - truelon
            m = truelon
        else:
            truelon = None

        # ------------ find mean anomaly for all orbits -----------

        if typeorbit.startswith('e'):
            _, m = smu.newtonnu(ecc, nu)

        if sh.show:
            print("ecc =", ecc)
            print("nu =", nu)

    else:
       p=None
       a=None
       ecc=None
       incl=None
       raan=None
       argp=None
       nu=None
       m=None
       arglat=None
       truelon=None
       lonper=None

    return p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper



# ----------------------------------------------------------------------------
#
#                           function ecef2eci
#
#  this function transforms a vector from the earth fixed (itrf) frame, to
#    the eci mean equator mean equinox (j2000).
#
#  author        : david vallado                  719-573-2600    4 jun 2002
#
#  revisions
#    vallado     - add terms for ast calculation                 30 sep 2002
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    eqeterms    - terms for ast calculation      0, 2
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#
#  locals        :
#    deltapsi    - nutation angle                 rad
#    trueeps     - true obliquity of the ecliptic rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       -                                rad
#    prec        - matrix for mod - eci
#    nut         - matrix for tod - mod
#    st          - matrix for pef - tod
#    stdot       - matrix for pef - tod rate
#    pm          - matrix for ecef - pef
#
#  coupling      :
#   precess      - rotation for precession
#   nutation     - rotation for nutation
#   sidereal     - rotation for sidereal time
#   polarm       - rotation for polar motion
#
#  references    :
#    vallado       2013, 223-229
#
# [reci, veci, aeci] = ecef2eci  (recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)
# ----------------------------------------------------------------------------


def ecef2eci(recef: np.ndarray, vecef: np.ndarray, aecef: np.ndarray,
             ttt: float, jdut1: float, lod: float, xp: float, yp: float,
             eqeterms: int, ddpsi: float, ddeps: float):
    """ this function transforms a vector from the earth fixed (itrf) frame, to
    the eci mean equator mean equinox (j2000).

    Parameters
    ----------
    recef : ndarray
        position vector earth fixed: km
    vecef : ndarray
        velocity vector earth fixed: km/s
    aecef : ndarray
        acceleration vector earth fixed: km/s2
    ttt : float
        julian centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    eqeterms : int
        terms for ast calculation: 0, 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta eps correction to gcrf: rad

    Returns
    -------
    reci : ndarray
        position vector eci: km
    veci : ndarray
        velocity vector eci: km/s
    aeci : ndarray
        acceleration vector eci: km/s2
    """

    # ---- find matrices
    prec, _, _, _, _ = obu.precess(ttt, '80')

    deltapsi, _, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)

    st, _ = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

    pm = smu.polarm(xp, yp, ttt, '80')

    # ---- perform transformations
    thetasa = earthrot * (1.0  - lod/86400.0)
    #        omegaearth = np.array([[0], [0], [thetasa]])
    #omegaearth = np.array([[0.0], [0.0], [thetasa]])
    omegaearth = np.array([0.0, 0.0, thetasa])
    #print("in ecef2eci:")

    tmpmat = prec @ nut @ st
    rpef = pm @ recef
    reci = tmpmat @ rpef

    vpef = pm @ vecef
    temp = np.cross(omegaearth, rpef.T) #turn from column to row vectors for cross product
    veci = tmpmat @ (vpef + temp.T)

    # veci1 = prec*nut * (stdot*recef + st*pm*vecef)  % alt approach using sidereal rate

    # two additional terms not needed if satellite is not on surface
    # of the Earth
    c1 = np.cross(omegaearth, temp)
    c2 = 2.0*np.cross(omegaearth, vpef.T)
    aeci = tmpmat @ (pm @ aecef + c1.T + c2.T)
    #print("aeci: ", aeci)

    return reci, veci, aeci



# ------------------------------------------------------------------------------
#
#                           function ecef2lle
#
#  these subroutines convert a geocentric equatorial (ijk) position vector into
#    latitude and longitude.  geodetic and geocentric latitude are found.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    r           - ijk position vector            km
#    jd          - julian date                    days from 4713 bc
#
#  outputs       :
#    latgc       - geocentric latitude            -pi to pi rad
#    latgd       - geodetic latitude              -pi to pi rad
#    lon         - longitude (west -)             -2pi to 2pi rad
#    hellp       - height above the ellipsoid     km
#
#  locals        :
#    rc          - range of site wrt earth center er
#    height      - height above earth wrt site    er
#    alpha       - angle from iaxis to point, lst rad
#    olddelta    - previous value of deltalat     rad
#    deltalat    - diff between delta and
#                  geocentric lat                 rad
#    delta       - declination angle of r in ijk  rad
#    rsqrd       - magnitude of r squared         er2
#    sintemp     - sine of temp                   rad
#    c           -
#
#  coupling      :
#    mag         - magnitude of a vector
#    gstime      - greenwich sidereal time
#    gcgd        - converts between geocentric and geodetic latitude
#
#  references    :
#    vallado       2001, 174-179, alg 12 and alg 13, ex 3-3
#
# [latgc, latgd, lon, hellp] = ecef2ll (r, jd)
# ------------------------------------------------------------------------------


def ecef2lle(r, jd):
#  from astropy.time import Time

#  mjd_ob = Time(dtg)
#  jd = mjd_ob.mjd
#  print(dtg, jd)

  # -------------------------  implementation   -----------------
  small = 0.00000001         # small value for tolerances

  # -------------------  initialize values   --------------------
  magr = smu.mag(r)
  oneminuse2 = 1.0  - eccearthsqrd

  # ---------------- find longitude value  ----------------------
  temp = np.sqrt(r[0]*r[0] + r[1]*r[1])
  if (math.fabs(temp) < small):
      rtasc = math.sign(r[2])*math.pi*0.5
  else:
      rtasc = math.atan2(r[1] , r[0])
  gst = stu.gstime(jd)
  lon = rtasc - gst

  if (math.fabs(lon) >= math.pi):
      if (lon < 0.0):
          lon = twopi + lon
      else:
          lon = lon - twopi

  # -------------- set up initial latitude value  ---------------
  decl = math.asin(r[2] / magr)
  latgc = decl
  deltalat = 100.0
  print(magr)
  rsqrd = magr*magr

  # ---- iterate to find geocentric and geodetic latitude  -----
  i = 1
  olddelta = deltalat*100 #have to initialize it to something ???
  while ((math.fabs(olddelta - deltalat) >= small) and (i < 10)):
      olddelta = deltalat
      rsite = np.sqrt(oneminuse2 / (1.0  - eccearthsqrd*(math.cos(latgc))**2))
      latgd = math.atan(math.tan(latgc) / oneminuse2)
      temp = latgd-latgc
      sintemp = math.sin(temp)
      hellp = (np.sqrt(rsqrd - rsite*rsite*sintemp*sintemp)
                  - rsite*math.cos(temp))
      deltalat = math.asin(hellp*sintemp / magr)
      latgc = decl - deltalat
      i = i + 1

  if (i >= 10):
      print('ijktolatlon did not converge\n ')

  latgc = latgc/math.pi*90.0
  latgd = latgd/math.pi*90.0
  lon = lon/math.tau*180.0


  return latgc, latgd, lon, hellp



# ------------------------------------------------------------------------------
#
#                           function ecef2ll
#
#  these subroutines convert a geocentric equatorial position vector into
#    latitude and longitude.  geodetic and geocentric latitude are found. the
#    inputs must be ecef.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#    vallado     - fix jdut1 var name, add clarifying comments   26 aug 2002
#    vallado     - fix documentation for ecef                    19 jan 2005
#
#  inputs          description                    range / units
#    r           - ecef position vector           km
#
#  outputs       :
#    latgc       - geocentric latitude            -pi to pi rad
#    latgd       - geodetic latitude              -pi to pi rad
#    lon         - longitude (west -)             -2pi to 2pi rad
#    hellp       - height above the ellipsoid     km
#
#  locals        :
#    temp        - diff between geocentric/
#                  geodetic lat                   rad
#    sintemp     - sine of temp                   rad
#    olddelta    - previous value of deltalat     rad
#    rtasc       - right ascension                rad
#    decl        - declination                    rad
#    i           - index
#
#  coupling      :
#    mag         - magnitude of a vector
#    gcgd        - converts between geocentric and geodetic latitude
#
#  references    :
#    vallado       2001, 174-179, alg 12 and alg 13, ex 3-3
#
# [latgc, latgd, lon, hellp] = ecef2ll (r)
# ------------------------------------------------------------------------------

def ecef2ll (r):

        # -------------------------  implementation   -----------------
        magr = np.linalg.norm(r)

        # ----------------- find longitude value  ---------------------
        temp = np.sqrt(r[0]*r[0] + r[1]*r[1])
        if (abs(temp) < small):
            rtasc = np.sign(r[2])*math.pi*0.5
        else:
            rtasc = math.atan2(r[1], r[0])
        lon = rtasc
        if (abs(lon) >= math.pi):   # mod it ?
            if (lon < 0.0):
                lon = twopi + lon
            else:
                lon = lon - twopi
        decl = math.asin(r[2] / magr)
        latgd = decl
   #     print('rd %11.7f rtasc %11.7f lon %11.7f latgd %11.7f \n', temp, rtasc * rad2deg, lon * rad2deg, latgd * rad2deg)

        # ------------- iterate to find geodetic latitude -------------
        i = 1
        olddelta = latgd + 10.0

        while ((abs(olddelta-latgd)>= small) and (i<10)):
            olddelta = latgd
            sintemp = math.sin(latgd)
            c = re  / (np.sqrt(1.0 -eccearthsqrd*sintemp*sintemp))
            latgd = math.atan((r[2]+c*eccearthsqrd*sintemp)/temp)
           # print('%3i  c %11.7f gd %11.7f  \n', i, c, latgd * rad2deg)
            i = i + 1

        # Calculate height
        if (math.pi*0.5 - abs(latgd)) < (math.pi/180.0):  # 1 deg
            hellp = (temp/math.cos(latgd)) - c
        elif latgd ==0: #we added this to deal with z ==0
            hellp = magr-re
        else:
            s = c * (1.0 - eccearthsqrd)
            hellp = r[2]/math.sin(latgd) - s

        latgc = math.asin(r[2]/magr)   # all locations
        #latgc = gd2gc(latgd)  # surface of the Earth locations
        return latgc, latgd, lon, hellp



# ------------------------------------------------------------------------------
#
#                           function ecef2llb
#
#  these subroutines convert a geocentric equatorial (ijk) position vector into
#    latitude and longitude.  geodetic and geocentric latitude are found.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    r           - ijk position vector            km
#
#  outputs       :
#    latgc       - geocentric latitude            -pi to pi rad
#    latgd       - geodetic latitude              -pi to pi rad
#    lon         - longitude (west -)             -2pi to 2pi rad
#    hellp       - height above the ellipsoid     km
#
#  locals        :
#    rc          - range of site wrt earth center er
#    height      - height above earth wrt site    er
#    alpha       - angle from iaxis to point, lst rad
#    olddelta    - previous value of deltalat     rad
#    deltalat    - diff between delta and
#                  geocentric lat                 rad
#    delta       - declination angle of r in ijk  rad
#    rsqrd       - magnitude of r squared         er2
#    sintemp     - sine of temp                   rad
#    c           -
#
#  coupling      :
#    mag         - magnitude of a vector
#    gcgd        - converts between geocentric and geodetic latitude
#
#  references    :
#    vallado       2001, 174-179, alg 12 and alg 13, ex 3-3
#
# [latgc, latgd, lon, hellp] = ecef2llb (r)
# ------------------------------------------------------------------------------

#this function breaks at math.atan when given [100000, 200000, 0.0]
def ecef2llb (r):

    magr = np.linalg.norm(r)
    # -------------------------  implementation   -------------------------
    # ---------------- find longitude value  ----------------------
    temp = np.sqrt(r[0]*r[0] + r[1]*r[1])
    if (abs(temp) < small):
        rtasc = np.sign(r[2])*math.pi*0.5
    else:
        rtasc = math.atan2(r[1] , r[0])
    lon = rtasc
    if (abs(lon) >= math.pi):
        if (lon < 0.0):
            lon = twopi + lon
        else:
            lon = lon - twopi

    a = 6378.1363
    # b = np.sign(r[2]) * 6356.75160056
    b = np.sign(r[2]) * a * np.sqrt(1.0 - 0.006694385)  # find semiminor axis of the earth

    # -------------- set up initial latitude value  ---------------
    atemp = 1.0 /(a*temp)
    e = (b*r[2]-a*a+b*b)*atemp
    f = (b*r[2]+a*a-b*b)*atemp
    third = 1.0 /3.0
    p = 4.0 *third*(e*f + 1.0)
    q = 2.0 *(e*e - f*f)
    d = p*p*p + q*q

    if (d > 0.0):
        nu = (np.sqrt(d)-q)**third - (np.sqrt(d)+q)**third
    else:
        sqrtp = np.sqrt(-p)
        nu = 2.0 *sqrtp*math.cos(third*np.arccos(q/(p*sqrtp)))
    g = 0.5 *(np.sqrt(e*e + nu) + e)
    t = np.sqrt(g*g + (f-nu*g)/(2.0 *g-e)) - g

    latgd = math.atan(a*(1.0 -t*t)/(2.0 *b*t))
    hellp = (temp-a*t)*math.cos(latgd) + (r[2]-b)*math.sin(latgd)

    latgc = math.asin(r[2]/magr)   # all locations
    #latgc = gd2gc(latgd)  # surface of the Earth locations
    return latgc, latgd, lon, hellp



# ----------------------------------------------------------------------------
#
#                           function ecef2mod
#
#  this function transforms a vector from the earth fixed (itrf) frame, to
#    the mean of date (mod) frame.
#
#  author        : david vallado                  719-573-2600    4 jun 2002
#
#  revisions
#
#  inputs          description                    range / units
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    eqeterms    - terms for ast calculation      0, 2
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    rmod        - position vector mod            km
#    vmod        - velocity vector mod            km/s
#    amod        - acceleration vector mod        km/s2
#
#  locals        :
#    deltapsi    - nutation angle                 rad
#    trueeps     - true obliquity of the ecliptic rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       -                                rad
#    nut         - matrix for tod - mod
#    st          - matrix for pef - tod
#    stdot       - matrix for pef - tod rate
#    pm          - matrix for ecef - pef
#
#  coupling      :
#   nutation     - rotation for nutation
#   sidereal     - rotation for sidereal time
#   polarm       - rotation for polar motion
#
#  references    :
#    vallado       2004, 219-228
#
# [rmod, vmod, amod] = ecef2mod (recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)
# ----------------------------------------------------------------------------


def ecef2mod(recef: np.ndarray, vecef: np.ndarray, aecef: np.ndarray,
             ttt: float, jdut1: float, lod: float, xp: float, yp: float,
             eqeterms: int, ddpsi: float, ddeps: float):
    """this function transforms a vector from the earth fixed (itrf) frame, to
    the mean of date (mod) frame.

    Parameters
    ----------
    recef : ndarray
        position vector earth fixed: km
    vecef : ndarray
        velocity vector earth fixed: km/s
    aecef : ndarray
        acceleration vector earth fixed: km/s2
    ttt : float
        julian centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    eqeterms : int
        terms for ast calculation: 0, 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta psi correction to gcrf: rad

    Returns
    -------
    rmod : ndarray
        position vector mod: km
    vmod : ndarray
        velocity vector mod: km/s
    amod : ndarray
        acceleration vector mod: km/s2
    """

    # ---- find matrices
    deltapsi, _, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)

    st, _ = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

    pm = smu.polarm(xp, yp, ttt, '80')

    # ---- perform transformations
    thetasa = earthrot * (1.0  - lod/86400.0)

    omegaearth = np.array([0.0, 0.0, thetasa])

    # trueeps-meaneps
    # deltapsi
    # nut
    rpef = pm@recef
    rmod = nut@st@rpef

    vpef = pm@vecef
    vmod = nut@st@(vpef + np.cross(omegaearth, rpef.T).T)

    temp = np.cross(omegaearth, rpef.T)
    amod = nut @ st @ (pm @ aecef + np.cross(omegaearth, temp).T \
            + 2.0*np.cross(omegaearth, vpef.T).T)

    return rmod, vmod, amod

def mod2ecef(rmod: np.ndarray, vmod: np.ndarray, amod: np.ndarray,
             ttt: float, jdut1: float, lod: float, xp: float, yp: float,
             eqeterms: int, ddpsi: float, ddeps: float):
    """this function transforms a vector from the earth fixed (itrf) frame, to
    the mean of date (mod) frame.

    Parameters
    ----------
    rmod : ndarray
        position vector mod: km
    vmod : ndarray
        velocity vector mod: km/s
    amod : ndarray
        acceleration vector mod: km/s2
    ttt : float
        julian centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    eqeterms : int
        terms for ast calculation: 0, 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta psi correction to gcrf: rad

    Returns
    -------
    recef: ndarray
        position vector ecef: km
    vecef: ndarray
        velocity vector ecef: km/s
    aecef: ndarray
        acceleration vector ecef: km/s2
    """

    # ---- find matrices
    deltapsi, _, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)

    st, _ = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

    pm = smu.polarm(xp, yp, ttt, '80')

    # ---- perform transformations
    thetasa = earthrot * (1.0  - lod/86400.0)

    omegaearth = np.array([0.0, 0.0, thetasa])

    rpef = st.T @ nut.T @ rmod
    recef = pm.T @ rpef

    vpef = st.T @ nut.T @ vmod
    vecef = pm.T @ (vpef - np.cross(omegaearth, rpef.T).T)

    temp = np.cross(omegaearth, rpef.T)
    aecef = pm.T @ (st.T @ nut.T @ amod - np.cross(omegaearth, temp).T
                    - 2.0*np.cross(omegaearth, vpef.T).T)
    return recef, vecef, aecef


# ----------------------------------------------------------------------------
#
#                           function ecef2tod
#
#  this function transforms a vector from the earth fixed (itrf) frame, to
#    the true of date (tod).
#
#  author        : david vallado                  719-573-2600    4 jun 2002
#
#  revisions
#
#  inputs          description                    range / units
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    eqeterms    - terms for ast calculation      0, 2
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    rtod        - position vector tod            km
#    vtod        - velocity vector tod            km/s
#    atod        - acceleration vector tod        km/s2
#
#  locals        :
#    deltapsi    - nutation angle                 rad
#    trueeps     - true obliquity of the ecliptic rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       -                                rad
#    prec        - matrix for mod - eci
#    nut         - matrix for tod - mod
#    st          - matrix for pef - tod
#    stdot       - matrix for pef - tod rate
#    pm          - matrix for ecef - pef
#
#  coupling      :
#   sidereal     - rotation for sidereal time
#   polarm       - rotation for polar motion
#
#  references    :
#    vallado       2013, 223-231
#
#  [rtod, vtod, atod] = ecef2tod(recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)
# ----------------------------------------------------------------------------


def ecef2tod(recef: np.ndarray, vecef: np.ndarray, aecef: np.ndarray,
             ttt: float, jdut1: float, lod: float, xp: float, yp: float,
             eqeterms: int, ddpsi: float, ddeps: float):
    """this function transforms a vector from the earth fixed (itrf) frame, to
    the true of date (tod) frame.

    Parameters
    ----------
    recef : ndarray
        position vector earth fixed: km
    vecef : ndarray
        velocity vector earth fixed: km/s
    aecef : ndarray
        acceleration vector earth fixed: km/s2
    ttt : float
        julian centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    eqeterms : int
        terms for ast calculation: 0, 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta psi correction to gcrf: rad

    Returns
    -------
    rtod : ndarray
        position vector tod: km
    vtod : ndarray
        velocity vector tod: km/s
    atod : ndarray
        acceleration vector tod: km/s2
    """
    # ---- find matrices - note nut is only needed for st argument inputs
    deltapsi, _, meaneps, omega, _ = obu.nutation(ttt, ddpsi, ddeps)

    st, _ = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

    pm = smu.polarm(xp, yp, ttt, '80')

    # ---- perform transformations
    thetasa = earthrot * (1.0  - lod/86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])

    rpef = pm @ recef
    rtod = st @ rpef

    vpef = pm @ vecef
    vtod = st @ (vpef + np.cross(omegaearth, rpef.T).T)

    temp = np.cross(omegaearth, rpef.T)
    atod = st @ (pm @ aecef + np.cross(omegaearth, temp).T \
            + 2.0*np.cross(omegaearth, vpef.T).T)

    return rtod, vtod, atod

def tod2ecef(rtod: np.ndarray, vtod: np.ndarray, atod: np.ndarray,
             ttt: float, jdut1: float, lod: float, xp: float, yp: float,
             eqeterms: int, ddpsi: float, ddeps: float):
    """this function transforms a vector from the true of date (tod) frame, to
    the earth fixed (itrf) frame.

    Parameters
    ----------
    rtod : ndarray
        position vector tod: km
    vtod : ndarray
        velocity vector tod: km/s
    atod : ndarray
        acceleration vector tod: km/s2
    ttt : float
        julian centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    eqeterms : int
        terms for ast calculation: 0, 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta psi correction to gcrf: rad

    Returns
    -------
    recef: ndarray
        position vector ecef: km
    vecef: ndarray
        velocity vector ecef: km/s
    aecef: ndarray
        acceleration vector ecef: km/s2
    """

    # ---- find matrices - note nut is only needed for st argument inputs
    deltapsi, _, meaneps, omega, _ = obu.nutation(ttt, ddpsi, ddeps)

    st, _ = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

    pm = smu.polarm(xp, yp, ttt, '80')

    # ---- perform transformations
    thetasa = earthrot * (1.0  - lod/86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])

    rpef = st.T @ rtod
    recef = pm.T @ rpef

    vpef = st.T @ vtod
    vecef = pm.T @ (vpef - np.cross(omegaearth, rpef.T).T)

    temp = np.cross(omegaearth, rpef.T)
    aecef = pm.T @ (st.T @ atod - np.cross(omegaearth, temp).T
                    - 2.0 * np.cross(omegaearth, vpef.T).T)

    return recef, vecef, aecef

# ----------------------------------------------------------------------------
#
#                           function ecef2teme
#
#  this function trsnforms a vector from the earth fixed (ITRF) frame to the
#    true equator mean equniox frame (teme). the results take into account
#    the effects of sidereal time, and polar motion.
#
#  author        : david vallado                  719-573-2600   30 oct 2017
#
#  revisions
#
#  inputs          description                    range / units
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    eqeterms    - use extra two terms (kinematic) after 1997  0, 2
#
#  outputs       :
#    rteme       - position vector teme           km
#    vteme       - velocity vector teme           km/s
#    ateme       - acceleration vector teme       km/s2
#
#  locals        :
#    st          - matrix for pef - tod
#    pm          - matrix for ecef - pef
#
#  coupling      :
#   gstime       - greenwich mean sidereal time   rad
#   polarm       - rotation for polar motion      pef - ecef
#
#  references    :
#    vallado       2013, 231-233
#
# [rteme, vteme, ateme] = ecef2teme(recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms)
# ----------------------------------------------------------------------------


def ecef2teme(recef: np.ndarray, vecef: np.ndarray, aecef: np.ndarray,
              ttt: float, jdut1: float, lod: float, xp: float, yp: float,
              eqeterms: int):
    """this function transforms a vector from the earth fixed (itrf) frame, to
    the mean of date (mod) frame.

    Parameters
    ----------
    recef : ndarray
        position vector earth fixed: km
    vecef : ndarray
        velocity vector earth fixed: km/s
    aecef : ndarray
        acceleration vector earth fixed: km/s2
    ttt : float
        julian centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    eqeterms : int
        terms for ast calculation: 0, 2

    Returns
    -------
    rteme : ndarray
        position vector teme: km
    vteme : ndarray
        velocity vector teme: km/s
    ateme : ndarray
        acceleration vector teme: km/s2
    """

    # ------------------------ find gmst --------------------------
    gmst = stu.gstime(jdut1)

    # find omeage from nutation theory
    omega = 125.04452222  + (-6962890.5390 *ttt + \
            7.455 *ttt*ttt + 0.008 *ttt*ttt*ttt) / 3600.0
    omega = np.fmod(omega, 360.0) * deg2rad

    # teme does not include the geometric terms here
    # after 1997, kinematic terms apply
    if (jdut1 > 2450449.5) and (eqeterms > 0):
        gmstg = gmst + 0.00264 * arcsec2rad * math.sin(omega) \
                + 0.000063 * arcsec2rad * math.sin(2.0 *omega)
    else:
        gmstg = gmst

    gmstg = np.fmod (gmstg, 2.0*math.pi)

    st = np.zeros((3, 3))
    st[0, 0] = math.cos(gmstg)
    st[0, 1] = -math.sin(gmstg)
    st[0, 2] = 0.0
    st[1, 0] = math.sin(gmstg)
    st[1, 1] = math.cos(gmstg)
    st[1, 2] = 0.0
    st[2, 0] = 0.0
    st[2, 1] = 0.0
    st[2, 2] = 1.0

    pm = smu.polarm(xp, yp, ttt, '80')

    rpef = pm@recef
    rteme = st@rpef

    thetasa = earthrot * (1.0  - lod/86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])

    vpef = pm@vecef
    temp = np.cross(omegaearth, rpef.T)
    vteme = st@(vpef + temp.T)

    tt1 = np.cross(omegaearth, temp).T
    print(tt1)
    tt2 = np.cross(omegaearth, vpef.T).T
    print(tt2)

    ateme = st@(pm@aecef + tt1 + 2.0*tt2)

    return rteme, vteme, ateme


# ----------------------------------------------------------------------------
#
#                           function ecef2pef
#
#  this function transforms a vector from the earth fixed itrf frame
#    (itrf), to the pseudo earth fixed frame (pef).
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#
#  inputs          description                    range / units
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#    opt         - arg to pass through to polarm  -added jmb
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    ttt         - julian centuries of tt         centuries
#
#  outputs       :
#    rpef        - position pseudo earth fixed    km
#    vpef        - velocity pseudo earth fixed    km/s
#    apef        - acceleration pseudo earth fixedkm/s2
#
#  locals        :
#
#  coupling      :
#   precess      - rotation for precession        mod - eci
#
#  references    :
#    vallado       2001, 219, eq 3-65 to 3-66
#
# [rpef, vpef, apef] = ecef2pef  (recef, vecef, aecef, opt, xp, yp, ttt)
# ----------------------------------------------------------------------------


def ecef2pef(recef: np.ndarray, vecef: np.ndarray, aecef: np.ndarray,
             opt: str, xp: float, yp: float, ttt: float):
    """this function transforms a vector from the earth fixed itrf frame
    (itrf), to the pseudo earth fixed frame (pef).

    Parameters
    ----------
    recef : ndarray
        position vector earth fixed: km
    vecef : ndarray
        velocity vector earth fixed: km/s
    aecef : ndarray
        acceleration vector earth fixed: km/s2
    opt : str
        polarm method option: '01', '02', '80'
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    ttt : float
        julian centuries of tt: centuries

    Returns
    -------
    rpef : ndarray
        position vector pef: km
    vpef : ndarray
        velocity vector pef: km/s
    apef : ndarray
        acceleration vector pef: km/s2
    """

    pm = smu.polarm(xp, yp, ttt, opt)

    rpef = pm@recef

    vpef = pm@vecef

    apef = pm@aecef

    return rpef, vpef, apef

def pef2ecef(rpef: np.ndarray, vpef: np.ndarray, apef: np.ndarray,
             opt: str, xp: float, yp: float, ttt: float):
    """this function transforms a vector from the pseudo earth fixed frame
    (pef) to the earth fixed itrf frame (itrf).

    Parameters
    ----------
    rpef : ndarray
        position vector pef: km
    vpef : ndarray
        velocity vector pef: km/s
    apef : ndarray
        acceleration vector pef: km/s2
    opt : str
        polarm method option: '01', '02', '80'
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    ttt : float
        julian centuries of tt: centuries

    Returns
    -------
    recef: ndarray
        position vector ecef: km
    vecef: ndarray
        velocity vector ecef: km/s
    aecef: ndarray
        acceleration vector ecef: km/s2
    """

    pm = smu.polarm(xp, yp, ttt, opt)

    recef = pm.T@rpef

    vecef = pm.T@vpef

    aecef = pm.T@apef

    return recef, vecef, aecef


# ----------------------------------------------------------------------------
#
#                           function eci2ecef
#
#  this function trsnforms a vector from the mean equator mean equniox frame
#    (j2000), to an earth fixed (ITRF) frame.  the results take into account
#    the effects of precession, nutation, sidereal time, and polar motion.
#
#  author        : david vallado                  719-573-2600   27 jun 2002
#
#  revisions
#    vallado     - add terms for ast calculation                 30 sep 2002
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    eqeterms    - terms for ast calculation      0, 2
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#
#  locals        :
#    deltapsi    - nutation angle                 rad
#    trueeps     - true obliquity of the ecliptic rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       -                                rad
#    prec        - matrix for mod - eci
#    nut         - matrix for tod - mod
#    st          - matrix for pef - tod
#    stdot       - matrix for pef - tod rate
#    pm          - matrix for ecef - pef
#
#  coupling      :
#   precess      - rotation for precession        eci - mod
#   nutation     - rotation for nutation          mod - tod
#   sidereal     - rotation for sidereal time     tod - pef
#   polarm       - rotation for polar motion      pef - ecef
#
#  references    :
#    vallado       2013, 223-229
#
# [recef, vecef, aecef] = eci2ecef  (reci, veci, aeci, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps)
# ----------------------------------------------------------------------------


def eci2ecef(reci: np.ndarray, veci: np.ndarray, aeci: np.ndarray, ttt: float,
             jdut1: float, lod: float, xp: float, yp: float, eqeterms: int,
             ddpsi: float, ddeps: float):
    """this function trsnforms a vector from the mean equator mean equniox frame
    (j2000), to an earth fixed (ITRF) frame.  the results take into account
    the effects of precession, nutation, sidereal time, and polar motion.

    Parameters
    ----------
    reci : ndarray
        eci position vector: km
    veci : ndarray
        eci velocity vector: km/s
    aeci : ndarray
        eci acceleration vector: km/s2
    latgd : float
        geodetic latitude of site: rad
    lon : float
        longituge of site: rad
    alt : float
        altitude of site: km
    ttt : float
        julian centuries: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess lenght of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    terms : int
        number of terms for ast calculation: 0, 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta eps correction to gcrf: rad

    Returns
    -------
    recef: ndarray
        position vector earth fixed: km
    vecef: ndarray
        velocity vector earth fixed: km/s
    aecef: ndarray
        acceleration vector earth fixed: km/s2
    """
    prec, _, _, _, _ = obu.precess(ttt, '80')

    deltapsi, _, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)

    st, _ = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

    pm = smu.polarm(xp, yp, ttt, '80')

    thetasa = earthrot * (1.0  - lod/86400.0)
    # omegaearth = np.array([[0], [0], [thetasa]])

    omegaearth = np.array([0.0, 0.0, thetasa])

    # [r]T[n]T[p]T
    tmpmat = st.T@nut.T@prec.T

    rpef = tmpmat@reci
    recef = pm.T@rpef

    temp = np.cross(omegaearth, rpef.T) #turn from column to row vectors for cross product
    vpef = tmpmat @ veci - temp.T
    vecef = pm.T @ vpef

    # two additional terms not needed if satellite is not on surface
    # of the Earth
    c1 = np.cross(omegaearth, temp)
    c2 = 2.0*np.cross(omegaearth, vpef.T)
    aecef = pm.T @ (tmpmat @ aeci - c1.T - c2.T)
    return recef, vecef, aecef


# ----------------------------------------------------------------------------
#
#                           function eci2mod
#
#  this function transfroms a vector from the mean equator, mean equinox frame
#    (j2000), to the mean equator mean equinox of date (mod).
#
#  author        : david vallado                  719-573-2600   27 jun 2002
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    reci        - position vector eci          km
#    veci        - velocity vector eci          km/s
#    aeci        - acceleration vector eci      km/s2
#    ttt         - julian centuries of tt         centuries
#
#  outputs       :
#    rmod        - position vector of date
#                    mean equator, mean equinox   km
#    vmod        - velocity vector of date
#                    mean equator, mean equinox   km/s
#    amod        - acceleration vector of date
#                    mean equator, mean equinox   km/s2
#
#  locals        :
#    none.
#
#  coupling      :
#   precess      - rotation for precession        eci - mod
#
#  references    :
#    vallado       2001, 214-215, eq 3-57
#
# [prec] = precession  (ttt)
# ----------------------------------------------------------------------------


def eci2mod(reci: np.ndarray, veci: np.ndarray, aeci: np.ndarray, ttt: float):
    """this function transfroms a vector from the mean equator, mean equinox
    frame (j2000), to the mean equator mean equinox of date (mod).

    Parameters
    ----------
    reci : ndarray
        position vector eci: km
    veci : ndarray
        velocity vector eci: km
    aeci : ndarray
        acceleration vector eci: km
    ttt : float
        julian centuries of tt

    Returns
    -------
    rmod : ndarray
        position vector mod: km
    vmod : ndarray
        velocity vector mod: km/s
    amod : ndarray
        acceleration vector mod: km/s2
    """

    prec, _, _, _, _ = obu.precess(ttt, '80')

    rmod = prec.T@reci

    vmod = prec.T@veci

    amod = prec.T@aeci

    return rmod, vmod, amod



# ----------------------------------------------------------------------------
#
#                           function eci2tod
#
#  this function transforms a vector from the mean equator mean equinox frame
#    (j2000) to the true equator true equinox of date (tod).
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#    opt         - calculation option             '80, '6a', '6b', '6c'
#    ttt         - julian centuries of tt         centuries
#    ddpsi       - correction for iau2000         rad
#    ddeps       - correction for iau2000         rad
#    ddx         - eop correction for x           rad
#    ddy         - eop correction for y           rad
#
#  outputs       :
#    rtod        - position vector of date
#                    true equator, true equinox   km
#    vtod        - velocity vector of date
#                    true equator, true equinox   km/s
#    atod        - acceleration vector of date
#                    true equator, true equinox   km/s2
#
#  locals        :
#    prec        - matrix for eci - mod
#    deltapsi    - nutation angle                 rad
#    trueeps     - true obliquity of the ecliptic rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       -                                rad
#    nut         - matrix for mod - tod
#
#  coupling      :
#   precess      - rotation for precession        mod - eci
#   nutation     - rotation for nutation          tod - mod
#
#  references    :
#    vallado       2001, 216-219, eq 3-654
#
# [rtod, vtod, atod] = eci2tod  (reci, veci, aeci, opt, ttt, ddpsi, ddeps, ddx, ddy)
# ----------------------------------------------------------------------------


def eci2tod(reci: np.ndarray, veci: np.ndarray, aeci: np.ndarray, ttt: float,
            ddpsi: float, ddeps: float):
    """this function transforms a vector from the mean equator mean equinox
    frame (j2000) to the true equator true equinox of date (tod).

    Parameters
    ----------
    reci : ndarray
        position vector eci: km
    veci : ndarray
        velocity vector eci: km/s
    aeci : ndarray
        acceleration vector eci: km/s2
    opt : str
        calculation option: '80, '6a', '6b', '6c'
    ttt : float
        julian centuries of tt: centuries
    ddpsi : float
        correction for iau2000: rad
    ddeps : float
        correction for iau2000: rad

    Returns
    -------
    rtod : ndarray
        position vector tod: km
    vtod : ndarray
        velocity vector tod: km/s
    atod : ndarray
        acceleration vector tod: km/s2
    """


    prec, _, _, _, _ = obu.precess(ttt, '80')
    deltapsi, trueeps, meaneps, _, nut = obu.nutation(ttt, ddpsi, ddeps)
    if sh.show:
        print('dpsi %11.7f trueeps %11.7f mean eps %11.7f deltaeps %11.7f \n'
              %(deltapsi * rad2arcsec, trueeps * rad2arcsec, meaneps * rad2arcsec,
                (trueeps-meaneps) * rad2arcsec))
        print('nut iau 76 \n')
        print('%20.14f %20.14f %20.14f \n'% (nut[0], nut[1], nut[2]))

    rtod = nut.T @ prec.T @ reci

    vtod = nut.T @ prec.T @ veci

    atod = nut.T @ prec.T @ aeci

    if sh.show:
        print('pn iau 76 \n')
        tmp = nut @ prec
        print('%20.14f %20.14f %20.14f \n' % (tmp[0], tmp[1], tmp[2]))

    return rtod, vtod, atod


# ----------------------------------------------------------------------------
#
#                           function eci2teme
#
#  this function transforms a vector from the mean equator mean equinox (j2000)
#    system (eci), to the true equator mean equinox (teme) system.
#
#  author        : david vallado                  719-573-2600   30 oct 2017
#
#  inputs          description                    range / units
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#    ttt         - julian centuries of tt         centuries
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    rteme       - position vector of date
#                    true equator, mean equinox   km
#    vteme       - velocity vector of date
#                    true equator, mean equinox   km/s
#    ateme       - acceleration vector of date
#                    true equator, mean equinox   km/s2
#
#  locals        :
#    prec        - matrix for eci - mod
#    nutteme     - matrix for mod - teme - an approximation for nutation
#    eqe         - rotation for equation of equinoxes (geometric terms only)
#    tm          - combined matrix for teme2eci
#
#  coupling      :
#   precess      - rotation for precession        eci - mod
#   nutation     - rotation for nutation          eci - tod
#
#  references    :
#    vallado       2013, 231-233
#
# [reci, veci, aeci] = teme2eci  (rteme, vteme, ateme, ttt, ddpsi, ddeps)
# ----------------------------------------------------------------------------


def eci2teme(reci: np.ndarray, veci: np.ndarray, aeci: np.ndarray,
             ttt: float, ddpsi: float, ddeps: float):
    """this function transforms a vector from the mean equator mean equinox frame
    (j2000) to the to the true equator mean equinox (teme) system.

    Parameters
    ----------
    reci : ndarray
        position vector eci: km
    veci : ndarray
        velocity vector eci: km/s
    aeci : ndarray
        acceleration vector eci: km/s2
    ttt : float
        julian centuries of tt: centuries
    ddpsi : float
        correction for iau2000: rad
    ddeps : float
        correction for iau2000: rad

    Returns
    -------
    rteme : ndarray
        position vector teme: km
    vteme : ndarray
        velocity vector teme: km/s
    ateme : ndarray
        acceleration vector teme: km/s2
    """

    prec, _, _, _, _ = obu.precess (ttt, '80')

    deltapsi, _, meaneps, _, nut = obu.nutation  (ttt, ddpsi, ddeps)

    # ------------------------ find eqeg ----------------------
    # rotate teme through just geometric terms
    eqeg = deltapsi* math.cos(meaneps)

    eqeg = np.fmod(eqeg, 2.0*math.pi)

    eqe = np.zeros((3, 3))
    eqe[0, 0] = math.cos(eqeg)
    eqe[0, 1] = math.sin(eqeg)
    eqe[0, 2] = 0.0
    eqe[1, 0] = -math.sin(eqeg)
    eqe[1, 1] = math.cos(eqeg)
    eqe[1, 2] = 0.0
    eqe[2, 0] = 0.0
    eqe[2, 1] = 0.0
    eqe[2, 2] = 1.0

    tm = eqe @ nut.T @ prec.T

    rteme = tm @ reci
    vteme = tm @ veci
    ateme = tm @ aeci

    return rteme, vteme, ateme

# ------------------------------------------------------------------------------
#
#                           function gd2gc
#
#  this function converts from geodetic to geocentric latitude for positions
#    on the surface of the earth.  notice that (1-f) squared = 1-esqrd.
#
#  author        : david vallado                  719-573-2600   30 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    latgd       - geodetic latitude              -pi to pi rad
#
#  outputs       :
#    latgc       - geocentric latitude            -pi to pi rad
#
#  locals        :
#    none.
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2001, 146, eq 3-11
#
# [latgc] = gd2gc (latgd)
# ------------------------------------------------------------------------------

def gd2gc (latgd: float):
    """this function converts from geodetic to geocentric latitude for positions
    on the surface of the earth.  notice that (1-f) squared = 1-esqrd.

    Parameters
    ----------
    latgd : float
        geodetic latitude: -pi to pi rad

    Returns
    -------
    latgc: float
        geocentric latitude: -pi to pi rad
    """
    latgc = math.atan((1.0  - eccearthsqrd) * math.tan(latgd))
    return latgc


# ------------------------------------------------------------------------------
#
#                           function gc2gd
#
#  this function converts from geodetic to geocentric latitude for positions
#    on the surface of the earth.  notice that (1-f) squared = 1-esqrd.
#
#  author        : david vallado                  719-573-2600   21 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    latgd       - geodetic latitude              -pi to pi rad
#
#  outputs       :
#    latgc       - geocentric latitude            -pi to pi rad
#
#  locals        :
#    none.
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2001, 146, eq 3-11
#
# [latgd] = gc2gd (latgc)
# ------------------------------------------------------------------------------

def gc2gd(latgc: float):
    """this function converts from geodetic to geocentric latitude for positions
    on the surface of the earth.  notice that (1-f) squared = 1-esqrd.

    Parameters
    ----------
    latgc : float
        geocentric latitude: -pi to pi rad

    Returns
    -------
    latgd: float
        geodetic latitude: -pi to pi rad
    """
    latgd = math.atan(math.tan(latgc)/(1.0  - eccearthsqrd))

    return latgd


# -----------------------------------------------------------------------------
#
#                           function dms2rad
#
#  this function converts degrees, minutes and seconds into radians.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    deg         - degrees                        0 .. 360
#    min         - minutes                        0 .. 59
#    sec         - seconds                        0.0 .. 59.99
#
#  outputs       :
#    dms         - result                         rad
#
#  locals        :
#    temp        - temporary variable
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 203, alg 17 alg 18, ex 3-8
#
# [dms] = dms2rad(deg, min, sec)
# -----------------------------------------------------------------------------

def dms2rad(deg: float, min: float, sec: float):
    """this function converts degrees, minutes and seconds into radians.

    Parameters
    ----------
    deg : float
        degrees: 0 - 360
    min : float
        minutes: 0 - 59
    sec : float
        seconds: 0 - 59.99

    Returns
    -------
    dms: float
        result: rad
    """
    dms = (deg + min/60.0 + sec/3600.0) * deg2rad
    return dms


# -----------------------------------------------------------------------------
#
#                           function rad2dms
#
#  this function converts radians to degrees, minutes and seconds.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    dms         - result                         rad
#
#  outputs       :
#    deg         - degrees                        0 .. 360
#    min         - minutes                        0 .. 59
#    sec         - seconds                        0.0 .. 59.99
#
#  locals        :
#    temp        - temporary variable
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2001, 199, alg 17 alg 18, ex 3-8
#
# [deg, min, sec] = rad2dms(dms)
# -----------------------------------------------------------------------------

def rad2dms(dms: float):
    """this function converts radians to degrees, minutes and seconds.

    Parameters
    ----------
    dms : float
        radians

    Returns
    -------
    deg: float
        degrees: 0-360
    min: float
        minutes: 0-59
    sec: float
        seconds: 0-59.99
    """
    temp = dms * rad2deg
    deg = np.fix(temp)
    min = np.fix((temp - deg)*60.0)
    sec = (temp - deg - min/60.0) * 3600.0
    return deg, min, sec

# ------------------------------------------------------------------------------
#
#                           function radec2rv
#
#  this function converts the right ascension and declination values with
#    position and velocity vectors of a satellite. uses velocity vector to
#    find the solution of singular cases.
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    rr          - radius of the satellite        er
#    rtasc       - right ascension                rad
#    decl        - declination                    rad
#    drr         - radius of the satellite rate   er/tu
#    drtasc      - right ascension rate           rad/tu
#    ddecl       - declination rate               rad/tu
#
#  outputs       :
#    r           -  position vector            er
#    v           -  velocity vector            er/tu
#
#  locals        :
#    temp        - temporary position vector
#    temp1       - temporary variable
#
#  coupling      :
#    none
#
#  references    :
#    vallado       2001, 246-248, alg 25
#
# [r, v] = radec2rv(rr, rtasc, decl, drr, drtasc, ddecl)
# ------------------------------------------------------------------------------

def radec2rv(rr: float, rtasc: float, decl: float, drr: float, drtasc: float,
             ddecl: float):
    """this function converts the right ascension and declination values with
    position and velocity vectors of a satellite. uses velocity vector to
    find the solution of singular cases.


    Parameters
    ----------
    rr : float
        radius of the satellite: er
    rtasc : float
        right ascension: rad
    decl : float
        declination: rad
    drr : float
        radius of the satellite rate: er/tu
    drtasc : float
        right ascension rate: rad/tu
    ddecl : float
        declination rate: rad/tu

    Returns
    -------
    r : ndarray
        position vector: er
    v : ndarray
        velocity vector: er/tu
    """

    sinrt = math.sin(rtasc)
    cosrt = math.cos(rtasc)
    sindec = math.sin(decl)
    cosdec = math.cos(decl)

    r = np.zeros(3)
    # r[0] = rr*math.cos(decl)*math.cos(rtasc)
    # r[1] = rr*math.cos(decl)*math.sin(rtasc)
    # r[2] = rr*math.sin(decl)
    # r = r.T
    r[0] = rr * cosdec * cosrt
    r[1] = rr * cosdec * sinrt
    r[2] = rr * sindec

    v = np.zeros(3)
    # v[0] = drr*math.cos(decl)*math.cos(rtasc) \
    #     - rr*math.sin(decl)*math.cos(rtasc)*ddecl \
    #     - rr*math.cos(decl)*math.sin(rtasc)*drtasc
    # v[1] = drr*math.cos(decl)*math.sin(rtasc) \
    #     - rr*math.sin(decl)*math.sin(rtasc)*ddecl \
    #     + rr*math.cos(decl)*math.cos(rtasc)*drtasc
    # v[2] = drr*math.sin(decl) + rr*math.cos(decl)*ddecl
    # v = v.T
    v[0] = drr * cosdec * cosrt - rr * sindec * cosrt * ddecl \
        - rr * cosdec * sinrt * drtasc
    v[1] = drr * cosdec * sinrt - rr* sindec * sinrt \
        + rr * cosdec * cosrt * drtasc
    v[2] = drr * sindec + rr * cosdec * ddecl
    return r, v


# ------------------------------------------------------------------------------
#
#                           function rv2radec
#
#  this function converts the right ascension and declination values with
#    position and velocity vectors of a satellite. uses velocity vector to
#    find the solution of singular cases.
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#    vallado     - fix rtasc tests                               29 sep 2002
#
#  inputs          description                    range / units
#    r           -  position vector               km
#    v           -  velocity vector               km/s
#
#  outputs       :
#    rr          - radius of the satellite        km
#    rtasc       - right ascension                rad
#    decl        - declination                    rad
#    drr         - radius of the satellite rate   km/s
#    drtasc      - right ascension rate           rad/s
#    ddecl       - declination rate               rad/s
#
#  locals        :
#    temp        - temporary position vector
#    temp1       - temporary variable
#
#  coupling      :
#    none
#
#  references    :
#    vallado       2001, 246-248, alg 25
#
# [rr, rtasc, decl, drr, drtasc, ddecl] = rv2radec(r, v)
# ------------------------------------------------------------------------------

def rv2radec(r: np.ndarray, v: np.ndarray):
    """this function converts the right ascension and declination values with
    position and velocity vectors of a satellite. uses velocity vector to
    find the solution of singular cases.

    Parameters
    ----------
    r : ndarray
        position vector: km
    v : ndarray
        velocity vector: km/s

    Returns
    -------
    rr: float
        radius of the satellite: km
    rtasc: float
        right ascension: rad
    decl: float
        declination: rad
    drr: float
        radius of the satellite rate: km/s
    drtasc: float
        right ascension rate: rad/s
    ddecl: float
        declination rate: rad/s
    """
    # ------------- calculate angles and rates ----------------
    rr = smu.mag(r)
    temp = np.sqrt(r[0]*r[0] + r[1]*r[1])
    if (temp < small):
        rtasc = math.atan2(v[1] , v[0])
    else:
        rtasc = math.atan2(r[1] , r[0])
    if (rtasc < 0.0):
        rtasc = rtasc + twopi
    decl = math.asin(r[2]/rr)

    temp1 = -r[1]*r[1] - r[0]*r[0]  # different now
    drr = np.dot(r, v)/rr
    if (abs(temp1) > small):
        drtasc = (v[0]*r[1] - v[1]*r[0]) / temp1
    else:
        drtasc = 0.0
    if (abs(temp) > small):
        ddecl = (v[2] - drr*math.sin(decl)) / temp
    else:
        ddecl = 0.0

    return rr, rtasc, decl, drr, drtasc, ddecl


# ------------------------------------------------------------------------------
#
#                           function rv2tradec
#
#  this function converts geocentric equatorial (eci) position and velocity
#    vectors into range, topcentric right acension, declination, and rates.
#    notice the value of small as it can affect the rate term calculations.
#    the solution uses the velocity vector to find the singular cases. also,
#    the right acension and declination rate terms are not observable unless
#    the acceleration vector is available.
#
#  author        : david vallado           davallado@gmail.com    19 jul 2004
#
#  inputs          description                              range / units
#    reci        - eci position vector                      km
#    veci        - eci velocity vector                      km/s
#    latgd       - geodetic latitude                        -pi/2 to pi/2 rad
#    lon         - longitude of site                        -2pi to 2pi rad
#    alt         - altitude                                 km
#    ttt         - julian centuries of tt                   centuries
#    jdut1       - julian date of ut1                       days from 4713 bc
#    lod         - excess length of day                     sec
#    xp          - polar motion coefficient                 rad
#    yp          - polar motion coefficient                 rad
#    terms       - number of terms for ast calculation      0, 2
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    rho         - satellite range from site                km
#    trtasc      - topocentric right ascension              0.0 to 2pi rad
#    tdecl       - topocentric declination                  -pi/2 to pi/2 rad
#    drho        - range rate                               km/s
#    dtrtasc     - topocentric rtasc rate                   rad / s
#    dtdecl      - topocentric decl rate                    rad / s
#
#  locals        :
#    rhoveci     - eci range vector from site               km
#    drhoveci    - eci velocity vector from site            km / s
#
#  coupling      :
#    mag         - magnitude of a vector
#    rot3        - rotation about the 3rd axis
#    rot2        - rotation about the 2nd axis
#
#  references    :
#    vallado       2022, 257, alg 26
#
#  [rho, trtasc, tdecl, drho, dtrtasc, dtdecl] = rv2tradec (reci, veci, latgd, lon, alt, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ------------------------------------------------------------------------------

def rv2tradec(reci: np.ndarray, veci: np.ndarray, latgd: float, lon: float,
              alt: float, ttt: float, jdut1: float, lod: float, xp: float,
              yp: float, terms, ddpsi: float, ddeps: float):
    """this function converts geocentric equatorial (eci) position and velocity
    vectors into range, topcentric right acension, declination, and rates.
    notice the value of small as it can affect the rate term calculations.
    the solution uses the velocity vector to find the singular cases. also,
    the right acension and declination rate terms are not observable unless
    the acceleration vector is available.

    Parameters
    ----------
    reci : ndarray
        position vector eci: km
    veci : ndarray
        velocity vector eci: km/s
    latgd : float
        geodetic latitude: rad
    lon : float
        longitude: rad
    alt : float
        altitude: km
    ttt : float
        centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days from 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    terms :
        NOT USED
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta eps correction to gcrf: rad

    Returns
    -------
    rho : float
        satellite range from site: km
    trtasc : float
        topocentric right ascension: 0.0 to 2pi rad
    tdecl : float
        topocentric declination: -pi/2 to pi/2 rad
    drho : float
        range rate: km/s
    dtrtasc : float
        topocentric rtasc rate: rad/s
    dtdecl : float
        topocentric decl rate: rad/s
    """

    # ----------------- get site vector in ecef -------------------
    rsecef, vsecef = obu.site (latgd, lon, alt)

    #rs
    #vs
    # -------------------- convert ecef to eci --------------------
    a = np.zeros(3)
    rseci, vseci, aeci = ecef2eci(rsecef, vsecef, a, ttt, jdut1, lod, xp, yp,
                                  2, ddpsi, ddeps)
    #rseci
    #vseci

    #rseci = rs
    #vseci = vs
    #[recef, vecef, aecef] = eci2ecef(reci, veci, aeci, ttt, jdut1, lod, xp, yp, 2, 0, 0)
    #reci = recef
    #veci = vecef

    # ------- find eci slant range vector from site to satellite ---------
    rhoveci = reci - rseci
    drhoveci = veci - vseci
    rho = smu.mag(rhoveci)

    # --------------- calculate topocentric rtasc and decl ---------------
    temp = np.sqrt(rhoveci[0] * rhoveci[0] + rhoveci[1] * rhoveci[1])
    if (temp < small):
        trtasc = math.atan2(drhoveci[1], drhoveci[0])
    else:
        trtasc = math.atan2(rhoveci[1], rhoveci[0])

    # directly over the north pole
    if (temp < small):
        tdecl = np.sign(rhoveci[2]) * halfpi   # +- 90 deg
    else:
        magrhoeci = smu.mag(rhoveci)
        tdecl = math.asin(rhoveci[2] / magrhoeci)
    if (trtasc < 0.0):
        trtasc = trtasc + twopi

    # ---------- calculate topcentric rtasc and decl rates -------------
    temp1 = -rhoveci[1] * rhoveci[1] - rhoveci[0] * rhoveci[0]
    drho = np.dot(rhoveci, drhoveci) / rho
    if (abs(temp1) > small):
        dtrtasc = (drhoveci[0]*rhoveci[1] - drhoveci[1] * rhoveci[0]) / temp1
    else:
        dtrtasc = 0.0

    if (abs(temp) > small):
        dtdecl = (drhoveci[2] - drho * math.sin(tdecl)) / temp
    else:
        dtdecl = 0.0

    return rho, trtasc, tdecl, drho, dtrtasc, dtdecl

# ------------------------------------------------------------------------------
#
#                           function rv2razel
#
#  this function converts geocentric equatorial (eci) position and velocity
#    vectors into range, azimuth, elevation, and rates.  notice the value
#    of small as it can affect the rate term calculations. the solution uses
#    the velocity vector to find the singular cases. also, the elevation and
#    azimuth rate terms are not observable unless the acceleration vector is
#    available.
#
#  author        : david vallado                  719-573-2600   22 jun 2002
#
#  revisions
#    vallado     - add terms for ast calculation                 30 sep 2002
#    vallado     - update for site fixes                          2 feb 2004
#
#  inputs          description                    range / units
#    reci        - eci position vector            km
#    veci        - eci velocity vector            km/s
#    rs          - eci site position vector       km
#    latgd       - geodetic latitude              -pi/2 to pi/2 rad
#    lon         - longitude of site              -2pi to 2pi rad
#    alt         - altitude                       km
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    terms       - number of terms for ast calculation 0, 2
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    rho         - satellite range from site      km
#    az          - azimuth                        0.0 to 2pi rad
#    el          - elevation                      -pi/2 to pi/2 rad
#    drho        - range rate                     km/s
#    daz         - azimuth rate                   rad / s
#    del         - elevation rate                 rad / s
#
#  locals        :
#    rhoveci     - eci range vector from site     km
#    drhoveci    - eci velocity vector from site  km / s
#    rhosez      - sez range vector from site     km
#    drhosez     - sez velocity vector from site  km
#    wcrossr     - cross product result           km / s
#    earthrate   - eci earth's rotation rate vec  rad / s
#    tempvec     - temporary vector
#    temp        - temporary real*8 value
#    temp1       - temporary real*8 value
#    i           - index
#
#  coupling      :
#    mag         - magnitude of a vector
#    rot3        - rotation about the 3rd axis
#    rot2        - rotation about the 2nd axis
#
#  references    :
#    vallado       2007, 268-269, alg 27
#
# [rho, az, el, drho, daz, del] = rv2razel (reci, veci, latgd, lon, alt, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ------------------------------------------------------------------------------


def rv2razel(reci: np.ndarray, veci: np.ndarray, latgd: float, lon: float,
             alt: float, ttt: float, jdut1: float, lod: float, xp: float,
             yp: float, terms: int, ddpsi: float, ddeps: float):
    """this function converts geocentric equatorial (eci) position and velocity
    vectors into range, azimuth, elevation, and rates.  notice the value
    of small as it can affect the rate term calculations. the solution uses
    the velocity vector to find the singular cases. also, the elevation and
    azimuth rate terms are not observable unless the acceleration vector is
    available.

    Parameters
    ----------
    reci : ndarray
        eci position vector: km
    veci : ndarray
        eci velocity vector: km/s
    latgd : float
        geodetic latitude of site: rad
    lon : float
        longituge of site: rad
    alt : float
        altitude of site: km
    ttt : float
        julian centuries: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess lenght of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    terms : int
        number of terms for ast calculation: 0, 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta eps correction to gcrf: rad

    Returns
    -------
    rho
        satellite range from site: km
    az
        azimuth: 0.0 to 2pi rad
    el
        elevation: -pi/2 to pi/2 rad
    drho
        range rate: km/s
    daz
        azimuth rate: rad/s
    del
        elevation rate: rad/s
    """

    # ----------------- get site vector in ecef -------------------
    rsecef, _ = obu.site(latgd, lon, alt)
    #print('rsecef    %14.7f %14.7f %14.7f \n', rsecef)

    # -------------------- convert eci to ecef --------------------
    a = np.array([[0],[0],[0]])
    recef, vecef, _ = eci2ecef(reci, veci, a, ttt, jdut1, lod, xp, yp, terms,
                                 ddpsi, ddeps)

    #print('sat recef    %14.7f %14.7f %14.7f \n', recef)
    # simplified - just use sidereal time rotation
    # thetasa = earthrot * (1.0  - 0.0/86400.0)
    # omegaearth = [0; 0; thetasa;]
    # [deltapsi, trueeps, meaneps, omega, nut] = nutation(ttt, ddpsi, ddeps)
    # [st, stdot] = sidereal(jdut1, deltapsi, meaneps, omega, 0, 0)
    #  recef = st'*reci
    #  vecef = st'*veci - cross(omegaearth, recef)


    # ------- find ecef range vector from site to satellite -------
    rhoecef = recef - rsecef
    drhoecef = vecef
    rho = smu.mag(rhoecef)

    # ------------- convert to sez for calculations ---------------
    # tempvec = smu.rot3(rhoecef, lon)
    # rhosez = smu.rot2(tempvec, halfpi-latgd)

    # tempvec = smu.rot3(drhoecef, lon)
    # drhosez = smu.rot2(tempvec, halfpi-latgd)


    # alternate (faster?) sez conversion
    sinlat = math.sin(latgd)
    coslat = math.cos(latgd)
    sinlon = math.sin(lon)
    coslon = math.cos(lon)
    ecef2sez = np.array([[sinlat*coslon, sinlat*sinlon, -coslat],
                         [-sinlon, coslon, 0],
                         [coslat*coslon, coslat*sinlon, sinlat]])
    rhosez = ecef2sez @ rhoecef
    drhosez = ecef2sez @ drhoecef


    # ------------- calculate azimuth and elevation ---------------
    temp = math.sqrt(rhosez[0]*rhosez[0] + rhosez[1]*rhosez[1])

    if (temp < small):           # directly over the north pole
        el = np.sign(rhosez[2])*halfpi   # +- 90 deg
    else:
        magrhosez = smu.mag(rhosez)
        el = math.asin(rhosez[2] / magrhosez)

    if (temp < small):
        dtemp = math.sqrt(drhosez[0]**2 + drhosez[1]**2)
        az = math.atan2(drhosez[1]/dtemp, -drhosez[0]/dtemp)
    else:
        az = math.atan2(rhosez[1]/temp, -rhosez[0]/temp)
    if (az < 0.0):
        az += twopi

    # ------ calculate range, azimuth and elevation rates ---------
    drho = np.dot(rhosez, drhosez)/rho
    if (abs(temp*temp) > small):
        daz = (drhosez[0]*rhosez[1] - drhosez[1]*rhosez[0]) / (temp*temp)
    else:
        daz = 0.0

    if (abs(temp) > small):
        del_ = (drhosez[2] - drho*math.sin(el)) / temp
    else:
        del_ = 0.0
    return rho, az, el, drho, daz, del_

# position and velocity to ecliptic latitude longitude
# dav 28 mar 04
#
# [rr, ecllon, ecllat, drr, decllon, decllat] = rv2ell (rijk, vijk)


def rv2ell (rijk: np.ndarray, vijk: np.ndarray):
    """this function takes position and velocity vectors and turns them into
    ecliptic latitude and longitude

    Parameters
    ----------
    rijk : ndarray
        position vector: km
    vijk : ndarray
        velocity vector: km/s

    Returns
    -------
    rr : float
        range: km
    ecllon : float
        ecliptic longitude: rad
    ecllat: float
        ecliptic latitude: rad
    drr : float
        rate of range: km/sec
    decllon: float
        rate of ecliptic longitude: rad
    decllat: float
        rate of ecliptic latitude: rad
    """

    obliquity = 0.40909280   #23.439291 /rad

    r = smu.rot1 (rijk, obliquity)
    v = smu.rot1 (vijk, obliquity)

    # ------------- calculate angles and rates ----------------
    rr = smu.mag(r)
    temp = np.sqrt(r[0]*r[0] + r[1]*r[1])
    if (temp < small):
        temp1 = np.sqrt(v[0]*v[0] + v[1]*v[1])
        if (abs(temp1) > small):
            ecllon = math.atan2(v[1] , v[0])
        else:
            ecllon = 0.0
    else:
        ecllon = math.atan2(r[1] , r[0])
    ecllat = math.asin(r[2]/rr)

    temp1 = -r[1]*r[1] - r[0]*r[0]  # different now
    drr = np.dot(r, v)/rr
    if (abs(temp1) > small):
        decllon = (v[0]*r[1] - v[1]*r[0]) / temp1
    else:
        decllon = 0.0
    if (abs(temp) > small):
        decllat = (v[2] - drr*math.sin(ecllat)) / temp
    else:
        decllat = 0.0
    return rr, ecllon, ecllat, drr, decllon, decllat

# ecliptic latitude longitude to position and velocity
# dav 28 mar 04
#
# [rijk, vijk] = ell2rv (rr, ecllon, ecllat, drr, decllon, decllat)

def ell2rv(rr: float, ecllon: float, ecllat: float, drr: float,
           decllon: float, decllat: float):
    """this function takes the range and ecliptic latitude and longitude
    and converts them into position and velocity vectors

    Parameters
    ----------
    rr : float
        range: km
    ecllon : float
        ecliptic longitude: rad
    ecllat : float
        ecliptic latidue: rad
    drr : float
        rate of range: km/s
    decllon : float
        rate of ecliptic longitude: rad/s
    decllat : float
        rate of ecliptic latitude: rad/s

    Returns
    -------
    rijk : ndarray
        position vector: km
    vijk : ndarray
        velocity vector: km/s
    """

    obliquity = 0.40909280   #23.439291 /rad

    r = np.zeros((3))
    r[0] = rr*math.cos(ecllat)*math.cos(ecllon)
    r[1] = rr*math.cos(ecllat)*math.sin(ecllon)
    r[2] = rr*math.sin(ecllat)

    v = np.zeros((3))
    v[0] = drr*math.cos(ecllat)*math.cos(ecllon) \
                - rr*math.sin(ecllat)*math.cos(ecllon)*decllat \
                - rr*math.cos(ecllat)*math.sin(ecllon)*decllon
    v[1] = drr*math.cos(ecllat)*math.sin(ecllon) \
                - rr*math.sin(ecllat)*math.sin(ecllon)*decllat \
                + rr*math.cos(ecllat)*math.cos(ecllon)*decllon
    v[2] = drr*math.sin(ecllat) + rr*math.cos(ecllat)*decllat

    rijk = smu.rot1(r, -obliquity)
    vijk = smu.rot1(v, -obliquity)

    return rijk, vijk


# ----------------------------------------------------------------------------
#
#                           function mod2eci
#
#  this function transforms a vector from the mean equator mean equinox of
#    date (mod) to the mean equator mean equinox (j2000) frame.
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    rmod        - position vector of date
#                    mean equator, mean equinox   km
#    vmod        - velocity vector of date
#                    mean equator, mean equinox   km/s
#    amod        - acceleration vector of date
#                    mean equator, mean equinox   km/s2
#    ttt         - julian centuries of tt         centuries
#
#  outputs       :
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#
#  locals        :
#    none.
#
#  coupling      :
#   precess      - rotation for precession        mod - eci
#
#  references    :
#    vallado       2001, 219-220, eq 3-68
#
# [reci, veci, aeci] = mod2eci  (rmod, vmod, amod, ttt)
# ----------------------------------------------------------------------------


def mod2eci(rmod: np.ndarray, vmod: np.ndarray, amod: np.ndarray, ttt: float):
    """this function transforms a vector from the mean equator mean equinox of
    date (mod) to the mean equator mean equinox (j2000) frame.

    Parameters
    ----------
    rmod : ndarray
        position vector mod: km
    vmod : ndarray
        velocity vector mod: km/s
    amod : ndarray
        acceleration vector mod: km/s2
    ttt : float
        julian centuries of tt: centuries

    Returns
    -------
    reci : ndarray
        position vector eci            km
    veci : ndarray
        velocity vector eci            km/s
    aeci : ndarray
        acceleration vector eci        km/s2
    """

    prec, _, _, _, _ = obu.precess(ttt, '80')

    reci = prec@rmod
    veci = prec@vmod
    aeci = prec@amod

    return reci, veci, aeci


# ----------------------------------------------------------------------------
#
#                           function tod2eci
#
#  this function transforms a vector from the true equator true equinox frame
#    of date (tod), to the mean equator mean equinox (j2000) frame.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    rtod        - position vector of date
#                    true equator, true equinox   km
#    vtod        - velocity vector of date
#                    true equator, true equinox   km/s
#    atod        - acceleration vector of date
#                    true equator, true equinox   km/s2
#    ttt         - julian centuries of tt         centuries
#
#  outputs       :
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#
#  locals        :
#    deltapsi    - nutation angle                 rad
#    trueeps     - true obliquity of the ecliptic rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       -                                rad
#    nut         - matrix for mod - tod
#
#  coupling      :
#   precess      - rotation for precession        mod - eci
#   nutation     - rotation for nutation          tod - mod
#
#  references    :
#    vallado       2001, 219-220, eq 3-68
#
# [reci, veci, aeci] = tod2eci  (rtod, vtod, atod, ttt, ddpsi, ddeps)
# ----------------------------------------------------------------------------

def tod2eci(rtod: np.ndarray, vtod: np.ndarray, atod: np.ndarray, ttt: float,
            ddpsi: float, ddeps: float):
    """this function transforms a vector from the true equator true equinox frame
    of date (tod), to the mean equator mean equinox (j2000) frame.

    Parameters
    ----------
    rtod : ndarray
        position vector tod: km
    vtod : ndarray
        velocity vector tod: km/s
    atod : ndarray
        acceleration vector tod: km/s2
    ttt : float
        julian centuries of tt: centuries
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta eps correction to gcrf: rad

    Returns
    -------
    reci : ndarray
        position vector eci: km
    veci: ndarray
        velocity vector eci: km/s
    aeci: ndarray
        acceleration vector eci: km/s2
    """

    prec, _, _, _, _ = obu.precess(ttt, '80')

    _, _, _, _, nut = obu.nutation(ttt, ddpsi, ddeps)

    reci = prec @ nut @ rtod

    veci = prec @ nut @ vtod

    aeci = prec @ nut @ atod

    return reci, veci, aeci


# ----------------------------------------------------------------------------
#
#                           function flt2rv.m
#
#  this function transforms  the flight elements - latgc, lon, fpav, az,
#    position and velocity magnitude into an eci position and velocity vector.
#
#  author        : david vallado                  719-573-2600   17 jun 2002
#
#  revisions
#    vallado     - fix extra terms in rtasc calc                  8 oct 2002
#
#  inputs          description                    range / units
#    magr        - eci position vector magnitude  km
#    magv        - eci velocity vector magnitude  km/sec
#    latgc       - geocentric latitude            rad
#    lon         - longitude                      rad
#    fpa         - sat flight path angle          rad
#    az          - sat flight path az             rad
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    terms       - number of terms for ast calculation 0, 2
#    ddpsi, ddeps - corrections for fk5 to gcrf    rad
#
#  outputs       :
#    r           - eci position vector            km
#    v           - eci velocity vector            km/s
#
#  locals        :
#    fpav        - sat flight path anglefrom vert rad
#
#  coupling      :
#    none        -
#
#  references    :
#    vallado       2013, xx
#    escobal            397
#    chobotov            67
#
# [reci, veci] = flt2rv (magr, magv, latgc, lon, fpa, az, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ----------------------------------------------------------------------------

def flt2rv(magr: float, magv: float, latgc: float, lon:float, fpa: float,
           az: float, ttt: float, jdut1: float, lod: float, xp: float,
           yp: float, terms: int, ddpsi: float, ddeps: float):
    """_this function transforms  the flight elements - latgc, lon, fpav, az,
    position and velocity magnitude into an eci position and velocity vector.

    Parameters
    ----------
    magr : float
        position magnitude: km
    magv : float
        velocity magnitude: km/s
    latgc : float
        geocentric latitude: rad
    lon : float
        longidute: rad
    fpa : float
        satellite flight path angle: rad
    az : float
        satellite azimuth: rad
    ttt : float
        julian centuries of tt: rad
    jdut1 : float
        julian date of ut1: days from 4712 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    terms : int
        # of terms for ast calculation: 0, 2
    ddpsi : float
        delta psi corrections for fk5 to gcrf
    ddeps : float
        delta eps corrections for fk5 to gcrf

    Returns
    -------
    r: ndarray
        eci position vector: km
    v: ndarray
        eci velocity vector: km/s
    """

    small = 1e-08
    # -------- form position vector
    recef = np.zeros(3)
    recef[0] = magr * math.cos(latgc) * math.cos(lon)
    recef[1] = magr * math.cos(latgc) * math.sin(lon)
    recef[2] = magr * math.sin(latgc)
    recef = np.transpose(recef)
    # -------- convert r to eci
    vecef = np.array([[0], [0], [0]])
    aecef = np.array([[0], [0], [0]])
    reci, veci, _ = ecef2eci(recef, vecef, aecef, ttt,
                              jdut1, lod, xp, yp, terms, ddpsi, ddeps)
    # ------------- calculate rtasc and decl ------------------
    temp = math.sqrt(reci[0] * reci[0] + reci[1] * reci[1])
    if (temp < small):
        # v needs to be defined herexxxxxxxxx
        rtasc = math.atan2(veci[1], veci[0])
    else:
        rtasc = math.atan2(reci[1], reci[0])

    decl = math.asin(reci[2] / magr)
    # -------- form velocity vector
    fpav = math.pi * 0.5 - fpa
    veci = np.zeros(3)
    veci[0] = magv * (-math.cos(rtasc) * math.sin(decl)
                      * (math.cos(az) * math.cos(fpav) - math.sin(rtasc)
                         * math.sin(az) * math.cos(fpav))
                      + math.cos(rtasc) * math.sin(decl) * math.sin(fpav))
    veci[1] = magv * (-math.sin(rtasc) * math.sin(decl)
                      * (math.cos(az) * math.cos(fpav) + math.cos(rtasc)
                         * math.sin(az) * math.cos(fpav))
                      + math.sin(rtasc) * math.cos(decl) * math.sin(fpav))
    veci[2] = magv * (math.sin(decl) * math.sin(fpav)
                      + math.cos(decl) * math.cos(az) * math.cos(fpav))
    reci = reci.T
    veci = veci.T
    return reci, veci

# ------------------------------------------------------------------------------
#
#                           function rvs2raz
#
#  this function converts range, azimuth, and elevation values with slant
#    range and velocity vectors for a satellite from a radar site in the
#    topocentric horizon (sez) system.
#
#  author        : david vallado                  719-573-2600   22 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    rhovec      - sez satellite range vector     km
#    drhovec     - sez satellite velocity vector  km / s
#
#  outputs       :
#    rho         - satellite range from site      km
#    az          - azimuth                        0.0 to 2pi rad
#    el          - elevation                      -pi/2 to pi/2 rad
#    drho        - range rate                     km / s
#    daz         - azimuth rate                   rad / s
#    del         - elevation rate                 rad / s
#
#  locals        :
#    sinel       - variable for sin(el)
#    cosel       - variable for cos(el)
#    sinaz       - variable for sin(az)
#    cosaz       - variable for cos(az)
#    temp        -
#    temp1       -
#
#  coupling      :
#    mag         - magnitude of a vector
#
#  references    :
#    vallado       2001, 250-251, eq 4-4, eq 4-5
#
# [rho, az, el, drho, daz, del] = rvs2raz (rhosez, drhosez)
# ------------------------------------------------------------------------------

def rvs2raz(rhosez: np.ndarray, drhosez: np.ndarray):
    """this function converts range, azimuth, and elevation values with slant
    range and velocity vectors for a satellite from a radar site in the
    topocentric horizon (sez) system.

    Parameters
    ----------
    rhosez : ndarray
        position vector sez: km
    drhosez : ndarray
        velocity vector sez: km/s

    Returns
    -------

    rho : float
        satellite range from site: km
    az : float
        azimuth: 0.0 to 2pi rad
    el : float
        elevation: -pi/2 to pi/2 rad
    drho : float
        range rate: km / s
    daz : float
        azimuth rate: rad / s
    del_ : float
        elevation rate: rad / s
    """

    small = 1e-08
    # ------------- calculate azimuth and elevation ---------------
    temp = math.sqrt(rhosez[0] * rhosez[0] + rhosez[1] * rhosez[1])
    if (abs(rhosez[1]) < small):
        if (temp < small):
            az = math.atan2(drhosez[1], -drhosez[0])
        elif (rhosez[0] > 0.0):
            az = math.pi
        else:
            az = 0.0
    else:
        az = math.atan2(rhosez[1], -rhosez[0])

    rho = smu.mag(rhosez)


    if ((temp < small)):
        el = np.sign(rhosez[2]) * halfpi
    else:
        el = math.asin(rhosez[2] / rho)

    # -------  calculate range, azimuth and elevation rates -------
    drho = np.dot(rhosez, drhosez.T) / rho
    if (abs(temp * temp) > small):
        daz = (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) / (temp * temp)
    else:
        daz = 0.0

    if (abs(temp) > small):
        del_ = (drhosez[2] - drho * math.sin(el)) / temp
    else:
        del_ = 0.0

    return rho, az, el, drho, daz, del_


#
# convert orbital elments and find nu
#  input all angles in rad
#  output is in rad
#  dav 6 may 2011
#
# jdut1 = jday(2011, 3, 22, 12, 30, 0)
# incl = 0.04801366433780/rad
# raan = 330.034263262230/rad
# argp = 91.3766761083524/rad
# lon = -7.2935047164727/rad
#
# jdut1 = jday(2014, 7, 26, 0, 0, 0.0)
# nu = lon2nu (jdut1, 7.020438698/rad, 0.070273056/rad, 19.90450011/rad, 352.5056022/rad)
#

def lon2nu(jdut1=None, lon=None, incl=None, raan=None, argp=None):
    # fprintf(' jd #16.8f lon #11.5f  incl #11.5f raan #11.5f argp #11.5f \n', jdut1, lon * rad2deg, incl * rad2deg, raan * rad2deg, argp * rad2deg)
    # need to use their GMST calculation
    ed = jdut1 + 0.0 - 2451544.5

    gmst = 99.96779469 + 360.985647366286 * ed + 2.9079e-13 * ed * ed

    gmst = np.fmod(gmst * deg2rad, 2.0 * np.pi)
    # ------------------------ check quadrants --------------------
    if (gmst < 0.0):
        gmst = gmst + 2.0 * np.pi

    lambdau = gmst + lon - raan
    # make sure lambdau is 0 to 360 deg
    if lambdau < 0.0:
        lambdau = lambdau + 2.0 * np.pi

    if lambdau > twopi:
        lambdau = lambdau - 2.0 * np.pi

    arglat = np.arctan(np.tan(lambdau) / np.cos(incl))
    # find nu
    if (lambdau >= 0.5 * np.pi) and (lambdau < 1.5 * np.pi):
        arglat = arglat + np.pi

    temp = arglat - argp
    nu = temp
    # fprintf(' #11.5f #11.5f #11.5f #11.5f  #16.10f ', lambdau * rad2deg, argp * rad2deg, lon * rad2deg, gmst * rad2deg, temp * rad2deg)
    print(' lu %11.5f argp %11.5f lon %11.5f gmst %11.5f arglat %11.5f nu %16.10f '
          % (lambdau * rad2deg, argp * rad2deg, lon * rad2deg, gmst * rad2deg, arglat * rad2deg, temp * rad2deg))
    #     fprintf(' nu = #11.5f deg \n', nu * rad2deg)
    return nu

#
# convert orbital elments and find lon
#  input all angles in rad
#  output is in rad
#  dav 6 may 2011
#
# jdut1 = jday(2011, 3, 22, 12, 30, 0)
# incl = 0.04801366433780/rad
# raan = 330.034263262230/rad
# argp = 91.3766761083524/rad
# lon = -7.2935047164727/rad
#

def nu2lon(jdut1=None, nu=None, incl=None, raan=None, argp=None):
    #    fprintf(' jd #16.8f lon #11.5f  incl #11.5f raan #11.5f argp #11.5f \n', jdut1, lon * rad2deg, incl * rad2deg, raan * rad2deg, argp * rad2deg)
#    need to use their GMST calculation
    ed = jdut1 - 2451544.5

    gmst = 99.96779469 + 360.985647366286 * ed + 2.9079e-13 * ed * ed

    gmst = np.fmod(gmst * deg2rad, 2.0 * np.pi)
    # ------------------------ check quadrants --------------------
    if (gmst < 0.0):
        gmst = gmst + 2.0 * np.pi

    arglat = nu + argp
    # make sure lambdau is 0 to 360 deg
    if arglat < 0.0:
        arglat = arglat + 2.0 * np.pi

    if arglat > twopi:
        arglat = arglat - 2.0 * np.pi

    lambdau = np.arctan(np.tan(arglat) * np.cos(incl))
    if (arglat >= 0.5 * np.pi) and (arglat < 1.5 * np.pi):
        lambdau = lambdau + np.pi

    temp = lambdau - gmst + raan
    # fprintf(' xx #11.5f  #11.5f #11.5f #11.5f ', lambdau * rad2deg, arglat * rad2deg, raan * rad2deg, temp * rad2deg)
# make sure lambdau is 0 to 360 deg
    temp = np.fmod(temp, 2.0 * np.pi)
    lon = temp
    # fprintf(' #11.5f  #11.5f ', nu * rad2deg, raan * rad2deg)

    #     fprintf(' nu = #11.5f deg \n', nu * rad2deg)

    return lon

#
# ----------------------------------------------------------------------------
#
#                           function pef2eci
#
#  this function trsnforms a vector from the pseudo earth fixed frame (pef),
#    to the mean equator mean equinox (j2000) frame.
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#    vallado     - add terms for ast calculation                 30 sep 2002
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    rpef        - position pseudo earth fixed    km
#    vpef        - velocity pseudo earth fixed    km/s
#    apef        - acceleration pseudo earth fixedkm/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    terms       - number of terms for ast calculation 0, 2
#
#
#
#  outputs       :
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#
#  locals        :
#    prec        - matrix for eci - mod
#    deltapsi    - nutation angle                 rad
#    trueeps     - true obliquity of the ecliptic rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       -                                rad
#    nut         - matrix for mod - tod
#    st          - matrix for tod - pef
#    stdot       - matrix for tod - pef rate
#
#  coupling      :
#   precess      - rotation for precession        mod - eci
#   nutation     - rotation for nutation          tod - mod
#   sidereal     - rotation for sidereal time     pef - tod
#
#  references    :
#    vallado       2001, 219-220, eq 3-68
#
# [reci, veci, aeci] = pef2eci  (rpef, vpef, apef, ttt, jdut1, lod, eqeterms, ddpsi, ddeps)
# ----------------------------------------------------------------------------

def pef2eci(rpef: np.ndarray, vpef: np.ndarray, apef: np.ndarray, ttt: float,
            jdut1: float, lod: float, eqeterms: int, ddpsi: float,
            ddeps: float):
    """this function trsnforms a vector from the pseudo earth fixed frame (pef),
    to the mean equator mean equinox (j2000) frame.

    Parameters
    ----------
    rpef : ndarray
        position vector pef: km
    vpef : ndarray
        velocity vector pef: km/s
    apef : np.ndarray
        acceleration vector pef: km/s2
    ttt : float
        julian centuries from tt: centuries
    jdut1 : float
        julian date from ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    eqeterms : int
        number of terms for ast calculation: 0 or 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta eps correction to gcrf: rad

    Returns
    -------
    reci: ndarray
        position vector eci: km
    veci: ndarray
        velocity vector eci: km/s
    aeci: ndarray
        acceleration vector eci: km/s2
    """

    prec, _, _, _, _ = obu.precess(ttt, '80')
    deltapsi, _, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)
    st, _ = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)
    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    reci = prec @ nut @ st @ rpef
    veci = prec @ nut @ st @ (vpef + np.cross(omegaearth, rpef.T).T)
    temp = np.cross(omegaearth, rpef.T)
    aeci = prec @ nut @ st @ (apef + np.cross(omegaearth, temp).T + 2.0
                              * np.cross(omegaearth, vpef.T).T)
    return reci, veci, aeci

# -----------------------------------------------------------------------------
#
#                           function tradec2rv
#
#  this function converts range, topcentric right acension, declination, and rates
#    into geocentric equatorial (eci) position and velocity vectors.
#
#  author        : david vallado           davallado@gmail.com    4 nov 2022
#
#  revisions
#
#  inputs          description                              range / units
#    rho         - satellite range from site                km
#    trtasc      - topocentric right ascension              0.0 to 2pi rad
#    tdecl       - topocentric declination                  -pi/2 to pi/2 rad
#    drho        - range rate                               km/s
#    dtrtasc     - topocentric rtasc rate                   rad / s
#    dtdecl      - topocentric decl rate                    rad / s
#    rseci       - eci site position vector                 km
#    vseci       - eci site velocity vector                 km/s
#    lod         - excess length of day                     sec
#
#  outputs       :
#    reci        - eci position vector                      km
#    veci        - eci velocity vector                      km/s
#
#  locals        :
#    rhov        - eci range vector from site               km
#    drhov       - eci velocity vector from site            km / s
#    omegaearth  - eci earth's rotation rate vec            rad / s
#    tempvec     - temporary vector
#    latgc       - site geocentric latitude                 rad
#
#  coupling      :
#    mag         - magnitude of a vector
#    rot3        - rotation about the 3rd axis
#    rot2        - rotation about the 2nd axis
#
#  references    :
#    vallado       2022, 254, eq 4-1 to 4-2
#
# [reci, veci] = tradec2rv (rho, trtasc, tdecl, drho, dtrtasc, dtdecl, rseci, lod)
# ------------------------------------------------------------------------------

def tradec2rv(rho: float, trtasc: float, tdecl: float, drho: float,
              dtrtasc: float, dtdecl: float, rseci: np.ndarray,
              vseci: np.ndarray, lod: float):
    """this function converts range, topcentric right acension, declination,
    and rates into geocentric equatorial (eci) position and velocity vectors.

    Parameters
    ----------
    rho : float
        satellite range from site: km
    trtasc : float
        topocentric right angle of ascension: rad
    tdecl : float
        topocentric declination: rad
    drho : float
        satellite range rate: km/s
    dtrtasc : float
        topocentric right angle of ascension rate: rad/s
    dtdecl : float
        topocentric declination rate: rad/s
    rseci : ndarray
        eci position vector of site: km
    vseci : ndarray
        eci velocity vector of site: km
    lod : float
        excess length of day: sec

    Returns
    -------
    reci: ndarray
        position vector of satellite eci: km
    veci: ndarray
        velocity vector of satellite eci: km/s
    """

    latgc = math.asin(rseci(3) / smu.mag(rseci))
    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    np.cross(omegaearth, rseci, vseci)

    # --------  calculate topocentric slant range vectors ------------------
    rhov = np.zeros(3)
    drhov = np.zeros(3)

    rhov[0] = rho * math.cos(tdecl) * math.cos(trtasc)
    rhov[1] = rho * math.cos(tdecl) * math.sin(trtasc)
    rhov[2] = rho * math.sin(tdecl)

    drhov[0] = (drho * math.cos(tdecl) * math.cos(trtasc) - rho * math.sin(tdecl)
                * math.cos(trtasc) * dtdecl
                - rho * math.cos(tdecl) * math.sin(trtasc) * dtrtasc)
    drhov[1] = (drho * math.cos(tdecl) * math.sin(trtasc)
                - rho * math.sin(tdecl) * math.sin(trtasc) * dtdecl
                + rho * math.cos(tdecl) * math.cos(trtasc) * dtrtasc)
    drhov[2] = drho * math.sin(tdecl) + rho * math.cos(tdecl) * dtdecl

    # ------ find eci range vector from site to satellite ------
    reci = rhov + rseci
    veci = drhov + math.cos(latgc) * vseci
    return reci, veci

# ----------------------------------------------------------------------------
#
#                           function teme2eci
#
#  this function transforms a vector from the true equator mean equinox system,
#    (teme) to the mean equator mean equinox (j2000) system.
#
#  author        : david vallado                  719-573-2600   30 oct 2017
#
#  inputs          description                    range / units
#    rteme       - position vector of date
#                    true equator, mean equinox   km
#    vteme       - velocity vector of date
#                    true equator, mean equinox   km/s
#    ateme       - acceleration vector of date (if available, else set to 0)
#                    true equator, mean equinox   km/s2
#    ttt         - julian centuries of tt         centuries
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#
#  locals        :
#    prec        - matrix for eci - mod
#    nutteme     - matrix for mod - teme - an approximation for nutation
#    eqeg        - rotation for equation of equinoxes (geometric terms only)
#    tm          - combined matrix for teme2eci
#
#  coupling      :
#   precess      - rotation for precession        eci - mod
#   nutation     - rotation for nutation          eci - tod
#
#  references    :
#    vallado       2013, 231-233
#
# [reci, veci, aeci] = teme2eci  (rteme, vteme, ateme, ttt, ddpsi, ddeps)
# ----------------------------------------------------------------------------

def teme2eci(rteme: np.ndarray, vteme: np.ndarray, ateme: np.ndarray,
             ttt: float, ddpsi: float, ddeps: float):
    """this function transforms a vector from the true equator mean equinox system,
    (teme) to the mean equator mean equinox (j2000) system.

    Parameters
    ----------
    rteme : ndarray
        position vector of date
        true equator, mean equinox: km
    vteme : ndarray
        velocity vector of date
        true equator, mean equinox: km
    ateme : ndarray
        acceleration vector of date
        true equator, mean equinox: km
    ttt : float
        julian centuries of tt: centuries
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta eps correction to gcrf: rad

    Returns
    -------
    reci: ndarray
        position vector eci: km
    veci: ndarray
        velocity vector eci: km/s
    aeci: ndarray
        acceleration vector eci: km/s2
    """
    prec, _, _, _, _ = obu.precess(ttt, '80')
    deltapsi, _, meaneps, _, nut = obu.nutation(ttt, ddpsi, ddeps)
    # ------------------------ find eqeg ----------------------
# rotate teme through just geometric terms
    eqeg = deltapsi * np.cos(meaneps)
    eqeg = np.fmod(eqeg, 2.0 * np.pi)
    eqe = np.zeros((3, 3))
    eqe[0, 0] = np.cos(eqeg)
    eqe[0, 1] = np.sin(eqeg)
    eqe[0, 2] = 0.0
    eqe[1, 0] = - np.sin(eqeg)
    eqe[1, 1] = np.cos(eqeg)
    eqe[1, 2] = 0.0
    eqe[2, 0] = 0.0
    eqe[2, 1] = 0.0
    eqe[2, 2] = 1.0
    tm = prec @ nut @ eqe.T
    reci = tm @ rteme
    veci = tm @ vteme
    aeci = tm @ ateme
    return reci, veci, aeci

# ----------------------------------------------------------------------------
#
#                           function teme2ecef
#
#  this function trsnforms a vector from the true equator mean equniox frame
#    (teme), to an earth fixed (ITRF) frame.  the results take into account
#    the effects of sidereal time, and polar motion.
#
#  author        : david vallado                  719-573-2600   30 oct 2017
#
#  revisions
#
#  inputs          description                    range / units
#    rteme       - position vector teme           km
#    vteme       - velocity vector teme           km/s
#    ateme       - acceleration vector teme (if available, else set to 0)      km/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    eqeterms    - use extra two terms (kinematic) after 1997  0, 2
#
#  outputs       :
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#
#  locals        :
#    st          - matrix for pef - tod
#    pm          - matrix for ecef - pef
#
#  coupling      :
#   gstime       - greenwich mean sidereal time   rad
#   polarm       - rotation for polar motion      pef - ecef
#
#  references    :
#    vallado       2013, 231-233
#
# [recef, vecef, aecef] = teme2ecef(rteme, vteme, ateme, ttt, jdut1, lod, xp, yp, eqeterms)
# ----------------------------------------------------------------------------

def teme2ecef(rteme: np.ndarray, vteme: np.ndarray, ateme: np.ndarray,
              ttt: float, jdut1: float, lod: float, xp: float, yp: float,
              eqeterms: int):
    """this function trsnforms a vector from the true equator mean equniox frame
    (teme), to an earth fixed (ITRF) frame.  the results take into account
    the effects of sidereal time, and polar motion.

    Parameters
    ----------
    rteme : ndarray
        position vector teme: km
    vteme : ndarray
        velocity vector teme: km/s
    ateme : ndarray
        acceleration vector teme: km/s2
    ttt : float
        julain centuries of tt: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    eqeterms : int
        use 2 extra terms (kinematic) after 1997: 0, 2

    Returns
    -------
    recef: ndarray
        position vector earth fixed: km
    vecef: ndarray
        velocity vector earth fixed: km/s
    aecef: ndarray
        acceleration vector earth fixed: km/s2
    """

    # ------------------------ find gmst --------------------------
    gmst = stu.gstime(jdut1)
    # find omega from nutation theory
    omega = 125.04452222 + (- 6962890.539 * ttt + 7.455 * ttt
                            * ttt + 0.008 * ttt * ttt * ttt) / 3600.0
    omega = np.fmod(omega, 360.0) * deg2rad
    # ------------------------ find mean ast ----------------------
    # teme does not include the geometric terms here
    # after 1997, kinematic terms apply
    if (jdut1 > 2450449.5) and (eqeterms > 0):
        gmstg = (gmst + 0.00264 * arcsec2rad * np.sin(omega)
                 + 6.3e-05 * arcsec2rad * np.sin(2.0 * omega))
    else:
        gmstg = gmst

    gmstg = np.fmod(gmstg, 2.0 * np.pi)
    st = np.zeros((3, 3))
    st[0, 0] = np.cos(gmstg)
    st[0, 1] = - np.sin(gmstg)
    st[0, 2] = 0.0
    st[1, 0] = np.sin(gmstg)
    st[1, 1] = np.cos(gmstg)
    st[1, 2] = 0.0
    st[2, 0] = 0.0
    st[2, 1] = 0.0
    st[2, 2] = 1.0
    pm = smu.polarm(xp, yp, ttt, '80')
    rpef = st.T@rteme
    recef = pm.T@rpef
    thetasa = earthrot * (1.0  - lod/86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])

    vpef = st.T@vteme - np.cross(omegaearth, rpef.T).T
    vecef = pm.T@vpef
    temp = np.cross(omegaearth, rpef.T)
    aecef = pm.T@(st.T@ateme - np.cross(omegaearth, temp).T \
            - 2.0*np.cross(omegaearth, vpef.T))

    #fprintf(1, 'st gmst #11.8f ast #11.8f ome  #11.8f \n', gmst*180/pi, ast*180/pi, omegaearth*180/pi)
    return recef, vecef, aecef

# ----------------------------------------------------------------------------
#
#                           function eci2pef
#
#  this function transforms a vector from the mean equator, mean equinox frame
#    (j2000), to the pseudo earth fixed frame (pef).
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#    vallado     - add terms for ast calculation                 30 sep 2002
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    terms       - number of terms for ast calculation 0, 2
#    ddpsi
#    ddeps
#
#  outputs       :
#    rpef        - position pseudo earth fixed    km
#    vpef        - velocity pseudo earth fixed    km/s
#    apef        - acceleration pseudo earth fixedkm/s2
#
#  locals        :
#    prec        - matrix for eci - mod
#    deltapsi    - nutation angle                 rad
#    trueeps     - true obliquity of the ecliptic rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       -                                rad
#    nut         - matrix for mod - tod
#    st          - matrix for tod - pef
#    stdot       - matrix for tod - pef rate
#
#  coupling      :
#   precess      - rotation for precession        mod - eci
#   nutation     - rotation for nutation          tod - mod
#   sidereal     - rotation for sidereal time     pef - tod
#
#  references    :
#    vallado       2001, 219, eq 3-65 to 3-66
#
# [rpef, vpef, apef] = eci2pef  (reci, veci, aeci, opt, ttt, jdut1, lod, eqeterms, ddpsi, ddeps)
# ----------------------------------------------------------------------------

def eci2pef(reci: np.ndarray, veci: np.ndarray, aeci: np.ndarray,
            ttt: float, jdut1: float, lod: float, terms: int, ddpsi: float,
            ddeps: float):
    """this function transforms a vector from the mean equator, mean equinox frame
    (j2000), to the pseudo earth fixed frame (pef).

    Parameters
    ----------
    reci : ndarray
        position vector eci: km
    veci : ndarray
        velocity vector eci: km/s
    aeci : ndarray
        acceleration vector eci: kms/2
    ttt : float
        julian centuries of  tt: centuries
    jdut1 : float
        julian date of ut1: days from 4713 bc
    lod : float
        excess length of day: sec
    terms : int
        # of terms for ast calculation: 0, 2
    ddpsi : float
        delta psi correction to gcrf: rad
    ddeps : float
        delta psi correction to gcrf: rad

    Returns
    -------
    rpef: ndarray
        position vector pseudo earth fixed: km
    vpef: ndarray
        velocity vector pseudo earth fixed: km/s
    apef: ndarray
        acceleration vector pseudo earth fixed: km/s2
    """

    prec, _, _, _, _ = obu.precess(ttt, '80')
    deltapsi, _, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)
    st, _ = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, terms)


    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    rpef = st.T@nut.T@prec.T@reci
    vpef = st.T@nut.T@prec.T@veci - np.cross(omegaearth, rpef.T).T

    temp = np.cross(omegaearth, rpef.T)
    apef = st.T@nut.T@prec.T@aeci - np.cross(omegaearth, temp).T \
        - 2.0*np.cross(omegaearth, vpef.T).T

    return rpef, vpef, apef

def eci2tirsiau06(reci: np.ndarray, veci: np.ndarray, aeci: np.ndarray,
                  opt: str, ttt: float, jdut1: float, lod: float,
                  ddx: float = None, ddy: float = None):
    """this function transforms a vector from the mean equator, mean equinox
    frame (GCRF) to the Terrestrial Intermediate Reference System (TIRS).

    Parameters
    ----------
    reci : ndarray
        position vector eci: km
    veci : ndarray
        velocity vector eci: km/s
    aeci : ndarray
        acceleration vector eci: km/s2
    opt : str
        "c": cio based, iau2006
        "a": class equinox based, 2000a
        "b": class equinox based, 2000b
    ttt : float
        julian centuries of date: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    ddx : float, optional
        coordinates of the Celestial Intermediate Pole,
        needed for opt 'c'
    ddy : float, optional
        coordinates of the Celestial Intermediate Pole,
        needed for opt 'c'

    Returns
    -------
    rtirs: ndarray
        position vector tirs: km
    vtirs: ndarray
        velocity vector tirs: km/s
    atirs: ndarray
        acceleration vector tirs: km/s2
    """

     # ---- cio based, iau2006
    if opt == 'c':
        _, _, _, pnb = obu.iau06xys(ttt, ddx, ddy)
        st = obu.iau06era(jdut1)
    # ---- class equinox based, 2000a
    if opt == 'a':
        deltapsi, pnb, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ \
            = obu.iau06pna(ttt)
        _, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')
    # ---- class equinox based, 2000b
    if opt == 'b':
        deltapsi, pnb, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ \
            = obu.iau06pnb(ttt)
        _, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')

    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    rtirs = st.T @ pnb.T @ reci
    vtirs = st.T @ pnb.T @ veci - np.cross(omegaearth, rtirs.T).T

    temp = np.cross(omegaearth, rtirs.T)
    atirs = st.T @ pnb.T @ aeci - np.cross(omegaearth, temp).T \
        - 2.0*np.cross(omegaearth, vtirs.T).T

    return rtirs, vtirs, atirs

def tirs2eciiau06(rtirs: np.ndarray, vtirs: np.ndarray, atirs: np.ndarray,
                  opt: str, ttt: float, jdut1: float, lod: float,
                  ddx: float = None, ddy: float = None):
    """this function transforms a vector from the mean equator, mean equinox
    frame (GCRF) to the Terrestrial Intermediate Reference System (TIRS).

    Parameters
    ----------
    rtirs : ndarray
        position vector eci: km
    vtirs : ndarray
        velocity vector eci: km/s
    atirs : ndarray
        acceleration vector eci: km/s2
    opt : str
        "c": cio based, iau2006
        "a": class equinox based, 2000a
        "b": class equinox based, 2000b
    ttt : float
        julian centuries of date: centuries
    jdut1 : float
        julian date of ut1: days since 4713 bc
    lod : float
        excess length of day: sec
    ddx : float, optional
        coordinates of the Celestial Intermediate Pole,
        needed for opt 'c'
    ddy : float, optional
        coordinates of the Celestial Intermediate Pole,
        needed for opt 'c'

    Returns
    -------
    reci: ndarray
        position vector eci: km
    veci: ndarray
        velocity vector eci: km/s
    aeci: ndarray
        acceleration vector eci: km/s2
    """

     # ---- cio based, iau2006
    if opt == 'c':
        _, _, _, pnb = obu.iau06xys(ttt, ddx, ddy)
        st = obu.iau06era(jdut1)
    # ---- class equinox based, 2000a
    if opt == 'a':
        deltapsi, pnb, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ \
            = obu.iau06pna(ttt)
        _, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')
    # ---- class equinox based, 2000b
    if opt == 'b':
        deltapsi, pnb, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ \
            = obu.iau06pnb(ttt)
        _, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')

    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    reci = pnb @ st @ rtirs
    veci = pnb @ st @ (vtirs + np.cross(omegaearth, rtirs.T).T)

    temp = np.cross(omegaearth, rtirs.T)
    aeci = pnb @ st @ (atirs + np.cross(omegaearth, temp).T \
        + 2.0*np.cross(omegaearth, vtirs.T).T)

    return reci, veci, aeci

# ----------------------------------------------------------------------------
#
#                           function cirs2ecefiau06
#
#  this function trsnforms a vector from the cirs
#    (gcrf), to an earth fixed (itrf) frame.  the results take into account
#    the effects of  sidereal time, and polar motion.
#
#  author        : david vallado                  719-573-2600   2 may 2020
#
#  revisions
#
#  inputs          description                    range / units
#    rcirs       - position vector cirs            km
#    vcirs       - velocity vector cirs            km/s
#    acirs       - acceleration vector cirs        km/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    option      - which approach to use          a-2000a, b-2000b, c-2000xys
#
#  outputs       :
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#
#  locals        :
#    pm          - transformation matrix for itrf-pef
#    st          - transformation matrix for pef-ire
#
#  coupling      :
#   iau00pm      - rotation for polar motion      itrf-pef
#   iau00era     - rotation for earth rotation    pef-ire
#
#  references    :
#    vallado       2004, 205-219
#
# [recef, vecef, aecef] = cirs2ecefiau06  (rcirs, vcirs, acirs, ttt, jdut1, lod, xp, yp, option, ddx, ddy)
# ----------------------------------------------------------------------------

def cirs2ecefiau06(rcirs: np.ndarray, vcirs:np.ndarray, acirs:np.ndarray,
                   ttt: float, jdut1: float, lod: float, xp: float, yp: float,
                   option: str, ddx: float, ddy: float):
    """this function trsnforms a vector from the cirs
    (gcrf), to an earth fixed (itrf) frame.  the results take into account
    the effects of  sidereal time, and polar motion.

    Parameters
    ----------
    rcirs : ndarray
        position vector cirs: km
    vcirs : ndarray
        velocity vector cirs: km/s
    acirs : ndarray
        acceleration vector cirs: km/s2
    ttt : float
        julian centuries of date: centuries
    jdut1 : float
        julian date of ut1: days from 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    option : str
        approach to use: 'a' - classical equinox 2000a, 'b' - classical equinox
        2000b, 'c' - cio iau2006
    ddx : float
        eop correction for x: rad
    ddy : float
        eop correction for y: rad

    Returns
    -------
    recef: ndarray
        position vector ecef: km
    vecef: ndarray
        velocity vector ecef: km/s
    aecef: ndarray
        acceleration vector ecef: km/s2
    """

    # ---- cio based, iau2006
    if option == 'c':
        x, y, s, pnb = obu.iau06xys(ttt, ddx, ddy)
        st = obu.iau06era(jdut1)
    # ---- class equinox based, 2000a
    if option == 'a':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pna(ttt)
        gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')
    # ---- class equinox based, 2000b
    if option == 'b':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pnb(ttt)
        gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')

    pm = smu.polarm(xp, yp, ttt, '06')
    # ---- setup parameters for velocity transformations
    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])

    rpef = st.T@rcirs
    recef = pm.T@rpef

    vpef = st.T@vcirs - np.cross(omegaearth, rpef.T).T
    vecef = pm.T@vpef

    temp = np.cross(omegaearth, rpef.T)
    aecef = pm.T@(st.T@acirs - np.cross(omegaearth, temp).T - 2.0*np.cross(omegaearth, vpef.T).T)

    return recef, vecef, aecef

# ----------------------------------------------------------------------------
#
#                           function cirs2eciiau06
#
#  this function transforms a vector from the cirs frame, to
#    the eci mean equator mean equinox (gcrf).
#
#  author        : david vallado                  719-573-2600   2 may 2020
#
#  revisions
#
#  inputs          description                    range / units
#    rcirs       - position vector earth fixed    km
#    vcirs       - velocity vector earth fixed    km/s
#    acirs       - acceleration vector earth fixedkm/s2
#    ttt         - julian centuries of tt         centuries
#    option      - which approach to use          a-2000a, b-2000b, c-2000xys
#    ddx         - eop correction for x           rad
#    ddy         - eop correction for y           rad
#
#  outputs       :
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#
#  locals        :
#    nut         - transformation matrix for ire-gcrf
#
#  coupling      :
#   iau00era     - rotation for earth rotyation   pef-ire
#   iau00xys     - rotation for prec/nut          ire-gcrf

#
#  references    :
#    vallado       2004, 205-219
#
# [reci, veci, aeci] = cirs2eciiau06 (rcirs, vcirs, acirs, ttt, option, ddx, ddy)
# ----------------------------------------------------------------------------

def cirs2eciiau06(rcirs: np.ndarray, vcirs: np.ndarray, acirs: np.ndarray,
                  ttt: float, option: str, ddx: float, ddy: float):
    """this function transforms a vector from the cirs frame, to
    the eci mean equator mean equinox (gcrf).

    Parameters
    ----------
    rcirs : ndarray
        position vector cirs: km
    vcirs : ndarray
        velocity vector cirs: km/s
    acirs : ndarray
        acceleration vector cirs: km/s2
    ttt : float
        julian centuries of date: centuries
    option : str
        approach to use: 'a' - classical equinox 2000a, 'b' - classical equinox
        2000b, 'c' - cio iau2006
    ddx : float
        eop correction for x: rad
    ddy : float
        eop correction for y: rad

    Returns
    -------
    reci: ndarray
        position vector eci: km
    veci: ndarray
        velocity vector eci: km/s
    aeci: ndarray
        acceleration vector eci: km/s2
    """
    # ---- ceo based, iau2006
    if option == 'c':
        x, y, s, pnb = obu.iau06xys(ttt, ddx, ddy)

    # ---- class equinox based, 2000a
    if option == 'a':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pna(ttt)

    # ---- class equinox based, 2000b
    if option == 'b':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pnb(ttt)

    # ---- perform transformations
    reci = pnb @ rcirs
    veci = pnb @ vcirs
    aeci = pnb @ acirs
    return reci, veci, aeci

# ----------------------------------------------------------------------------
#
#                           function eci2cirsiau06
#
#  this function trsnforms a vector from the mean equator mean equniox frame
#    (gcrf), to the CIRS frame.  the results take into account
#    the effects of precession, nutation.
#
#  author        : david vallado                  719-573-2600    2 may 2020
#
#  revisions
#
#  inputs          description                    range / units
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#    ttt         - julian centuries of tt         centuries
#    option      - which approach to use          a-2000a, b-2000b, c-2000xys
#    ddx         - eop correction for x           rad
#    ddy         - eop correction for y           rad
#
#  outputs       :
#    rcirs       - position vector earth fixed    km
#    vcirs       - velocity vector earth fixed    km/s
#    acirs       - acceleration vector earth fixedkm/s2
#
#  locals        :
#    nut         - transformation matrix for ire-gcrf
#
#  coupling      :
#   iau00era     - rotation for earth rotation    pef-ire
#   iau00xys     - rotation for prec/nut          ire-gcrf
#
#  references    :
#    vallado       2004, 205-219
#
# [rcirs, vcirs, acirs] = eci2cirsiau06  (reci, veci, aeci, ttt, option, ddx, ddy)
# ----------------------------------------------------------------------------

def eci2cirsiau06(reci: np.ndarray, veci: np.ndarray, aeci:np.ndarray,
                  ttt: float, option: str, ddx: float = None,
                  ddy: float = None):
    """this function transforms a vector from the mean equator mean equniox
    frame (gcrf), to the CIRS frame. the results take into account
    the effects of precession, nutation.

    Parameters
    ----------
    reci : ndarray
        position array eci: km
    veci : ndarray
        velocity vector eci: km/s
    aeci : ndarray
        acceleration vector eci: km/s2
    ttt : float
        julian centuries of date: centuries
    option : str
        approach to use: 'a' - classical equinox 2000a, 'b' - classical equinox
        2000b, 'c' - cio iau2006
    ddx : float, optional
        eop correction for x: rad
            only needed for option 'c'
    ddy : float, optional
        eop correction for y: rad
            only needed for option 'c'

    Returns
    -------
    rcirs: ndarray
        position vector cirs: km
    vcirs: ndarray
        velocity vector cirs: km/s
    acirs: ndarray
        acceleration vector: km/s2
    """

    # ---- cio based, iau2000
    if option == 'c':
        x, y, s, pnb = obu.iau06xys(ttt, ddx, ddy)

    # ---- class equinox based, 2000a
    if option == 'a':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, \
            lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pna(ttt)

    # ---- class equinox based, 2000b
    if option == 'b':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, \
            lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pnb(ttt)
    rcirs = pnb.T @ reci
    vcirs = pnb.T @ veci
    acirs = pnb.T @ aeci
    return rcirs, vcirs, acirs

# ----------------------------------------------------------------------------
#
#                           function eci2ecefiau06
#
#  this function trsnforms a vector from the mean equator mean equniox frame
#    (gcrf), to an earth fixed (itrf) frame.  the results take into account
#    the effects of precession, nutation, sidereal time, and polar motion.
#
#  author        : david vallado                  719-573-2600   16 jul 2004
#
#  revisions
#
#  inputs          description                    range / units
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    option      - which approach to use          a-2000a, b-2000b, c-2000xys
#    ddx         - eop correction for x           rad
#    ddy         - eop correction for y           rad
#
#  outputs       :
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#
#  locals        :
#    pm          - transformation matrix for itrf-pef
#    st          - transformation matrix for pef-ire
#    nut         - transformation matrix for ire-gcrf
#
#  coupling      :
#   iau00pm      - rotation for polar motion      itrf-pef
#   iau00era     - rotation for earth rotation    pef-ire
#   iau00xys     - rotation for prec/nut          ire-gcrf
#
#  references    :
#    vallado       2004, 205-219
#
# [recef, vecef, aecef] = eci2ecefiau06  (reci, veci, aeci, ttt, jdut1, lod, xp, yp, option, ddx, ddy)
# ----------------------------------------------------------------------------

def eci2ecefiau06(reci: np.ndarray, veci: np.ndarray, aeci:np.ndarray,
                  ttt: float, jdut1: float, lod: float, xp: float,
                  yp: float, option: str, ddx: float, ddy:float):
    """this function trsnforms a vector from the mean equator mean equniox frame
    (gcrf), to an earth fixed (itrf) frame.  the results take into account
    the effects of precession, nutation, sidereal time, and polar motion.

    Parameters
    ----------
    reci : ndarray
        position vector eci: km
    veci : ndarray
        velocity vector eci: km/s
    aeci : ndarray
        acceleration vector eci: km/s2
    ttt : float
        julian centuries of date: centuries
    jdut1 : float
        julian date of ut1: days from 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    option : str
        approach to use: 'a' - classical equinox 2000a, 'b' - classical equinox
        2000b, 'c' - cio iau2006
    ddx : float
        eop correction for x: rad
    ddy : float
        eop correction for y: rad

    Returns
    -------
    recef: ndarray
        position vector ecef: km
    vecef: ndarray
        velocity vector ecef: km/s
    aecef: ndarray
        acceleration vector ecef: km/s2
    """
    # ---- ceo based, iau2006
    if option == 'c':
        x, y, s, pnb = obu.iau06xys(ttt, ddx, ddy)
        st = obu.iau06era(jdut1)
    # ---- class equinox based, 2000a
    if option == 'a':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pna(ttt)
        gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')
    # ---- class equinox based, 2000b
    if option == 'b':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pnb(ttt)
        gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')

    pm = smu.polarm(xp, yp, ttt, '01')
    # ---- setup parameters for velocity transformations
    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])

    rtirs = st.T@pnb.T@reci
    recef = pm.T@rtirs

    vtirs = st.T@pnb.T@veci - np.cross(omegaearth, rtirs.T).T
    vecef = pm.T@vtirs
    if sh.iauhelp:
        print("spv")
        print(st.T@pnb.T@veci)
        print(np.cross(omegaearth, rtirs.T))
        print("pmt")
        print(pm.T)
        print("vtirs")
        print(vtirs)

    temp = np.cross(omegaearth, rtirs.T)
    aecef = pm.T@(st.T@pnb.T@aeci - np.cross(omegaearth, temp).T
                  - 2.0*np.cross(omegaearth, vtirs.T).T)


    if sh.iauhelp:
        rcirs = pnb.T @ reci
        vcirs = pnb.T @ veci
        if (option == 'a') or (option == 'b'):
            rmod20 = prec.T @ reci
            vmod20 = prec.T @ veci
            print('eci           IAU-2006 %c     ' % (option))
            print(reci)
            print('MOD           IAU-2006 %c     ' % (option))
            print(rmod20)
            print(' v ', (vmod20))
            print(' a', (aeci))
            print('ERS           IAU-2006 %c   ' % (option))
            print(rcirs)

        if option == 'c':
            print('CIRS          IAU-2006 CIO ', (rcirs))
            print(' v ', (vcirs))
            print('TIRS          IAU-2006 %c   ' % (option))
            print(rtirs)
            print(' v ', (vtirs))

    return recef, vecef, aecef


# ----------------------------------------------------------------------------
#
#                           function ecef2cirsiau06
#
#  this function transforms a vector from the earth fixed (itrf) frame, to
#    the cirs. Sidereal time and polar motion are taken into account.
#
#  author        : david vallado                  719-573-2600   2 may 2020
#
#  revisions
#
#  inputs          description                    range / units
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    option      - which approach to use          a-2000a, b-2000b, c-2000xys
#  outputs       :
#    rcirs       - position vector cirs            km
#    vcirs       - velocity vector cirs            km/s
#    acirs       - acceleration vector cirs        km/s2
#
#  locals        :
#    pm          - transformation matrix for itrf-pef
#    st          - transformation matrix for pef-ire
#
#  coupling      :
#   iau00pm      - rotation for polar motion      itrf-pef
#   iau00era     - rotation for earth rotyation   pef-ire
#   iau00xys     - rotation for prec/nut          ire-gcrf

#
#  references    :
#    vallado       2004, 205-219
#
# [rcirs, vcirs, acirs] = ecef2cirsiau06 (recef, vecef, aecef, ttt, jdut1, lod, xp, yp, option, ddx, ddy)
# ----------------------------------------------------------------------------

def ecef2cirsiau06(recef: np.ndarray, vecef: np.ndarray, aecef:np.ndarray,
                   ttt: float, jdut1: float, lod: float, xp: float, yp: float,
                   option: str, ddx: float, ddy: float):
    """this function transforms a vector from the earth fixed (itrf) frame, to
    the cirs. Sidereal time and polar motion are taken into account.

    Parameters
    ----------
    recef: ndarray
        position vector ecef: km
    vecef: ndarray
        velocity vector ecef: km/s
    aecef: ndarray
        acceleration vector ecef: km/s2
    ttt : float
        julian centuries of date: centuries
    jdut1 : float
        julian date of ut1: days from 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    option : str
        approach to use: 'a' - classical equinox 2000a, 'b' - classical equinox
        2000b, 'c' - cio iau2006
    ddx : float
        eop correction for x: rad
    ddy : float
        eop correction for y: rad

    Returns
    -------
    rcirs : ndarray
        position vector cirs: km
    vcirs : ndarray
        velocity vector cirs: km/s
    acirs : ndarray
        acceleration vector cirs: km/s2
    """

    # ---- cio based, iau2006
    if option == 'c':
        x, y, s, pnb = obu.iau06xys(ttt, ddx, ddy)
        st = obu.iau06era(jdut1)
    # ---- class equinox based, 2000a
    if option == 'a':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pna(ttt)
        gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')
    # ---- class equinox based, 2000b
    if option == 'b':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pnb(ttt)
        gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')

    pm = smu.polarm(xp, yp, ttt, '06')
    # ---- setup parameters for velocity transformations
    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    # ---- perform transformations
    rpef = pm @ recef
    rcirs = st @ rpef
    vpef = pm @ vecef
    vcirs = st @ (vpef + np.cross(omegaearth, rpef.T))
    temp = np.cross(omegaearth, rpef.T)
    acirs = st @ (pm @ aecef + np.cross(omegaearth, temp).T + 2.0
                  * np.cross(omegaearth, vpef.T))
    return rcirs, vcirs, acirs

# ----------------------------------------------------------------------------
#
#                           function ecef2eciiau06
#
#  this function transforms a vector from the earth fixed (itrf) frame, to
#    the eci mean equator mean equinox (gcrf).
#
#  author        : david vallado                  719-573-2600   16 jul 2004
#
#  revisions
#
#  inputs          description                    range / units
#    recef       - position vector earth fixed    km
#    vecef       - velocity vector earth fixed    km/s
#    aecef       - acceleration vector earth fixedkm/s2
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    option      - which approach to use          a-2000a, b-2000b, c-2000xys
#    ddx         - eop correction for x           rad
#    ddy         - eop correction for y           rad
#
#  outputs       :
#    reci        - position vector eci            km
#    veci        - velocity vector eci            km/s
#    aeci        - acceleration vector eci        km/s2
#
#  locals        :
#    pm          - transformation matrix for itrf-pef
#    st          - transformation matrix for pef-ire
#    nut         - transformation matrix for ire-gcrf
#
#  coupling      :
#   iau00pm      - rotation for polar motion      itrf-pef
#   iau00era     - rotation for earth rotyation   pef-ire
#   iau00xys     - rotation for prec/nut          ire-gcrf

#
#  references    :
#    vallado       2004, 205-219
#
# [reci, veci, aeci] = ecef2eciiau06 (recef, vecef, aecef, ttt, jdut1, lod, xp, yp, option, ddx, ddy)
# ----------------------------------------------------------------------------

def ecef2eciiau06(recef: np.ndarray, vecef:np.ndarray, aecef:np.ndarray,
                  ttt: float, jdut1: float, lod: float, xp: float, yp: float,
                  option: str, ddx: float, ddy: float):
    """this function transforms a vector from the earth fixed (itrf) frame, to
    the eci mean equator mean equinox (gcrf).

    Parameters
    ----------
    recef: ndarray
        position vector ecef: km
    vecef: ndarray
        velocity vector ecef: km/s
    aecef: ndarray
        acceleration vector ecef: km/s2
    ttt : float
        julian centuries of date: centuries
    jdut1 : float
        julian date of ut1: days from 4713 bc
    lod : float
        excess length of day: sec
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    option : str
        approach to use: 'a' - classical equinox 2000a, 'b' - classical equinox
        2000b, 'c' - cio iau2006
    ddx : float
        eop correction for x: rad
    ddy : float
        eop correction for y: rad

    Returns
    -------
    reci: ndarray
        position vector eci: km
    veci: ndarray
        velocity vector eci: km/s
    aeci: ndarray
        acceleration vector eci: km/s2
    """

    pnb = np.eye(3)
    st = np.eye(3)

    # ---- ceo based, iau2006
    if option == 'c':
        x, y, s, pnb = obu.iau06xys(ttt, ddx, ddy)
        st = obu.iau06era(jdut1)
    # ---- class equinox based, 2000a
    if option == 'a':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pna(ttt)
        gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')
    # ---- class equinox based, 2000b
    if option == 'b':
        deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
            lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pnb(ttt)
        gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')

    pm = smu.polarm(xp, yp, ttt, '06')
    # ---- setup parameters for velocity transformations
    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    # ---- perform transformations
    rpef = pm @ recef
    reci = pnb @ st @ rpef

    vpef = pm @ vecef
    veci = pnb @ st @ (vpef + np.cross(omegaearth, rpef.T).T)
    temp = np.cross(omegaearth, rpef.T)
    aeci = pnb @ st @ (pm @ aecef + np.cross(omegaearth, temp).T + 2.0 * np.cross(omegaearth, vpef.T).T)
    return reci, veci, aeci


# -----------------------------------------------------------------------------
#
#                           function rad2hms
#
#  this function converts radians to hours, minutes and seconds.  notice
#    the conversion 0.2617 is simply the radian equivalent of 15 degrees.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    hms         - result                         rad
#
#  outputs       :
#    hr          - hours                          0 .. 24
#    min         - minutes                        0 .. 59
#    sec         - seconds                        0.0 .. 59.99
#
#  locals        :
#    temp        - conversion from hours to rad   0.261799
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2001, 200, alg 19 alg 20, ex 3-9
#
# [hr, min, sec] = rad2hms(hms)
# -----------------------------------------------------------------------------

def rad2hms(hms: float):
    """this function converts radians to hours, minutes and seconds. notice
    the conversion 0.2617 is simply the radian equivalent of 15 degrees.

    Parameters
    ----------
    hms : float
        radians to convert

    Returns
    -------
    hr: float
        hours: 0...24
    min: float
        minutes: 0...59
    sec: float
        seconds: 0...59.99
    """

    temp = 15.0 * math.pi/180.0
    temp = hms / temp
    hr = np.fix(temp)
    min = np.fix((temp - hr)*60.0)
    sec = (temp - hr - min/60.0) * 3600.0
    return hr, min, sec


# -----------------------------------------------------------------------------
#
#                           function hms2rad
#
#  this function converts hours, minutes and seconds into radians.  notice
#    the conversion 0.2617 is simply the radian equivalent of 15 degrees.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    hr          - hours                          0 .. 24
#    min         - minutes                        0 .. 59
#    sec         - seconds                        0.0 .. 59.99
#
#  outputs       :
#    hms         - result                         rad
#
#  locals        :
#    temp        - conversion from hours to rad   0.261799
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 204, alg 19 alg 20, ex 3-9
#
# [hms] = hms2rad(hr, min, sec)
# -----------------------------------------------------------------------------

def hms2rad(hr: float, min: float, sec: float):
    """this function converts hours, minutes and seconds into radians. notice
    the conversion 0.2617 is simply the radian equivalent of 15 degrees.

    Parameters
    ----------
    hr: float
        hours: 0...24
    min: float
        minutes: 0...59
    sec: float
        seconds: 0...59.99

    Returns
    -------
    hms: float
        radians
    """

    temp = 15.0 * math.pi/180.0

    hms = (hr + min/60.0 + sec/3600.0)*temp
    return hms

#  -----------------------------------------------------------------------------
#
#                            procedure twoline2rv
#
#  this function converts the two line element set character string data to
#    variables and initializes the sgp4 variables. several intermediate varaibles
#    and quantities are determined. note that the result is a structure so multiple
#    satellites can be processed simultaneously without having to reinitialize. the
#    verification mode is an important option that permits quick checks of any
#    changes to the underlying technical theory. this option works using a
#    modified tle file in which the start, stop, and delta time values are
#    included at the end of the second line of data. this only works with the
#    verification mode. the catalog mode simply propagates from -1440 to 1440 min
#    from epoch and is useful when performing entire catalog runs.
#
# Author:
#   Jeff Beck
#   beckja@alumni.lehigh.edu
#   1.0  aug  6, 2006 - update for paper dav
#   2.0  mar  8, 2007 - misc fixes and manual operation updates
#   2.01 may  9, 2007 - fix for correction to year of 57
#   2.02 oct  8, 2007 - fix for manual jdstart jdstop matlab formats
# original comments from Vallado C++ version:
#   author        : david vallado                  719-573-2600    1 mar 2001
#
#   inputs        :
#   longstr1      - TLE character string
#   longstr2      - TLE character string
#   typerun       - character for mode of SGP4 execution
#                   'c' = catalog mode (propagates at 20 min timesteps from
#                           one day before epoch to one day after)
#                   'v' = verification mode (propagates according to start,
#                           stop, and timestep specified in longstr2)
#                   'm' = manual mode (prompts user for start, stop, and
#                           timestep for propagation)
#                   'u' = user imnput, whatever
#   typeinput     - type of manual input           mfe 'm', epoch 'e', dayofyr 'd'
#
#   outputs       :
#     satrec      - structure containing all the sgp4 satellite information
#
#   coupling      :
#     getgravconst
#     days2mdhms  - conversion of days to month, day, hour, minute, second
#     jday        - convert day month year hour minute second into julian date
#     sgp4init    - initialize the sgp4 variables
#
#   references    :
#     norad spacetrack report #3
#     vallado, crawford, hujsak, kelso  2006
#
# [startmfe, stopmfe, deltamin, satrec] = twoline2rv(longstr1, longstr2, ...
#          typerun, typeinput, opsmode, whichconst)
#  ----------------------------------------------------------------------------

def twoline2rv(longstr1: str, longstr2: str, typerun: str,
               typeinput, opsmode: str, whichconst: str):
    """this function converts the two line element set character string data to
    variables and initializes the sgp4 variables. several intermediate varaibles
    and quantities are determined. note that the result is a structure so multiple
    satellites can be processed simultaneously without having to reinitialize. the
    verification mode is an important option that permits quick checks of any
    changes to the underlying technical theory. this option works using a
    modified tle file in which the start, stop, and delta time values are
    included at the end of the second line of data. this only works with the
    verification mode. the catalog mode simply propagates from -1440 to 1440 min
    from epoch and is useful when performing entire catalog runs.

    Parameters
    ----------
    longstr1 : str
        TLE character string
    longstr2 : str
        TLE character string
    typerun : character for mode of SGP4 execution
        'c' = catalog mode (propagates at 20 min timesteps from
            one day before epoch to one day after)
        'v' = verification mode (propagates according to start,
            stop, and timestep specified in longstr2)
        'u' = user input, whatever
    typeinput : type of manual input
        used for commented out 'm' typerun - currently unused
    opsmode : str
        satrec opsmode
    whichconst : str
        which set of constants to use: 721, 72, 84
    """

    # sgp4fix no longer needed, put in satrec
    # global tumin radiusearthkm xke j2 j3 j4 j3oj2

    # Needs a constant in space_constants?
    xpdotp   =  1440.0 / (2.0*math.pi)   # 229.1831180523293  # [rev/day]/[rad/min]

    year   = 0
    satrec = {}
    satrec['error'] = 0


    longstr1 = list(longstr1)
    longstr2 = list(longstr2)

    # set the implied decimal points since doing a formated read
    # fixes for bad input data values (missing, ...)
    for j in range(10, 16):
        if (longstr1[j] == ' '):
            longstr1[j] = '_'



    if (longstr1[44] != ' '):
        longstr1[43] = longstr1[44]

    longstr1[44] = '.'

    if (longstr1[7] == ' '):
        longstr1[7] = 'U'


    if (longstr1[9] == ' '):
        longstr1[9] = '.'


    for j in range(45, 50):
        if (longstr1[j] == ' '):
            longstr1[j] = '0'


    if (longstr1[51] == ' '):
        longstr1[51] = '0'

    if (longstr1[53] != ' '):
        longstr1[52] = longstr1[53]

    longstr1[53] = '.'

    longstr2[25] = '.'

    for j in range(26, 33):
        if (longstr2[j] == ' '):
            longstr2[j] = '0'



    if (longstr1[62] == ' '):
        longstr1[62] = '0'


    if ((len(longstr1) < 68) or (longstr1[67] == ' ')):
        longstr1[67] = '0'

    longstr1 = "".join(longstr1)
    longstr2 = "".join(longstr2)

    # parse first line
    carnumb =float(longstr1[0])
    satrec['satnum'] = float(longstr1[2:7])
    satrec['classification'] = longstr1[7]
    satrec['intldesg'] = longstr1[9:17]
    satrec['epochyr'] = float(longstr1[18:20])
    satrec['epochdays'] = float(longstr1[20:32])
    satrec['ndot'] = float(longstr1[33:43])
    satrec['nddot'] = float(longstr1[43:50])
    nexp = float(longstr1[50:52])
    satrec['bstar'] = float(longstr1[52:59])
    ibexp = float(longstr1[59:61])
    numb = float(longstr1[62])
    satrec['elnum'] = float(longstr1[64:68])

    # parse second line
    cardnumb = float(longstr2[0])
    satrec['satnum'] = float(longstr2[2:7])
    satrec['inclo'] = float(longstr2[7:16])
    satrec['nodeo'] = float(longstr2[16:25])
    satrec['ecco'] = float(longstr2[25:33])
    satrec['argpo'] = float(longstr2[33:42])
    satrec['mo'] = float(longstr2[42:51])
    satrec['no_kozai'] = float(longstr2[51:63])
    satrec['revnum'] = float(longstr2[63:68])

    if (typerun == 'v'):
        startmfe = float(longstr2[69:81])
        stopmfe  = float(longstr2[82:96])
        deltamin = float(longstr2[96:105])
    # perform complete catalog evaluation
    if (typerun == 'c'):
        startmfe =  -1440.0
        stopmfe  =  1440.0
        deltamin = 20.0
    # user input
    if (typerun == 'u'):
        startmfe =  0.0
        stopmfe  =  14400.0
        deltamin = 1440.0

    # ---- find no, ndot, nddot ----
    satrec['no_kozai'] = satrec['no_kozai'] / xpdotp #//* rad/min
    satrec['nddot'] = satrec['nddot'] * 10.0**nexp
    # note the implied decimal is set when adjusting longstr1 above
    satrec['bstar'] = satrec['bstar'] * 10.0**ibexp

    # ---- convert to sgp4 units ----
    #    satrec['a']    = (satrec['no']*tumin)^(-2/3)                # [er]
    satrec['ndot'] = satrec['ndot'] / (xpdotp * 1440.0)          # [rad/min^2]
    satrec['nddot'] = satrec['nddot'] / (xpdotp * 1440.0 * 1440)     # [rad/min^3]

    # ---- find standard orbital elements ----
    satrec['inclo'] = satrec['inclo'] * deg2rad
    satrec['nodeo'] = satrec['nodeo'] * deg2rad
    satrec['argpo'] = satrec['argpo'] * deg2rad
    satrec['mo'] = satrec['mo'] * deg2rad

    # sgp4fix not needed here
    #    satrec['alta'] = satrec['a']*(1.0 + satrec['ecco']) - 1.0
    #    satrec['altp'] = satrec['a']*(1.0 - satrec['ecco']) - 1.0

    # ----------------------------------------------------------------
    # find sgp4epoch time of element set
    # remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
    # and minutes from the epoch (time)
    # --------------------------------------------------------------

    # ------------- temp fix for years from 1957-2056 ----------------
    # ------ correct fix will occur when year is 4-digit in 2le ------
    if (satrec['epochyr'] < 57):
        year = satrec['epochyr'] + 2000
    else:
        year = satrec['epochyr'] + 1900


    mon, day, hr, minute, sec = stu.days2mdh(year, satrec['epochdays'])
    satrec['jdsatepoch'], satrec['jdsatepochf'] = stu.jday(year, mon, day,
                                                           hr, minute, sec)

    # default values
    #startmfe = 0.0
    #stopmfe  = 1440.0
    #deltamin = 1.0

    # input start stop times manually
    # if ((typerun != 'v') and (typerun != 'c')  and (typerun != 'u')):
    #     # ------------- enter start/stop ymd hms values --------------------
    #     if (typeinput == 'e'):
    #         startyear = input('input start year')
    #         startmon  = input('input start mon')
    #         startday  = input('input start day')
    #         starthr   = input('input start hr')
    #         startmin  = input('input start min')
    #         startsec  = input('input start sec')
    #         jdstart, jdstartf = stu.jday(startyear, startmon, startday,
    #                                      starthr, startmin, startsec)

    #         stopyear = input('input stop year')
    #         stopmon  = input('input stop mon')
    #         stopday  = input('input stop day')
    #         stophr   = input('input stop hr')
    #         stopmin  = input('input stop min')
    #         stopsec  = input('input stop sec')
    #         jdstop, jdstopf = stu.jday(stopyear, stopmon, stopday, stophr,
    #                                    stopmin, stopsec)

    #         startmfe = (jdstart + jdstartf - satrec['jdsatepoch']
    #                     - satrec['jdsatepochf']) * 1440.0
    #         stopmfe  = (jdstop + jdstopf - satrec['jdsatepoch']
    #                     - satrec['jdsatepochf']) * 1440.0
    #         deltamin = input('input time step in minutes ')

    #     # -------- enter start/stop year and days of year values -----------
    #     if (typeinput == 'd'):
    #         startyear    = input('input start year')
    #         startdayofyr = input('input start dayofyr')
    #         stopyear     = input('input stop year')
    #         stopdayofyr  = input('input stop dayofyr')

    #         [mon,day,hr,minute,sec] = days2mdh ( startyear,startdayofyr)
    #         [jdstart,jdstartf] = jday( startyear,mon,day,hr,minute,sec)
    #         [mon,day,hr,minute,sec] = days2mdh ( stopyear,stopdayofyr)
    #         [jdstop, jdstopf] = jday( stopyear,mon,day,hr,minute,sec)

    #         startmfe = (jdstart + jdstartf - satrec['jdsatepoch'] - satrec['jdsatepochf']) * 1440.0
    #         stopmfe  = (jdstop + jdstopf - satrec['jdsatepoch'] - satrec['jdsatepochf']) * 1440.0
    #         deltamin = input('input time step in minutes ')

    #     # ------------------ enter start/stop mfe values -------------------
    #     if (typeinput == 'm'):
    #         startmfe = input('input start mfe: ')
    #         stopmfe  = input('input stop mfe: ')
    #         deltamin = input('input time step in minutes: ')




    # ------------- initialize the orbit at sgp4epoch --------------
    sgp4epoch = satrec['jdsatepoch'] + satrec['jdsatepochf'] - 2433281.5 # days since 0 Jan 1950
    satrec = obu.sgp4init(whichconst, opsmode, satrec,
                        satrec['jdsatepoch'] + satrec['jdsatepochf'] - 2433281.5,
                        satrec['bstar'], satrec['ndot'], satrec['nddot'],
                        satrec['ecco'], satrec['argpo'], satrec['inclo'],
                        satrec['mo'], satrec['no_kozai'],satrec['nodeo'])
    return startmfe, stopmfe, deltamin, satrec


def hilleqcm2eci(rtgt: np.ndarray, vtgt: np.ndarray, x: float, y: float,
                  z: float, dx: float, dy: float, dz: float):
    """this function takes the position and velocity of a target orbit, and the
    relative position and velocity of the interceptor to get the inertial
    position and velocity of the interceptor

    Parameters
    ----------
    rtgt : ndarray
        target position: km
    vtgt : ndarray
        target velocity: km/s
    x, y, z : float
        relative position of interceptor: km
    dx, dy, dz : float
        relative velocity of interceptor: km/s

    Returns
    -------
    rinteci: ndarray
        interceptor position: km
    vinteci: ndarray
        interceptor velocity: km/s
    """

    rrsw1, vrsw1, transmat1 = rv2rsw(rtgt, vtgt)

    rmag1 = smu.mag(rrsw1)
    vmag1 = smu.mag(vrsw1)

    # TODO: Add ebar to rv2coe outputs instead of computing it again here -zeg
    ptgt, atgt, ecc, _, _, _, _, _, _, _, _ = rv2coe(rrsw1, vrsw1)
    ebar = ((vmag1**2 - mu / rmag1) * rrsw1 - (rrsw1 @ vrsw1) * vrsw1) / mu
    eunit = ebar / ecc
    lambdap = math.atan(eunit[1] / eunit[0])
    nu1 = -lambdap


    arclength = y
    ea1, _ = smu.newtonnu(ecc, nu1)
    if abs(arclength) > 0.001:
        DE = arclength / atgt
        deltaea = inverselliptic2(DE, ecc **2)
        F1, E1, _ = elliptic12(ea1, ecc **2)
        ea2e = ea1 + deltaea
        i = 1
        arclength1a = arclength + 10
        while (i < 10) and (abs(arclength1a - arclength) > 0.001):
            F2, E2, _ = elliptic12(ea2e, ecc**2)
            arclength1a = atgt * (E2 - E1)
            corr = arclength / (ea2e - ea1)
            ea2e = ea2e - (arclength1a - arclength) / corr
            i = i + 1

        _, nu2 = smu.newtone(ecc, ea2e)

    cosnu2 = math.cos(nu2)
    sinnu2 = math.sin(nu2)

    rpqw2 = np.array([(ptgt * cosnu2) / (1 + ecc * cosnu2),
                      (ptgt * sinnu2) / (1 + ecc * cosnu2),
                      0])
    vpqw2 = np.array([(-math.sqrt(mu / ptgt) * sinnu2),
                      math.sqrt(mu / ptgt) * (ecc + cosnu2),
                      0])

    rrsw2, vrsw2, transmat2 = rv2rsw(rpqw2, vpqw2)
    dphi = math.asin(z / smu.mag(rrsw2))
    dlambda = nu2 - nu1

    sindphi = math.sin(dphi)
    cosdphi = math.cos(dphi)
    sindlambda = math.sin(dlambda)
    cosdlambda = math.cos(dlambda)

    rintrsw1unit = np.array([cosdphi * cosdlambda,
                             cosdphi * sindlambda,
                             sindphi])

    #a few dphi and dlambda mixed up - mjc
    rsw2sez = np.array([[sindphi * cosdlambda, sindphi * sindlambda, -cosdphi],
                          [-sindlambda, cosdlambda, 0],
                          [cosdphi * cosdlambda, cosdphi * sindlambda, sindphi]])
    rintsezunit = rsw2sez @ rintrsw1unit

    rintsezZ = x + rrsw2[0]
    rscale = rintsezZ / rintsezunit[2]

    rintrsw1 = rscale * rintrsw1unit

    magrtgt2 = smu.mag(rpqw2)
    vintsez = [-(dz/magrtgt2) * rscale,
               (dy + vrsw1[1]) * rscale * (cosdphi/magrtgt2),
               dx + vrsw2[0]]

    vintrsw1 = rsw2sez.T @ vintsez

    rinteci = transmat1.T @ rintrsw1
    vinteci = transmat1.T @ vintrsw1

    return rinteci, vinteci

def eci2hilleqcm(rtgt: np.ndarray, vtgt: np.ndarray, rint: np.ndarray,
                 vint: np.ndarray):
    """this function takes the eci position and velocity vectors of
    a target and interceptor, and returns the relative eqcm position of the
    interceptor

    Parameters
    ----------
    rtgt : ndarray
        position vector of the target: km
    vtgt : ndarray
        velocity vector of the target: km/s
    rint : ndarray
        positionv vector of interceptor: km
    vint : ndarray
        velocity vector of interceptor: km/s

    Returns
    -------
    x, y, z: float
        relative position of interceptor: km
    dx, dy, dz: float
        relative velocity of interceptor: km/s
    """
    rtgtrsw1, vtgtrsw1, transmat1 = rv2rsw(rtgt, vtgt)
    rintrsw1 = transmat1 @ rint
    vintrsw1 = transmat1 @ vint

    # TODO: Add ebar to rv2coe outputs instead of computing it again here -zeg
    ptgt, atgt, ecc, _, _, _, _, _, _, _, _ = rv2coe(rtgtrsw1, vtgtrsw1)
    rtgtmag1 = smu.mag(rtgtrsw1)
    vtgtmag1 = smu.mag(vtgtrsw1)
    rintmag1 = smu.mag(rintrsw1)

    dphi = math.asin(rintrsw1[2] / rintmag1)
    dlambda = math.atan(rintrsw1[1] / rintrsw1[0])

    ebar = ((vtgtmag1**2 - mu / rtgtmag1) * rtgtrsw1
            - (rtgtrsw1 @ vtgtrsw1) * vtgtrsw1) / mu
    eunit = smu.unit(ebar)
    lambdap = math.atan(eunit[1]/ eunit[0])
    nu1 = -lambdap
    nu2 = dlambda - lambdap

    sinnu2, cosnu2, sindphi, cosdphi, sindlambda, cosdlambda = \
        smu.getsincos(nu2, dphi, dlambda)
    rtgtpqw2 = np.array([(ptgt * cosnu2) / (1 + ecc * cosnu2),
                         (ptgt * sinnu2) / (1 + ecc * cosnu2),
                         0])

    vtgtpqw2 = np.array([-math.sqrt(mu / ptgt) * sinnu2,
                         math.sqrt(mu/ptgt) * (ecc + cosnu2),
                         0])
    rtgtrsw2, vtgtrsw2, transmat2 = rv2rsw(rtgtpqw2, vtgtpqw2)

    rsw2sez = np.array([[sindphi * cosdlambda, sindphi * sindlambda, -cosdphi],
                        [-sindlambda, cosdlambda, 0],
                        [cosdphi * cosdlambda, cosdphi * sindlambda, sindphi]])
    rintsez = rsw2sez @ rintrsw1
    vintsez = rsw2sez @ vintrsw1


    # arc1 = IEISK((E1, E2), ecc**2)
    btgt = math.sqrt(atgt * ptgt)
    ea1,_ = smu.newtonnu(ecc,nu1)
    ea2,_ = smu.newtonnu(ecc,nu2)

    if np.abs(ea2 - ea1) > np.pi:
        if ea1 < 0.0:
            ea1 = 2.0 * np.pi + ea1
        else:
            ea1 = 2.0 * np.pi - ea1

    arc1 = arclength_ellipse(atgt,btgt,ea1,ea2)[0]
    #print("testing arclength func: %f", arc1)

    # F1, E1, _ = elliptic12(ea1, ecc ** 2)
    # F2, E2, _ = elliptic12(ea2, ecc ** 2)

    # fixit = 0.0
    # if E2 - E1 < 0.0:
    #     fixit = np.pi

    # arc1 = atgt * (E2 - E1 + fixit)
    # print("testing elliptic12 func: %f", arc1)

    rtgtmag2 = smu.mag(rtgtpqw2)
    rinteqcm = np.array([rintsez[2] - rtgtrsw2[0],
                         arc1,
                         sindphi * rtgtmag2])

    dotlambda = vintsez[1] / (rintmag1 * cosdphi)
    dotphi = -vintsez[0] / rintmag1
    vinteqcm = np.array([vintsez[2] - vtgtrsw2[0],
                           dotlambda * rtgtmag2 - abs(vtgtrsw1[1]),
                           dotphi * rtgtmag2])

    return rinteqcm[0], rinteqcm[1], rinteqcm[2], \
        vinteqcm[0], vinteqcm[1], vinteqcm[2]


# ----------------------------------------------------------------------------
#
#                           function fk4
#
#    This function converts vectors from b1950 to j2000 epochs.
#    be aware that this process is not exact. there are different secular rates
#    for each system, and there are differences in the central location. the
#    matrices are multiplied directly for speed.
#
#  author        : david vallado                  719-573-2600   21 jun 2002
#
#  revisions     : michael courville                             2 jul 2024
#                -
#
#  inputs          description                              range / units
#    rb1950      - b1950 eci position vector (optional)      er, km, etc
#    vb1950      - b1950 eci velocity vector (optional)      er/s, km/s, etc
#                  (Note: If only trying to find vj2000,
#                  input an empty numpy list for rb1950)
#    option      - systems tool kit, original,              'stk' (default)
#                  or 6 dimension approach                  'org', '6d'
#                  (Note: if only trying to only return
#                  trans martix, input an empty numpy list
#                  for rb1950 and vb1950)
#
#  outputs       :
#    rj2000      - j2000 eci position vector                er, km, etc
#                - 3x3 -> 3x1 [rx, ry, rz]
#                - 6x6 -> 6x1 [rx, ry, rz, rxdot, rydot, rzdot]
#    vj2000      - j2000 eci velocity vector                er/s, km/s, etc
#                - 3x3 -> 3x1 [vx, vy, vz]
#                - 6x6 -> 6x1 [vx, vy, vz, vxdot, vydot, vzdot]
#    fk4m        - conversion matrix b1950 to j2000
#
#  locals        :
#
#  coupling      :
#
#  references    :
#    vallado       2001, 227-228
#
# rj2000 = fk4m * rb1950;
# ----------------------------------------------------------------------------
def fk4(rb1950: np.ndarray = np.array([]), vb1950: np.ndarray = np.array([]),
        option: str ='org'):
    """This function converts vectors from b1950 to j2000 epochs.
    be aware that this process is not exact. there are different secular rates
    for each system, and there are differences in the central location. the
    matrices are multiplied directly for speed.

    Parameters
    ----------
    rb1950 : np.ndarray, optional
        b1950 eci position vector, by default np.array([]): er, km, etc
    vb1950 : np.ndarray, optional
        b1950 eci velocity vector, by default np.array([]): er/s, km/s, etc
    option : str, optional
        transformation approach (systems tool kit, original, or 6 dimension), \
        by default 'org': 'org','stk','6d'

    Note: if only trying to only return transformation martix, input an empty \
        numpy array for rb1950 and vb1950

    Returns
    -------
    rj2000: np.ndarray
        j2000 eci position vector: er, km, etc
    vj2000: np.ndarray
        j2000 eci velocity vector: er/s, km/s, etc
    fk4m: np.ndarray
        conversion matrix for provided option
    """

    # Default Empty Outputs
    rj2000 = np.array([])
    vj2000 = np.array([])
    fk4m = np.array([])

    # Input Error Check:
    if not isinstance(rb1950, np.ndarray) or not isinstance(vb1950, np.ndarray):
        print('\nError: Position and velocity vectors must be numpy arrays\n')
        print('Returning empty arrays')
        return rj2000, vj2000, fk4m
    elif not isinstance(option, str):
        print('\nError: Invalid option\n')
        print('Returning empty arrays')
        return rj2000, vj2000, fk4m
    elif option not in ['org','stk','6d']:
        print('\nError: Invalid option\n')
        print('Returning empty arrays')
        return rj2000, vj2000, fk4m


    # Transformation martrix options:
     # in book (original)
    # 1950 - 2000 (Vallado 4th edition, pages 234-235)
    if option == 'org':
        fk4m = np.array([[0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
                            [0.0111814832391717, 0.9999374848933135, -0.0000271625947142],
                            [0.0048590037723143, -0.0000271702937440, 0.9999881946043742]])
    # stk approach (the difference is about 5 m between the various approaches)
    # New way is formed by multiplying the matrices on pages
    # 173 and 174 and adding in the correction to equinox given
    # on page 168 of the supplement to the astronomical almanac
    # 1950 - 2000
    elif option == 'stk':
        fk4m = np.array([[0.999925678612394, -0.011181874556714, -0.004858284812600],
                            [0.011181874524964, 0.999937480517880, -0.000027169816135],
                            [0.004858284884778, -0.000027156932874, 0.999988198095508]])
    # from Exp supp to Ast Almanac pg 185 6x6
    # 1950 - 2000
    elif option == '6d' or (rb1950.size == 6 or vb1950.size == 6):
        fk4m = np.array([[0.9999256782, -0.0111820611, -0.0048579477,
                            0.00000242395018, -0.00000002710663, -0.00000001177656],
                            [0.0111820610, 0.9999374784, -0.0000271765,
                            0.00000002710663, 0.00000242397878, -0.00000000006587],
                            [0.0048579479, -0.0000271474, 0.9999881997,
                            0.00000001177656, -0.00000000006582, 0.00000242410173],
                            [-0.000551, -0.238565, 0.435739,
                            0.99994704, -0.01118251, -0.00485767],
                            [0.238514, -0.002667, -0.008541,
                            0.01118251, 0.99995883, -0.00002718],
                            [-0.435623, 0.012254, 0.002117,
                            0.00485767, -0.00002714, 1.00000956]])

    # Matrix size verification:
    if rb1950.size == 0 and vb1950.size == 0:
        print('\nWarning: Position and velocity input are empty\n')
        print('Returning empty position/velocity arrays')
        print('Returning b1950 to j2000 transformation matrix for', option, 'option')
        return rj2000, vj2000, fk4m
    # Position matrix is size 0 if only trying to find velocity
    elif (rb1950.size != 3 and rb1950.size != 0) and (option == 'stk' or option == 'org'):
        print('\nError: Position matrix is not size 3\n')
        print('Returning empty position/velocity arrays')
        print('Returning b1950 to j2000 transformation matrix for', option, 'option')
        return rj2000, vj2000, fk4m
    # Velocity matrix is size 0 if only trying to find velocity
    elif (vb1950.size != 3 and vb1950.size != 0) and (option == 'stk' or option == 'org'):
        print('\nError: Velocity matrix is not size 3\n')
        print('Returning empty position/velocity arrays')
        print('Returning b1950 to j2000 transformation matrix for', option, 'option')
        return rj2000, vj2000, fk4m
    # Position matrix is size 0 if only trying to find velocity
    elif (rb1950.size != 6 and rb1950.size != 0) and option == '6d':
        print('\nError: Position matrix is not size 6\n')
        print('Returning empty position/velocity arrays')
        print('Returning b1950 to j2000 transformation matrix for', option, 'option')
        return rj2000, vj2000, fk4m
    # Velocity matrix is size 0 if only trying to find velocity
    elif (vb1950.size != 6 and vb1950.size != 0) and option == '6d':
        print('\nError: Velocity matrix is not size 6\n')
        print('Returning empty position/velocity arrays')
        print('Returning b1950 to j2000 transformation matrix for', option, 'option')
        return rj2000, vj2000, fk4m

    # Final transformation:
    else:
        if rb1950.size != 0:
            rj2000 = rb1950 @ fk4m
        if vb1950.size != 0:
            vj2000 = vb1950 @ fk4m

    return rj2000, vj2000, fk4m


# ----------------------------------------------------------------------------
#
#                           function fk4i
#
#    This function converts vectors from j2000 to b1950 epochs (i = inverse)
#    be aware that this process is not exact. there are different secular rates
#    for each system, and there are differences in the central location. the
#    matrices are multiplied directly for speed.
#
#  author        : david vallado                  719-573-2600   21 jun 2002
#
#  revisions     : michael courville                             2 jul 2024
#                -
#
#  inputs          description                              range / units
#    rj2000i      - j2000 eci position vector (optional)      er, km, etc
#    vj2000i      - j2000 eci velocity vector (optional)      er/s, km/s, etc
#                  (Note: If only trying to find vb1950i,
#                  input an empty numpy list for rj2000i)
#    option      - systems tool kit, original,              'stk' (default)
#                  or 6 dimension approach                  'org', '6d'
#                  (Note: if only trying to return trans
#                  matrix, input an empty numpy list for
#                  rj2000i and vj2000i)
#
#  outputs       :
#    rb1950i      - b1950 eci position vector                er, km, etc
#                - 3x3 -> 3x1 [rx, ry, rz]
#                - 6x6 -> 6x1 [rx, ry, rz, rxdot, rydot, rzdot]
#    vb1950i      - b1950 eci velocity vector                er/s, km/s, etc
#                - 3x3 -> 3x1 [vx, vy, vz]
#                - 6x6 -> 6x1 [vx, vy, vz, vxdot, vydot, vzdot]
#    fk4mi       - conversion matrix j2000 to b1950
#
#  locals        :
#
#  coupling      :
#
#  references    :
#    vallado       2001, 227-228
#
# rb1950i = fk4mi * rj2000i;
# ----------------------------------------------------------------------------
def fk4i(rj2000i: np.ndarray = np.array([]), vj2000i: np.ndarray = np.array([]),
         option: str ='org'):
    """This function converts vectors from j2000 to b1950 epochs (i = inverse)
    be aware that this process is not exact. there are different secular rates
    for each system, and there are differences in the central location. the
    matrices are multiplied directly for speed.

    Parameters
    ----------
    rj2000i : np.ndarray, optional
        j2000 eci position vector, by default np.array([]): er, km, etc
    vj2000i : np.ndarray, optional
        j2000 eci velocity vector, by default np.array([]): er/s, km/s, etc
    option : str, optional
        transformation approach (systems tool kit, original, or 6 dimension), \
        by default 'org': 'org','stk','6d'

    Note: if only trying to return transformation martix, input an empty \
        numpy array for rj2000i and vj2000i

    Returns
    -------
    rb1950i: np.ndarray
        b1950 eci position vector: er, km, etc
    vb1950i: np.ndarray
        b1950 eci velocity vector: er/s, km/s, etc
    fk4mi: np.ndarray
        conversion matrix for provided option
    """
    # Default Empty Outputs
    rb1950i = np.array([])
    vb1950i = np.array([])
    fk4mi = np.array([])

    # Input Error Check:
    if not isinstance(rj2000i, np.ndarray) or not isinstance(vj2000i, np.ndarray):
        print('\nError: Position and velocity vectors must be numpy arrays\n')
        print('Returning empty arrays')
        return rb1950i, vb1950i, fk4mi
    elif not isinstance(option, str):
        print('\nError: Invalid option\n')
        print('Returning empty arrays')
        return rb1950i, vb1950i, fk4mi
    elif option not in ['org','stk','6d']:
        print('\nError: Invalid option\n')
        print('Returning empty arrays')
        return rb1950i, vb1950i, fk4mi

    # Transformation martrix options:
    # in book (original) (Vallado 4th edition, pages 234-235)
    # 1950 - 2000
    if option == 'org':
        fk4mi = np.array([[9.99925679e-01, 1.11814832e-02, 4.85900377e-03],
                        [-1.11814832e-02, 9.99937485e-01, -2.71702937e-05],
                        [-4.85900382e-03, -2.71625947e-05, 9.99988195e-01]])
    # stk approach (the difference is about 5 m between the various approaches)
    # New way is formed by multiplying the matrices on pages
    # 173 and 174 and adding in the correction to equinox given
    # on page 168 of the supplement to the astronomical almanac
    # 1950 - 2000
    elif option == 'stk':
        fk4mi = np.array([[9.99925679e-01, 1.11818745e-02, 4.85828488e-03],
                        [-1.11818746e-02,  9.99937481e-01, -2.71569329e-05],
                        [-4.85828481e-03, -2.71698161e-05,  9.99988198e-01]])
    # from Exp supp to Ast Almanac pg 185 6x6
    # 1950 - 2000
    elif option == '6d' or (rj2000i.size == 6 or vj2000i.size == 6):
        fk4m = np.array([[0.9999256782, -0.0111820611, -0.0048579477,
                        0.00000242395018, -0.00000002710663, -0.00000001177656],
                        [0.0111820610, 0.9999374784, -0.0000271765,
                        0.00000002710663, 0.00000242397878, -0.00000000006587],
                        [0.0048579479, -0.0000271474, 0.9999881997,
                        0.00000001177656, -0.00000000006582, 0.00000242410173],
                        [-0.000551, -0.238565, 0.435739,
                        0.99994704, -0.01118251, -0.00485767],
                        [0.238514, -0.002667, -0.008541,
                        0.01118251, 0.99995883, -0.00002718],
                        [-0.435623, 0.012254, 0.002117,
                        0.00485767, -0.00002714, 1.00000956]])
        fk4mi = np.linalg.inv(fk4m)

    # Matrix size verification:
    if rj2000i.size == 0 and vj2000i.size == 0:
        print('\nWarning: Position and velocity input are empty\n')
        print('Returning empty position/velocity arrays')
        print (' Returning j2000 to b1950 transformation matrix for', option, 'option')
        return rb1950i, vb1950i, fk4mi
    # Position matrix is size 0 if only trying to find velocity
    elif (rj2000i.size != 3 and rj2000i.size != 0) and (option == 'stk' or option == 'org'):
        print('\nError: Position matrix is not size 3\n')
        print('Returning empty position/velocity arrays')
        print (' Returning j2000 to b1950 transformation matrix for', option, 'option')
        return rb1950i, vb1950i, fk4mi
    # Velocity matrix is size 0 if only trying to find velocity
    elif (vj2000i.size != 3 and vj2000i.size != 0) and (option == 'stk' or option == 'org'):
        print('\nError: Velocity matrix is not size 3\n')
        print('Returning empty position/velocity arrays')
        print (' Returning j2000 to b1950 transformation matrix for', option, 'option')
        return rb1950i, vb1950i, fk4mi
    # Position matrix is size 0 if only trying to find velocity
    elif (rj2000i.size != 6 and rj2000i.size != 0) and option == '6d':
        print('\nError: Position matrix is not size 6\n')
        print('Returning empty position/velocity arrays')
        print (' Returning j2000 to b1950 transformation matrix for', option, 'option')
        return rb1950i, vb1950i, fk4mi
    # Velocity matrix is size 0 if only trying to find velocity
    elif (vj2000i.size != 6 and vj2000i.size != 0) and option == '6d':
        print('\nError: Velocity matrix is not size 6\n')
        print('Returning empty position/velocity arrays')
        print (' Returning j2000 to b1950 transformation matrix for', option, 'option')
        return rb1950i, vb1950i, fk4mi

    # Final transformation:
    else:
        if rj2000i.size != 0:
            rb1950i = rj2000i @ fk4mi
        if vj2000i.size != 0:
            vb1950i = vj2000i @ fk4mi

    return rb1950i, vb1950i, fk4mi

if __name__ == '__main__':

    print('hill to ecqm and ecqm to hill conversion test\n')
    x = 500
    y = 10
    z = 250
    xd = 200
    yd = 0.01
    zd = 20

    rtgteci = np.array([-605.7904308, -5870.230407, 3493.052004])
    vtgteci = np.array([-1.568251615, -3.702348353, -6.479484915])

    rintecix, vintecix = hilleqcm2eci(rtgteci, vtgteci, x, y, z, xd, yd, zd)
    print(rintecix)
    print(vintecix)
    x, y, z, dx, dy, dz = eci2hilleqcm(rtgteci, vtgteci, rintecix, vintecix)
    print(x, y, z)
    print(dx, dy, dz)

    # print('\nConversion Function fk4 and fk4i Testing (Explanatory Supplement Astronomical Almanac, 1992)\n')
    # print('Original fk4 Test (Vallado)')
    # r1, v1, m1 = fk4(np.array([1,2,3]), np.array([4,5,6]), 'org')
    # print('r1')
    # print(r1)
    # print('v1')
    # print (v1)
    # print('m1')
    # print(m1)

    # r2, v2, m2 = fk4i(r1, v1, 'org')
    # print('r2')
    # print(r2)
    # print('v2')
    # print (v2)
    # print('m2')
    # print(m2)

    # print('STK fk4 Test')
    # r1, v1, m1 = fk4(np.array([1,2,3]), np.array([4,5,6]),'stk')
    # print('r1')
    # print(r1)
    # print('v1')
    # print (v1)
    # print('m1')
    # print(m1)

    # r2, v2, m2 = fk4i(r1, v1,'stk')
    # print('r2')
    # print(r2)
    # print('v2')
    # print (v2)
    # print('m2')
    # print(m2)

    # print('6 dimensional fk4 Test (Explanatory Supplement Astronomical \
    #       Almanac, 1992)')
    # r1, v1, m1 = fk4(np.array([1,2,3,4,5,6]), np.array([7,8,9,10,11,12]), '6d')
    # print('r1')
    # print(r1)
    # print('v1')
    # print (v1)
    # print('m1')
    # print(m1)

    # r2, v2, m2 = fk4i(r1, v1, '6d')
    # print('r2')
    # print(r2)
    # print('v2')
    # print (v2)
    # print('m2')
    # print(m2)


