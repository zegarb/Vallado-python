import math
import numpy as np
from numpy.polynomial import Polynomial
from pprint import pprint as pp
from space_constants import *
import orbit_utils
import spacemath_utils as smu
import spacetime_utils as stu
import orbit_utils as obu
from space_constants import sethelp as sh



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
#  inputs          description                    range / units
#    apin        - geomagnetic planetary amplitude
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
#  inputs          description                    range / units
#    kpin        - geomagnetic planetary index
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
#    r           - position vector                km
#    v           - velocity vector                km/s
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
    r: vector array
        position vector
    v: vector array
        velocity vector
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
        r, v = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper)
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
        r = X * fvec + Y * gvec
        v = XD * fvec + YD * gvec
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

    return r, v




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
#    constastro
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
#  inputs          description                    range / units
#    cartcov     - 6x6 cartesian covariance matrix
#    cartstate   - 6x1 cartesian orbit state      (x y z vx vy vz)
#    anom        - anomaly                        'latlon', 'radec'
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    terms       - number of terms for ast calculation 0, 2
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
#    constastro
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
#    constastro
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
#    mum         - gravitational paramater        m**3/s**2 NOTE Meters!
#
#  coupling      :
#    constastro
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
#    terms       - number of terms for ast calculation 0, 2
#
#  outputs       :
#    rho         - satellite range from site      km
#    rtasc       - topocentric right ascension    0.0 to 2pi rad
#    decl        - topocentric declination        -pi/2 to pi/2 rad
#    drho        - range rate                     km/s
#    daz         - xxazimuth rate                   rad / s
#    del         - xxelevation rate                 rad / s
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

def rv2tradc(reci=None, veci=None, latgd=None, lon=None,
             alt=None, ttt=None, jdut1=None, lod=None, xp=None,
             yp=None, terms=None, ddpsi=None, ddeps=None):
    # --------------------- implementation ------------------------
# ----------------- get site vector in ecef -------------------
    rsecef, vsecef = obu.site(latgd, lon, alt)
    #rs
#vs
# -------------------- convert ecef to eci --------------------
    a = np.zeros(3)
    rseci, vseci, aeci = ecef2eci(rsecef, vsecef, a, ttt,
                                jdut1, lod, xp, yp, 2, ddpsi, ddeps)
    #rseci
#vseci

    #rseci = rs
#vseci = vs
#[recef, vecef, aecef] = eci2ecef(reci, veci, aeci, ttt, jdut1, lod, xp, yp, 2, 0, 0)
#reci = recef
#veci = vecef

    # ------- find eci range vector from site to satellite -------
    rhoeci = reci - rseci
    drhoeci = veci - vseci
    rho = smu.mag(rhoeci)
    # ------------- calculate azimuth and elevation ---------------
    temp = np.sqrt(rhoeci[0] * rhoeci[0] + rhoeci[1] * rhoeci[1])
    if (temp < small):
        trtasc = math.atan2(drhoeci[1], drhoeci[0])
    else:
        trtasc = math.atan2(rhoeci[1], rhoeci[0])

    if ((temp < small)):
        tdecl = np.sign(rhoeci[2]) * halfpi
    else:
        magrhoeci = smu.mag(rhoeci)
        tdecl = np.arcsin(rhoeci[2] / magrhoeci)

    if (trtasc < 0.0):
        trtasc = trtasc + 2.0 * np.pi

    # ------ calculate range, azimuth and elevation rates ---------
    temp1 = - rhoeci[1] * rhoeci[1] - rhoeci[0] * rhoeci[0]
    drho = np.dot(rhoeci, drhoeci) / rho
    if (np.abs(temp1) > small):
        dtrtasc = (drhoeci[0] * rhoeci[1] - drhoeci[1] * rhoeci[0]) / temp1
    else:
        dtrtasc = 0.0

    if (np.abs(temp) > small):
        dtdecl = (drhoeci[2] - drho * np.sin(tdecl)) / temp
    else:
        dtdecl = 0.0

    return rho, trtasc, tdecl, drho, dtrtasc, dtdecl



# ------------------------------------------------------------------------------
#
#                           function rv2rsw
#
#  this function converts position and velocity vectors into radial, tangential (in-
#    track), and normal (cross-track) coordinates. note that there are numerous
#    nomenclatures for these systems. this is the rsw system of vallado. the reverse
#    values are found using the transmat transpose.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#                -
#  inputs          description                    range / units
#    reci        - position vector                km
#    veci        - velocity vector                km/s
#
#  outputs       :
#    rrsw        - position vector                km
#    vrsw        - velocity vector                km/s
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


def rv2rsw(reci=None, veci=None):
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
    rrsw = smu.matvecmult(transmat, reci, 3)
    vrsw = smu.matvecmult(transmat, veci, 3)
    #   alt approach
#       rrsw[0] = smu.mag(reci)
#       rrsw[1] = 0.0
#       rrsw[2] = 0.0
#       vrsw[0] = np.dot(reci, veci)/rrsw[0]
#       vrsw[1] = sqrt(veci[0]**2 + veci[1]**2 + veci[2]**2 - vrsw[0]**2)
#       vrsw[2] = 0.0
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
#    rntw        - position vector                km
#    vntw        - velocity vector                km/s
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


def rv2ntw(r=None, v=None):
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
#                           function rv2flt.m
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
#    r           - eci position vector            km
#    v           - eci velocity vector            km/s
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    terms       - number of terms for ast calculation 0, 2
#    ddpsi, ddeps - corrections for fk5 to gcrf    rad
#
#  outputs       :
#    magr        - eci position vector magnitude  km
#    magv        - eci velocity vector magnitude  km/sec
#    latgc       - geocentric latitude            rad
#    lon         - longitude                      rad
#    fpa         - sat flight path angle          rad
#    az          - sat flight path az             rad
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
# [lon, latgc, rtasc, decl, fpa, az, magr, magv] = rv2flt (r, v, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ----------------------------------------------------------------------------

def rv2flt(reci=None, veci=None, ttt=None, jdut1=None,
           lod=None, xp=None, yp=None, terms=None, ddpsi=None,
           ddeps=None):
    small = 1e-08
    magr = smu.mag(reci)
    magv = smu.mag(veci)
    # -------- convert r to ecef for lat/lon calculation
    avec = np.zeros(3)
    recef, vecef, aecef = eci2ecef(reci, veci, avec, ttt,
                                 jdut1, lod, xp, yp, terms, ddpsi, ddeps)
    # ----------------- find longitude value  ----------------- uses ecef
    temp = np.sqrt(recef[0] * recef[0] + recef[1] * recef[1])
    if (temp < small):
        lon = math.atan2(vecef[1], vecef[0])
    else:
        lon = math.atan2(recef[1], recef[0])

    #latgc = math.atan2(recef[2] , sqrt(recef[0]**2 + recef[1]**2))
    latgc = np.arcsin(recef[2] / magr)
    # ------------- calculate rtasc and decl ------------------ uses eci
    temp = np.sqrt(reci[0] * reci[0] + reci[1] * reci[1])
    if (temp < small):
        rtasc = math.atan2(veci[1], veci[0])
    else:
        rtasc = math.atan2(reci[1], reci[0])

    #decl = math.atan2(reci[2] , sqrt(reci[0]**2 + reci[1]**2))
    decl = np.arcsin(reci[2] / magr)
    h = np.cross(reci, veci)
    hmag = smu.mag(h)
    rdotv = np.dot(reci, veci)
    fpav = math.atan2(hmag, rdotv)
    fpa = np.pi * 0.5 - fpav
    hcrossr = np.cross(h, reci)
    az = math.atan2(reci[0] * hcrossr[1] - reci[1] * hcrossr[0],
                    hcrossr[2] * magr)
    return lon, latgc, rtasc, decl, fpa, az, magr, magv


# ----------------------------------------------------------------------------
#
#                           function rv2eq.m
#
#  this function transforms a position and velocity vector into the flight
#    elements - latgc, lon, fpa, az, position and velocity magnitude.
#
#  author        : david vallado                  719-573-2600    7 jun 2002
#
#  revisions
#    vallado     - fix special orbit types (ee)                   5 sep 2002
#    vallado     - add constant file use                         29 jun 2003
#
#  inputs          description                               range / units
#    r           - eci position vector                       km
#    v           - eci velocity vector                       km/s
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

def rv2eq(r=None, v=None):
    # -------- convert to classical elements ----------------------
    p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(r, v)
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
        magr = smu.mag(r)
        magv = smu.mag(v)
        a = 1.0 / (2.0 / magr - magv ** 2 / mu)
        n = np.sqrt(mu / (a * a * a))
        wvec = np.cross(r, v) / smu.mag(np.cross(r, v))
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
        evec = - r / magr + np.cross(v, np.cross(r, v)) / mu
        ag = np.dot(evec, gvec)
        af = np.dot(evec, fvec)
        X = np.dot(r, fvec)
        Y = np.dot(r, gvec)
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



#
# ----------------------------------------------------------------------------
#
#                           function rv2adbar.m
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
#    r           - eci position vector            km
#    v           - eci velocity vector            km/s
#
#  outputs       :
#    rmag        - eci position vector magnitude  km
#    vmag        - eci velocity vector magnitude  km/sec
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

def rv2adbar(r=None, v=None):
    small = 1e-08
    rmag = smu.mag(r)
    vmag = smu.mag(v)
    # ---------------- calculate rtasc and decl -------------------
    temp = np.sqrt(r[0] * r[0] + r[1] * r[1])
    if (temp < small):
        temp1 = np.sqrt(v[0] * v[0] + v[1] * v[1])
        if (np.abs(temp1) > small):
            rtasc = math.atan2(v[1], v[0])
        else:
            rtasc = 0.0
    else:
        rtasc = math.atan2(r[1], r[0])

    decl = np.arcsin(r[2] / rmag)
    h = np.cross(r, v)
    hmag = smu.mag(h)
    rdotv = np.dot(r, v)
    fpav = math.atan2(hmag, rdotv)
    hcrossr = np.cross(h, r)
    az = math.atan2(r[0] * hcrossr[1] - r[1] * hcrossr[0], hcrossr[2] * rmag)
    return rmag, vmag, rtasc, decl, fpav, az


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

def printcov(covin=None, covtype=None, cu=None, anom=None):
    print("in printcov, anom is ", anom)
    if (str(anom) == str('truea')) or (str(anom) == str('meana')):
        semi = 'a m  '
    else:
        if (str(anom) == str('truen')) or (str(anom) == str('meann')):
            semi = 'n rad'

    if str(covtype) == str('ct'):
        print('cartesian covariance \n' % ())
        print('        x  m            y m             z  m           xdot  m/s       ydot  m/s       zdot  m/s  \n' % ())

    if str(covtype) == str('cl'):
        print('classical covariance \n' % ())
        if (cu == 'm'):
            print('          %s          ecc           incl rad      raan rad         argp rad        ' % (semi))
            if (str(anom) == str('meana')) or (str(anom) == str('meann')):
                print(' m rad \n' % ())
            else:
                if (str(anom) == str('truea')) or (str(anom) == str('truen')):
                    print(' nu rad \n' % ())
        else:
            print('          %s           ecc           incl deg      raan deg         argp deg        ' % (semi))
            if (str(anom) == str('meana')) or (str(anom) == str('meann')):
                print(' m deg \n' % ())
            else:
                if (str(anom) == str('truea')) or (str(anom) == str('truen')):
                    print(' nu deg \n' % ())

    if str(covtype) == str('eq'):
        print('equinoctial covariance \n' % ())
        #            if (cu == 'm')
        if (str(anom) == str('meana')) or (str(anom) == str('meann')):
            print('         %5s           af              ag           chi             psi         meanlonM rad\n' % (semi))
        else:
            if (str(anom) == str('truea')) or (str(anom) == str('truen')):
                print('         %5s           af              ag           chi             psi         meanlonNu rad\n' % (semi))

    if str(covtype) == str('fl'):
        print('flight covariance \n' % ())
        if (cu == 'm'):
            print('       lon  rad      latgc rad        fpa rad         az rad           r  m           v  m/s  \n' % ())
        else:
            print('       lon  deg      latgc deg        fpa deg         az deg           r  m           v  m/s  \n' % ())

    if str(covtype) == str('sp'):
        print('spherical covariance \n' % ())
        if (cu == 'm'):
            print('      rtasc rad       decl rad        fpa rad         az rad           r  m           v  m/s  \n' % ())
        else:
            print('      rtasc deg       decl deg        fpa deg         az deg           r  m           v  m/s  \n' % ())

    print("covin:")
    print(covin)
    #        print('#16e#16e#16e#16e#16e#16e\n', covin')
    print('covin transpose\n', (np.transpose(covin)))


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

def setcov(reci=None, veci=None, year=None, mon=None, day=None,
           hr=None, min=None, sec=None, dut1=None, dat=None, ttt=None,
           jdut1=None, lod=None, xp=None, yp=None, terms=None,
           printopt=None, anom=None, anom1=None, ddpsi=None, ddeps=None):
    fr = 1.0
    cartstate = np.array([[reci[0]], [reci[1]], [reci[2]],
                          [veci[0]], [veci[1]], [veci[2]]])
    print("in setcov, reci aand veci are:")
    print(reci)
    print(veci)

    # -------- convert to a classical orbit state
    p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(reci, veci)
    classstate = np.zeros(6)
    classstate[0] = a * 1000
    classstate[1] = ecc
    classstate[2] = incl
    classstate[3] = omega
    classstate[4] = argp
    if (str(anom) == str('meana')) or (str(anom) == str('meann')):
        classstate[5] = m
    else:
        if (str(anom) == str('truea')) or (str(anom) == str('truen')):
            classstate[5] = nu

    # -------- convert to a flight orbit state
    lon, latgc, rtasc, decl, fpa, az, magr, magv = rv2flt(reci, veci, ttt, jdut1, lod,
                                                   xp, yp, terms, ddpsi, ddeps)
    flstate = np.zeros(6)
    if str(anom1) == str('radec'):
        flstate[0] = rtasc
        flstate[1] = decl
    elif str(anom1) == str('latlon'):
        flstate[0] = lon
        flstate[1] = latgc

    flstate[2] = fpa
    flstate[3] = az
    flstate[4] = magr * 1000

    flstate[5] = magv * 1000
    # test position and velocity going back
    avec = np.array([[0.0], [0.0], [0.0]])
    recef, vecef, aecef = eci2ecef(reci, veci, avec, ttt, jdut1,
                                 lod, xp, yp, terms, ddpsi, ddeps)
    vx = magv * (- np.cos(lon) * np.sin(latgc) * np.cos(az) * np.cos(fpa)
                 - np.sin(lon) * np.sin(az) * np.cos(fpa)
                 + np.cos(lon) * np.cos(latgc) * np.sin(fpa))
    vy = magv * (- np.sin(lon) * np.sin(latgc) * np.cos(az) * np.cos(fpa)
                 + np.cos(lon) * np.sin(az) * np.cos(fpa)
                 + np.sin(lon) * np.cos(latgc) * np.sin(fpa))
    vz = magv * (np.sin(latgc) * np.sin(fpa)
                 + np.cos(latgc) * np.cos(az) * np.cos(fpa))
    # correct:
    ve1 = magv * (- np.cos(rtasc) * np.sin(decl) * np.cos(az) * np.cos(fpa)
                  - np.sin(rtasc) * np.sin(az) * np.cos(fpa)
                  + np.cos(rtasc) * np.cos(decl) * np.sin(fpa))

    ve2 = magv * (- np.sin(rtasc) * np.sin(decl) * np.cos(az) * np.cos(fpa)
                  + np.cos(rtasc) * np.sin(az) * np.cos(fpa)
                  + np.sin(rtasc) * np.cos(decl) * np.sin(fpa))
    ve3 = magv * (np.sin(decl) * np.sin(fpa)
                  + np.cos(decl) * np.cos(az) * np.cos(fpa))
    # -------- convert to an equinoctial orbit state
    a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr = rv2eq(reci, veci)
    eqstate = np.zeros(6)
    if (str(anom) == str('meana')) or (str(anom) == str('truea')):
        eqstate[0] = a
    else:
        if (str(anom) == str('meann')) or (str(anom) == str('truen')):
            eqstate[0] = n

    eqstate[1] = af
    eqstate[2] = ag
    eqstate[3] = chi
    eqstate[4] = psi
    if (str(anom) == str('meana')) or (str(anom) == str('meann')):
        eqstate[5] = meanlonM
    else:
        if (str(anom) == str('truea')) or (str(anom) == str('truen')):
            eqstate[5] = meanlonNu

    if printopt == 'y':
        # --------------------- write out input data --------------------------
# test velocity prints
        vtecef = np.array([vx, vy, vz])
        vteci = np.array([ve1, ve2, ve3])

        if smu.mag(vteci - np.transpose(veci)) > 0.01:
            print('ERROR in test of vel in setcov %11.7f \n'
                  % (smu.mag(vteci - np.transpose(veci))))
        print('input data \n' % ())
        print(' re %8.6f km  \n' % (re))
        print(' mu %14.8f km3/s2  \n' % (mu))
        print('year %5i ' % (year))
        print('mon %4i ' % (mon))
        print('day %3i ' % (day))
        print('hr %3i:%2i:%8.6f\n' % (hr, min, sec))
        print('dut1 %8.6f s' % (dut1))
        print(' dat %3i s' % (dat))
        print(' xp %8.6f "' % (xp))
        print(' yp %8.6f "' % (yp))
        print(' lod %8.6f s\n' % (lod))
        print('r :')
        print(reci)
        print('v :')
        print(veci)
        print('          p km       a km      ecc      incl deg    ' % ())
        print(' raan deg     argp deg      nu deg      m deg \n' % ())
        print('coes %11.4f %11.4f %11.7f %11.5f %11.5f'
              % (p, a, ecc, incl * rad2deg, omega * rad2deg))
        print('%11.5f %11.5f %11.5f\n' % (argp * rad2deg, nu * rad2deg, m * rad2deg))
        print('          a (km)       af           ag' % ())
        print('           chi        psi            meanlonM     meanLonNu   fr\n' % ())
        print('eq   %14.7f %14.7f %14.7f %15.7f %14.7f %14.7f %14.7f %2.0f\n'
              % (a, af, ag, chi, psi, meanlonM * rad2deg, meanlonNu * rad2deg, fr))
        print('       lon deg       latgc deg     rtasc deg      decl deg      fpa deg       ' % ())
        print(' az deg       magr km      magv km/s\n' % ())
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
#    constastro
#
#  references    :
#    Vallado and Alfano 2015
#
#   [cartcov, tm] = covcl2ct(classcov, classstate, anom)
# ----------------------------------------------------------------------------

def covcl2ctnew(classcov=None, classstate=None, anom=None):
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
    if (str(anom) == str('meana')) or (str(anom) == str('meann')):
        mean = classstate[5]
        e, nu = smu.newtonm(ecc, mean)
    elif (str(anom) == str('truea')) or (str(anom) == str('truen')):
        # note that mean is not used in the partials, but nu is!
        nu = classstate[5]
        e, mean = smu.newtonnu(ecc, nu)

    p = a * (1 - ecc ** 2) / 1000

    r, v = coe2rv(p, ecc, incl, raan, argp, nu, 0.0, 0.0, 0.0)
    rx = r[0] * 1000
    ry = r[1] * 1000
    rz = r[2] * 1000
    vx = v[0] * 1000
    vy = v[1] * 1000
    vz = v[2] * 1000
    # assign trig values for efficiency
    sin_inc = np.sin(incl)
    cos_inc = np.cos(incl)
    sin_raan = np.sin(raan)
    cos_raan = np.cos(raan)
    sin_w = np.sin(argp)
    cos_w = np.cos(argp)
    sin_nu = np.sin(nu)
    cos_nu = np.cos(nu)
    # assign elements of PQW to ECI transformation (pg 168)
    p11 = cos_raan * cos_w - sin_raan * sin_w * cos_inc
    p12 = - cos_raan * sin_w - sin_raan * cos_w * cos_inc
    p13 = sin_raan * sin_inc
    p21 = sin_raan * cos_w + cos_raan * sin_w * cos_inc
    p22 = - sin_raan * sin_w + cos_raan * cos_w * cos_inc
    p23 = - cos_raan * sin_inc
    p31 = sin_w * sin_inc
    p32 = cos_w * sin_inc
    p33 = cos_inc
    # assign constants for efficiency
    p0 = np.sqrt(mum / (a * (1.0 - ecc * ecc)))
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
    if (str(anom) == str('meana')) or (str(anom) == str('meann')):
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
    if ((str(anom) == str('meana')) or (str(anom) == str('meann'))):
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
    if (str(anom) == str('meana')) or (str(anom) == str('meann')):
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
    if (str(anom) == str('meana')) or (str(anom) == str('meann')):
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
    if (str(anom) == str('meana')) or (str(anom) == str('meann')):
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
    if (str(anom) == str('meana')) or (str(anom) == str('meann')):
        #tm[5, 5] = tm[5, 5] / dMdnu + tm[5, 1] / dMde
#tm[5, 1] = tm[5, 5] * dMde - dMde*atm/dMdnu
        tm[5, 5] = tm[5, 5] / dMdnu
        tm[5, 1] = tm[5, 1] - tm[5, 5] * dMde

    # ---------- calculate the output covariance matrix -----------
    cartcov = tm @ classcov @ np.transpose(tm)
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
#    constastro
#    rv2coe      - position and velocity vectors to classical elements
#
#  references    :
#    Vallado and Alfano 2015
#
#   [classcov, tm] = covct2cl(cartcov, cartstate, anom)
# ----------------------------------------------------------------------------

def covct2clnew(cartcov=None, cartstate=None, anom=None):
    # -------- define gravitational constant

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
    p, a, ecc, incl, omega, argp, nu, mean, arglat, truelon, lonper = rv2coe(reci, veci)
    p = p * 1000.0
    a = a * 1000.0
    n = np.sqrt(mum / a ** 3)
    # -------- calculate common quantities
    sqrt1me2 = np.sqrt(1.0 - ecc * ecc)
    magr = np.sqrt(rx ** 2 + ry ** 2 + rz ** 2)
    magr3 = magr ** 3
    magv = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    # ----------  form pqw position and velocity vectors ----------
    r_dot_v = np.dot(reci, veci) * 1000 * 1000
    ecc_term = magv * magv - mum / magr
    ecc_x = (ecc_term * rx - r_dot_v * vx) / mum
    ecc_y = (ecc_term * ry - r_dot_v * vy) / mum
    ecc_z = (ecc_term * rz - r_dot_v * vz) / mum
    ecc_vec = np.transpose(np.array([ecc_x, ecc_y, ecc_z]))
    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    h_vec = np.transpose(np.array([hx, hy, hz]))
    h = smu.mag(h_vec)
    h_squared = h * h
    nx = - hy
    ny = hx
    nz = 0.0
    node_vec = np.transpose(np.array([nx, ny, nz]))
    node = smu.mag(node_vec)
    n_squared = node * node
    n_dot_e = np.dot(node_vec, ecc_vec)
    sign_anode = np.sign(ny)
    cos_anode = nx / node
    omega = sign_anode * np.arccos(cos_anode)
    sign_w = np.sign((magv ** 2 - mum / magr) * rz - r_dot_v * vz)
    cos_w = n_dot_e / (ecc * node)
    argp = sign_w * np.arccos(cos_w)
    w_scale = - sign_w / np.sqrt(1 - cos_w * cos_w)
    r_dot_e = np.dot(reci, ecc_vec) * 1000
    cos_nu = r_dot_e / (magr * ecc)
    sign_nu = np.sign(r_dot_v)
    nu = sign_nu * np.arccos(cos_nu)
    nu_scale = - sign_nu / np.sqrt(1 - cos_nu * cos_nu)
    px, ax, eccx, inclx, nodex, argpx, nux, mx, arglatx, truelonx, lonperx = rv2coe(reci, veci)
    print(' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f \n'
          % (ax, eccx, inclx * rad2deg, nodex * rad2deg, argpx * rad2deg, nux * rad2deg, mx * rad2deg))
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
    tm[1, 0] = - p0 * (((vx * vy - mum * rx * ry / magr3) * ecc_y)
                      + ((vx * vz - mum * rx * rz / magr3) * ecc_z)
                      - (vy * vy + vz * vz - mum / magr + mum * rx * rx / magr3)
                      * ecc_x)
    tm[1, 1] = - p0 * (((vx * vy - mum * rx * ry / magr3) * ecc_x)
                      + ((vy * vz - mum * ry * rz / magr3) * ecc_z)
                      - (vx * vx + vz * vz - mum / magr + mum * ry * ry / magr3)
                      * ecc_y)
    tm[1, 2] = - p0 * (((vx * vz - mum * rx * rz / magr3) * ecc_x)
                      + ((vy * vz - mum * ry * rz / magr3) * ecc_y)
                      - (vy * vy + vx * vx - mum / magr + mum * rz * rz / magr3)
                      * ecc_z)
    tm[1, 3] = - p0 * ((rx * vy - 2 * ry * vx) * ecc_y
                      + (ry * vy + rz * vz) * ecc_x
                      + (rx * vz - 2 * rz * vx) * ecc_z)
    tm[1, 4] = - p0 * ((ry * vx - 2 * rx * vy) * ecc_x
                      + (rx * vx + rz * vz) * ecc_y
                      + (ry * vz - 2 * rz * vy) * ecc_z)
    tm[1, 5] = - p0 * ((rx * vx + ry * vy) * ecc_z
                      + (rz * vx - 2 * rx * vz) * ecc_x
                      + (rz * vy - 2 * ry * vz) * ecc_y)
    # ---- partials of incl wrt (rx ry rz vx vy vz)
    p3 = 1.0 / node
    tm[2, 0] = - p3 * (vy - hz * (vy * hz - vz * hy) / h_squared)
    tm[2, 1] = p3 * (vx - hz * (vx * hz - vz * hx) / h_squared)
    tm[2, 2] = - p3 * (hz * (vy * hx - vx * hy) / h_squared)
    tm[2, 3] = p3 * (ry - hz * (ry * hz - rz * hy) / h_squared)
    tm[2, 4] = - p3 * (rx - hz * (rx * hz - rz * hx) / h_squared)
    tm[2, 5] = p3 * (hz * (ry * hx - rx * hy) / h_squared)

    # ---- partials of node wrt (rx ry rz vx vy vz)
    p4 = 1.0 / n_squared
    tm[3, 0] = - p4 * vz * ny
    tm[3, 1] = p4 * vz * nx
    tm[3, 2] = p4 * (vx * ny - vy * nx)
    tm[3, 3] = p4 * rz * ny
    tm[3, 4] = - p4 * rz * nx
    tm[3, 5] = p4 * (ry * nx - rx * ny)
    # ---- partials of argp wrt (rx ry rz vx vy vz)
    p5 = 1.0 / (node * a * a)
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

    if (str(anom) == str('meana')) or (str(anom) == str('meann')):
        # ---- partials of (rx ry rz vx vy vz) wrt mean anomaly
# then update for mean anomaly
        ecc = smu.mag(ecc_vec)
        dMdnu = (1.0 - ecc ** 2) ** 1.5 / ((1.0 + ecc * np.cos(nu)) ** 2)
        dMde = (-np.sin(nu) * ((ecc * np.cos(nu) + 1)
                               * (ecc + np.cos(nu))
                               / np.sqrt((ecc + np.cos(nu)) ** 2)
                               + 1.0 - 2.0 * ecc ** 2 - ecc ** 3 * np.cos(nu))
                / ((ecc * np.cos(nu) + 1.0) ** 2 * np.sqrt(1 - ecc ** 2)))
        # p6 = -sin(nu)*(sign((ecc*cos(nu) + 1)) + 1.0 - 2.0*ecc**2 - ecc**3*cos(nu)) / ((ecc*cos(nu) + 1.0)**2 * sqrt(1-ecc**2))  # dm/de
        tm[5, 0] = tm[5, 0] * dMdnu + tm[1, 0] * dMde
        tm[5, 1] = tm[5, 1] * dMdnu + tm[1, 1] * dMde
        tm[5, 2] = tm[5, 2] * dMdnu + tm[1, 2] * dMde
        tm[5, 3] = tm[5, 3] * dMdnu + tm[1, 3] * dMde
        tm[5, 4] = tm[5, 4] * dMdnu + tm[1, 4] * dMde
        tm[5, 5] = tm[5, 5] * dMdnu + tm[1, 5] * dMde

    # ---------- calculate the output covariance matrix -----------
    classcov = tm * cartcov * np.transpose(tm)
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
#
#  locals        :
#    r           - position vector                m
#    v           - velocity vector                m/s
#    tm          - transformation matrix
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

def covct2ntw(cartcov=None, cartstate=None):
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
    covntw = tm * cartcov * np.transpose(tm)
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
#
#  locals        :
#    r           - position vector                m
#    v           - velocity vector                m/s
#    tm          - transformation matrix
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

def covct2o2(cartcov=None, cartstate=None):
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
    covopntw = tm * cartcov * np.transpose(tm)
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
#
#  locals        :
#    r           - position vector                m
#    v           - velocity vector                m/s
#    tm          - transformation matrix
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

def covct2rsw(cartcov=None, cartstate=None):
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
    covoprsw = tm * cartcov * np.transpose(tm)
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
#
#  locals        :
#    r           - position vector                m
#    v           - velocity vector                m/s
#    tm          - transformation matrix
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

def covo22ct(covopntw=None, cartstate=None):
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
    tm = np.transpose(tm)
    cartcov = tm * covopntw * np.transpose(tm)
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
#    rmag        - eci position vector magnitude  km
#    vmag        - eci velocity vector magnitude  km/sec
#    rtasc       - right ascension of sateillite  rad
#    decl        - declination of satellite       rad
#    fpav        - sat flight path angle from vertrad
#    az          - sat flight path azimuth        rad
#
#  outputs       :
#    r           - eci position vector            km
#    v           - eci velocity vector            km/s
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

def adbar2rv(rmag=None, vmag=None, rtasc=None, decl=None, fpav=None, az=None):
    # -------- form position vector
    r = np.zeros(3)
    r[0] = rmag * np.cos(decl) * np.cos(rtasc)
    r[1] = rmag * np.cos(decl) * np.sin(rtasc)
    r[2] = rmag * np.sin(decl)
    # -------- form velocity vector
    v = np.zeros(3)
    v[0] = vmag * (np.cos(rtasc) * (- np.cos(az) * np.sin(fpav) * np.sin(decl)
                                    + np.cos(fpav) * np.cos(decl))
                   - np.sin(az) * np.sin(fpav) * np.sin(rtasc))
    v[1] = vmag * (np.sin(rtasc) * (- np.cos(az) * np.sin(fpav) * np.sin(decl)
                                    + np.cos(fpav) * np.cos(decl))
                   + np.sin(az) * np.sin(fpav) * np.cos(rtasc))
    v[2] = vmag * (np.cos(az) * np.cos(decl) * np.sin(fpav)
                   + np.cos(fpav) * np.sin(decl))
    r = np.transpose(r)
    v = np.transpose(v)
    return r, v



# azl2radc
#
# this function finds the rtasc decl values given the az-el
#
#
#

def azl2radc(az=None, el=None, lat=None, lst=None):
    decl = np.arcsin(np.sin(el) * np.sin(lat)
                     + np.cos(el) * np.cos(lat) * np.cos(az))
    slha1 = (- (np.sin(az) * np.cos(el) * np.cos(lat))
             / (np.cos(decl) * np.cos(lat)))
    clha1 = ((np.sin(el) - np.sin(lat) * np.sin(decl))
             / (np.cos(decl) * np.cos(lat)))
    lha1 = math.atan2(slha1, clha1)
    #    print(' lha1 #13.7f \n', lha1 * rad2deg)

    # alt approach
    slha2 = - (np.sin(az) * np.cos(el)) / (np.cos(decl))
    clha2 = ((np.cos(lat) * np.sin(el) - np.sin(lat) * np.cos(el) * np.cos(az))
             / (np.cos(decl)))
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
#    del         - elevation rate                 rad / s
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

def raz2rvs(rho=None, az=None, el=None, drho=None, daz=None, del_=None):
    # ----------------------- initialize values -------------------
    sinel = np.sin(el)
    cosel = np.cos(el)
    sinaz = np.sin(az)
    cosaz = np.cos(az)
    # ------------------- form sez range vector -------------------
    rhosez = np.zeros(3)
    rhosez[0] = - rho * cosel * cosaz
    rhosez[1] = rho * cosel * sinaz
    rhosez[2] = rho * sinel
    # ----------------- form sez velocity vector ------------------
    drhosez = np.zeros(3)
    drhosez[0] = - drho * cosel * cosaz + rhosez[2] * del_ * cosaz + rhosez[1] * daz
    drhosez[1] = drho * cosel * sinaz - rhosez[2] * del_ * sinaz - rhosez[0] * daz
    drhosez[2] = drho * sinel + rho * del_ * cosel
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
#    del         - elevation rate                 rad / s
#    rs          - ecef site position vector      km
#    latgd       - geodetic latitude              -pi/2 to pi/2 rad
#    lon         - longitude of site              -2pi to 2pi rad
#    alt         - altitude                       km
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    terms       - number of terms for ast calculation 0, 2
#
#  outputs       :
#    reci        - eci position vector            km
#    veci        - eci velocity vector            km/s
#
#  locals        :
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

def razel2rv(rho=None, az=None, el=None, drho=None, daz=None,
             del_=None, latgd=None, lon=None, alt=None, ttt=None,
             jdut1=None, lod=None, xp=None, yp=None, terms=None,
             ddpsi=None, ddeps=None):
    # -------------------------  implementation   -----------------
    # -----------  find sez range and velocity vectors ------------
    rhosez, drhosez = raz2rvs(rho, az, el, drho, daz, del_)
    # -----------  perform sez to ijk (ecef) transformation -------
    tempvec = smu.rot2(rhosez, latgd - halfpi)
    rhoecef = smu.rot3(tempvec, -lon)
    rhoecef = np.transpose(rhoecef)
    tempvec = smu.rot2(drhosez, latgd - halfpi)
    drhoecef = smu.rot3(tempvec, -lon)
    drhoecef = np.transpose(drhoecef)
    # ----------  find ecef range and velocity vectors -------------
    rs, vs = obu.site(latgd, lon, alt)
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
    reci, veci, aeci = ecef2eci(recef, vecef, a, ttt, jdut1,
                              lod, xp, yp, terms, ddpsi, ddeps)
    return reci, veci

#
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

def coe2rvS(p=None, ecc=None, incl=None, omega=None, argp=None,
            nu=None, arglat=None, truelon=None, lonper=None):
    # -------------------------------------------------------------
#       determine what type of orbit is involved and set up the
#       set up angles for the special cases.
# -------------------------------------------------------------
    if (ecc < small):
        # ----------------  circular equatorial  ------------------
        if ((incl < small) or (np.abs(incl - np.pi) < small)):
            argp = 0.0
            omega = 0.0
            nu = truelon
        else:
            # --------------  circular inclined  ------------------
            argp = 0.0
            nu = arglat
    else:
        # ---------------  elliptical equatorial  -----------------
        if (((incl < small) or (np.abs(incl - np.pi) < small))):
            argp = lonper
            omega = 0.0

    # ----------  form pqw position and velocity vectors ----------
    cosnu = np.cos(nu)
    sinnu = np.sin(nu)
    temp = p / (1.0 + ecc * cosnu)
    rpqw = np.zeros(3)
    rpqw[0] = temp * cosnu
    rpqw[1] = temp * sinnu
    rpqw[2] = 0.0
    rpqw
    if (np.abs(p) < 0.0001):
        p = 0.0001

    vpqw = np.zeros(3)
    vpqw[0] = - sinnu * np.sqrt(mu) / np.sqrt(p)
    vpqw[1] = (ecc + cosnu) * np.sqrt(mu) / np.sqrt(p)
    vpqw[2] = 0.0
    vpqw
    # ----------------  perform transformation to ijk  ------------
    tempvec = smu.rot3(rpqw, - argp)
    tempvec = smu.rot1(tempvec, - incl)
    r = smu.rot3(tempvec, - omega)
    tempvec = smu.rot3(vpqw, - argp)
    tempvec = smu.rot1(tempvec, - incl)
    v = smu.rot3(tempvec, - omega)
    r = np.transpose(r)
    v = np.transpose(v)
    return r, v


#
# ------------------------------------------------------------------------------
#
#                           function coe2rvh
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
#    mu          - gravitational parameter in km3/s2
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
# [r, v] = coe2rv (p, ecc, incl, omega, argp, nu, arglat, truelon, lonper, mu)
# ------------------------------------------------------------------------------

def coe2rvh(p=None, ecc=None, incl=None, omega=None, argp=None, nu=None,
            arglat=None, truelon=None, lonper=None, mu=None):
    # -------------------------  implementation   -----------------

    # -------------------------------------------------------------
#       determine what type of orbit is involved and set up the
#       set up angles for the special cases.
# -------------------------------------------------------------
    if (ecc < small):
        # ----------------  circular equatorial  ------------------
        if (incl < small) or (np.abs(incl - np.pi) < small):
            argp = 0.0
            omega = 0.0
            nu = truelon
        else:
            # --------------  circular inclined  ------------------
            argp = 0.0
            nu = arglat
    else:
        # ---------------  elliptical equatorial  -----------------
        if ((incl < small) or (np.abs(incl - np.pi) < small)):
            argp = lonper
            omega = 0.0

    # ----------  form pqw position and velocity vectors ----------
    cosnu = np.cos(nu)
    sinnu = np.sin(nu)
    temp = p / (1.0 + ecc * cosnu)
    rpqw = np.zeros(3)
    rpqw[0] = temp * cosnu
    rpqw[1] = temp * sinnu
    rpqw[2] = 0.0
    if (np.abs(p) < 0.0001):
        p = 0.0001

    vpqw = np.zeros(3)
    vpqw[0] = - sinnu * np.sqrt(mu) / np.sqrt(p)
    vpqw[1] = (ecc + cosnu) * np.sqrt(mu) / np.sqrt(p)
    vpqw[2] = 0.0
    # ----------------  perform transformation to ijk  ------------
    tempvec = smu.rot3(rpqw, - argp)
    tempvec = smu.rot1(tempvec, - incl)
    r = smu.rot3(tempvec, - omega)
    tempvec = smu.rot3(vpqw, - argp)
    tempvec = smu.rot1(tempvec, - incl)
    v = smu.rot3(tempvec, - omega)
    r = np.transpose(r)
    v = np.transpose(v)
    return r, v



# ------------------------------------------------------------------------------
#
#                           function rv2coeS
#
#  this function finds the classical orbital elements given the geocentric
#    equatorial position and velocity vectors. Show interm calcs
#
#  author        : david vallado                  719-573-2600   21 jun 2002
#
#  revisions
#    vallado     - fix special cases                              5 sep 2002
#    vallado     - delete extra check in inclination code        16 oct 2002
#    vallado     - add constant file use                         29 jun 2003
#    vallado     - add mu                                         2 apr 2007
#
#  inputs          description                    range / units
#    r           - ijk position vector            km
#    v           - ijk velocity vector            km / s
#    mu          - gravitational parameter        km3 / s2
#
#  outputs       :
#    p           - semilatus rectum               km
#    a           - semimajor axis                 km
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
#  locals        :
#    hbar        - angular momentum h vector      km2 / s
#    ebar        - eccentricity     e vector
#    nbar        - line of nodes    n vector
#    c1          - v**2 - u/r
#    rdotv       - r dot v
#    hk          - hk unit vector
#    sme         - specfic mechanical energy      km2 / s2
#    i           - index
#    e           - eccentric, parabolic,
#                  hyperbolic anomaly             rad
#    temp        - temporary variable
#    typeorbit   - type of orbit                  ee, ei, ce, ci
#
#  coupling      :
#    mag         - magnitude of a vector
#    angl        - find the angl between two vectors
#    newtonnu    - find the mean anomaly
#
#  references    :
#    vallado       2007, 121, alg 9, ex 2-5
#
# [p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper ] = rv2coeS (r, v)
# ------------------------------------------------------------------------------


def rv2coeS (r, v):

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
        c1 = magv*magv - mu /magr
        rdotv = np.dot(r, v)
        ebar = np.zeros((3))
        for i in range(3):
            ebar[i] = (c1*r[i] - rdotv*v[i])/mu
        #end
        ecc = smu.mag(ebar)

        # ------------  find a e and semi-latus rectum   ----------
        sme = (magv*magv*0.5) - (mu /magr)
        if (abs(sme) > small):
            a = -mu  / (2.0 *sme)
        else:
            a = infinite
        #end
        p = magh*magh/mu

        # -----------------  find inclination   -------------------
        hk = hbar[2]/magh
        incl = math.acos(hk)

        # --------  determine type of orbit for later use  --------
        # ------ elliptical, parabolic, hyperbolic inclined -------
        typeorbit = 'ei'
        if (ecc < small):
            # ----------------  circular equatorial ---------------
            if  (incl<small) or (abs(incl-math.pi)<small):
                typeorbit = 'ce'
            else:
                # --------------  circular inclined ---------------
                typeorbit = 'ci'
            #end
        else:
            # - elliptical, parabolic, hyperbolic equatorial --
            if  (incl<small) or (abs(incl-math.pi)<small):
                typeorbit = 'ee'
            #end
        #end


        print("Orbit type in rv2coeS is ", typeorbit)
        # ----------  find longitude of ascending node ------------
        if (magn > small):
            temp = nbar[0] / magn
            if (abs(temp) > 1.0):
                temp = math.sign(temp)
            #end
            omega = math.acos(temp)
            if (nbar[1] < 0.0):
                omega = twopi - omega
            #end
        else:
            omega=None
        #end

        # ---------------- find argument of perigee ---------------
        if (typeorbit == 'ei'):
            np.dot(nbar, ebar)
            argp = smu.angl(nbar, ebar)
            if (ebar[2] < 0.0):
                argp = twopi - argp
            #end
        else:
            argp=None
        #end

        # ------------  find true anomaly at epoch    -------------
        if typeorbit.startswith('e'):
            np.dot(ebar, r)
            nu = smu.angl(ebar, r)
            if (rdotv < 0.0):
                nu = twopi - nu
            #end
        else:
            nu=None
        #end

        # ----  find argument of latitude - circular inclined -----
        # and in general cases too
        if (typeorbit == 'ci') or (typeorbit == 'ei'):
            arglat = smu.angl(nbar, r)
            if (r[2] < 0.0):
                arglat = twopi - arglat
            #end
            m = arglat
        else:
            arglat=None
        #end

        # -- find longitude of perigee - elliptical equatorial ----
        if  (ecc>small) and (typeorbit == 'ee'):
            temp = ebar[0]/ecc
            if (abs(temp) > 1.0):
                temp = math.sign(temp)
            #end
            lonper = math.acos(temp)
            if (ebar[1] < 0.0):
                lonper = twopi - lonper
            #end
            if (incl > halfpi):
                lonper = twopi - lonper
            #end
        else:
            lonper=None
        #end

        # -------- find true longitude - circular equatorial ------
        if  (magr>small) and (typeorbit == 'ce'):
            temp = r[0]/magr
            if (abs(temp) > 1.0):
                temp = math.sign(temp)
            #end
            truelon = math.acos(temp)
            if (r[1] < 0.0):
                truelon = twopi - truelon
            #end
            if (incl > halfpi):
                truelon = twopi - truelon
            #end
            m = truelon
        else:
            truelon=None
        #end

        # ------------ find mean anomaly for all orbits -----------
        if typeorbit.startswith('e'):
          e, m = smu.newtonnu(ecc, nu)
       # end

    else:
       p=None
       a=None
       ecc=None
       incl=None
       omega=None
       argp=None
       nu=None
       m=None
       arglat=None
       truelon=None
       lonper=None
    #end
    return p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper




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
#  inputs          description                    range / units
#    r           - ijk position vector            km
#    v           - ijk velocity vector            km / s
#    mu          - gravitational parameter        km3 / s2
#
#  outputs       :
#    p           - semilatus rectum               km
#    a           - semimajor axis                 km
#    ecc         - eccentricity
#    incl        - inclination                    0.0  to pi rad
#    raan       - longitude of ascending node    0.0  to 2pi rad
#    argp        - argument of perigee            0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#    m           - mean anomaly                   0.0  to 2pi rad
#    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
#    truelon     - true longitude            (ce) 0.0  to 2pi rad
#    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
#
#  locals        :
#    hbar        - angular momentum h vector      km2 / s
#    ebar        - eccentricity     e vector
#    nbar        - line of nodes    n vector
#    c1          - v**2 - u/r
#    rdotv       - r dot v
#    hk          - hk unit vector
#    sme         - specfic mechanical energy      km2 / s2
#    i           - index
#    e           - eccentric, parabolic,
#                  hyperbolic anomaly             rad
#    temp        - temporary variable
#    typeorbit   - type of orbit                  ee, ei, ce, ci
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


def rv2coe (r, v):

    muin = mu # this is the km version
    small = 1.0e-12
    # -------------------------  implementation   -----------------
    magr = smu.mag(r)
    magv = smu.mag(v)
    # ------------------  find h n and e vectors   ----------------
    hbar = np.cross(r, v)
    magh = smu.mag(hbar)
    if (magh >= 0.0):
        nbar = np.zeros((3))
        nbar[0] = -hbar[1]
        nbar[1] = hbar[0]
        nbar[2] = 0.0
        magn = smu.mag(nbar)
        c1 = (magv*magv) - (muin /magr)
        rdotv = np.dot(r, v)
        ebar = np.zeros((3))
        for i in range(3):
            ebar[i] = (c1*r[i] - rdotv*v[i])/muin
        #end
        ecc = smu.mag(ebar)

        # ------------  find a e and semi-latus rectum   ----------
        sme = (magv*magv*0.5) - (muin /magr)
        if (abs(sme) > small):
            a = -muin  / (2.0 *sme)
        else:
            a = infinite
        #end
        p = magh*magh/muin

        # -----------------  find inclination   -------------------

        hk = hbar[2]/magh
        incl = math.acos(hk)

        # --------  determine type of orbit for later use  --------
        # ------ elliptical, parabolic, hyperbolic inclined -------
        typeorbit = 'ei'
        if (ecc < small):
            # ----------------  circular equatorial ---------------
            if  (incl<small) or (abs(incl-math.pi)<small):
                typeorbit = 'ce'
            else:
                # --------------  circular inclined ---------------
                typeorbit = 'ci'
            #end
        else:
            # - elliptical, parabolic, hyperbolic equatorial --
            if  (incl<small) or (abs(incl-math.pi)<small):
                typeorbit = 'ee'
            #end
        #end


        #print("Orbit type in rv2coe is ", typeorbit)


        # ----------  find right ascension of ascending node ------------
        if (magn > small):
            temp = nbar[0] / magn
            if (abs(temp) > 1.0):
                temp = math.sign(temp)
            #end
            raan = math.acos(temp)
            if (nbar[1] < 0.0):
                raan = twopi - raan
            #end
        else:
            raan=None
        #end

        # ---------------- find argument of perigee ---------------
        if (typeorbit == 'ei'):
            argp = smu.angl(nbar, ebar)
            if (argp and (ebar[2] < 0.0)):
                argp = twopi - argp
            #end
        else:
            argp=None
        #end

        # ------------  find true anomaly at epoch    -------------
        if typeorbit.startswith('e'):
            nu = smu.angl(ebar, r)
        #    print("nu init is ", nu)
            if (rdotv < 0.0):
                nu = twopi - nu
        #        print("nu now is ", nu)
            #end
        # "True anomaly is not defined for circular orbits because
        # they have no periapsis. We can overcome this limitation by selecting a
        # direction in the orbit to replace periapsis as the location for the
        # initial measurement. Computer-software routines must account for this
        # special case." (pg 18)
        else:
            nu = None ####in some cases matlab was counting on this value!!!!-jmb


        #end

        # ----  find argument of latitude - circular inclined -----
        # -- find in general cases too
        if (typeorbit == 'ci') or (typeorbit == 'ei'):
            arglat = smu.angl(nbar, r)
            if (r[2] < 0.0):
                arglat = twopi - arglat
            #end
            m = arglat
        else:
            arglat=None
        #end

        # -- find longitude of perigee - elliptical equatorial ----
        if  (ecc>small) and (typeorbit == 'ee'):
            temp = ebar[0]/ecc
            if (abs(temp) > 1.0):
                temp = math.sign(temp)
            #end
            lonper = math.acos(temp)
            if (ebar[1] < 0.0):
                lonper = twopi - lonper
            #end
            if (incl > halfpi):
                lonper = twopi - lonper
            #end
        else:
            lonper=None
        #end

        # -------- find true longitude - circular equatorial ------
        if  (magr>small) and (typeorbit == 'ce'):
            temp = r[0]/magr
            if (abs(temp) > 1.0):
                temp = math.sign(temp)
            #end
            truelon = math.acos(temp)
            if (r[1] < 0.0):
                truelon = twopi - truelon
            #end
            if (incl > halfpi):
                truelon = twopi - truelon
            #end
            m = truelon
        else:
            truelon=None
        #end

        # ------------ find mean anomaly for all orbits -----------
       # if (typeorbit(1:1) == 'e')
       # print("ecc =", ecc)
       # print("nu =", nu)
        e, m = smu.newtonnu(ecc, nu)
       # end

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
    #end
    return p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper




#
# ------------------------------------------------------------------------------
#
#                           function rv2coeh
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
#  inputs          description                    range / units
#    r           - ijk position vector            km
#    v           - ijk velocity vector            km / s
#    mu          - gravitational parameter        km3 / s2
#
#  outputs       :
#    p           - semilatus rectum               km
#    a           - semimajor axis                 km
#    ecc         - eccentricity
#    incl        - inclination                    0.0  to pi rad
#    raan       - longitude of ascending node    0.0  to 2pi rad
#    argp        - argument of perigee            0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#    m           - mean anomaly                   0.0  to 2pi rad
#    arglat      - argument of latitude      (ci) 0.0  to 2pi rad
#    truelon     - true longitude            (ce) 0.0  to 2pi rad
#    lonper      - longitude of periapsis    (ee) 0.0  to 2pi rad
#
#  locals        :
#    hbar        - angular momentum h vector      km2 / s
#    ebar        - eccentricity     e vector
#    nbar        - line of nodes    n vector
#    c1          - v**2 - u/r
#    rdotv       - r dot v
#    hk          - hk unit vector
#    sme         - specfic mechanical energy      km2 / s2
#    i           - index
#    e           - eccentric, parabolic,
#                  hyperbolic anomaly             rad
#    temp        - temporary variable
#    typeorbit   - type of orbit                  ee, ei, ce, ci
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

def rv2coeh(r=None, v=None, re=None, mu=None):
    small = 1e-10
    #undefined = 999999.1
    print("in rv2coeh, mu is ", mu)
    # -------------------------  implementation   -----------------
    magr = smu.mag(r)
    magv = smu.mag(v)
    # ------------------  find h n and e vectors   ----------------
    hbar = np.cross(r, v)
    magh = smu.mag(hbar)
    if (magh > small):
        nbar = np.zeros((3))
        nbar[0] = -hbar[1]
        nbar[1] = hbar[0]
        nbar[2] = 0.0
        magn = smu.mag(nbar)
        c1 = magv * magv - mu / magr
        rdotv = np.dot(r, v)
        ebar = np.zeros((3))
        for i in range(3):
            ebar[i] = (c1 * r[i] - rdotv * v[i]) / mu
        ecc = smu.mag(ebar)
        # ------------  find a e and semi-latus rectum   ----------
        sme = (magv * magv * 0.5) - (mu / magr)
        if (np.abs(sme) > small):
            a = - mu / (2.0 * sme)
        else:
            a = infinite
        p = magh * magh / mu
        # -----------------  find inclination   -------------------
        hk = hbar[2] / magh
        incl = np.arccos(hk)
        # --------  determine type of orbit for later use  --------
# ------ elliptical, parabolic, hyperbolic inclined -------
        typeorbit = 'ei'
        if (ecc < small):
            # ----------------  circular equatorial ---------------
            if (incl < small) or (np.abs(incl - np.pi) < small):
                typeorbit = 'ce'
            else:
                # --------------  circular inclined ---------------
                typeorbit = 'ci'
        else:
            # - elliptical, parabolic, hyperbolic equatorial --
            if (incl < small) or (np.abs(incl - np.pi) < small):
                typeorbit = 'ee'

        print("tyeporbit is ",typeorbit)
        # ----------  find right ascension of ascending node ------------
        if (magn > small):
            temp = nbar[0] / magn
            if (np.abs(temp) > 1.0):
                temp = np.sign(temp)
            raan = np.arccos(temp)
            if (nbar[1] < 0.0):
                raan = twopi - raan
        else:
            raan=None
        # ---------------- find argument of perigee ---------------
        if str(typeorbit) == str('ei'):
            argp = smu.angl(nbar, ebar)
            if (ebar[2] < 0.0):
                argp = twopi - argp
        else:
            argp=None
        # ------------  find true anomaly at epoch    -------------
        if typeorbit.startswith('e'):
            nu = smu.angl(ebar, r)
            if (rdotv < 0.0):
                nu = twopi - nu
        else:
            nu=None
        # ----  find argument of latitude - circular inclined -----
# -- find in general cases too
        if (str(typeorbit) == 'ci') or (str(typeorbit) == 'ei'):
            arglat = smu.angl(nbar, r)
            if (r[2] < 0.0):
                arglat = twopi - arglat
            m = arglat
        else:
            arglat=None
        # -- find longitude of perigee - elliptical equatorial ----
        if (ecc > small) and (str(typeorbit) == 'ee'):
            temp = ebar[0] / ecc
            if (np.abs(temp) > 1.0):
                temp = np.sign(temp)
            lonper = np.arccos(temp)
            if (ebar[1] < 0.0):
                lonper = twopi - lonper
            if (incl > 0.5 * np.pi):
                lonper = twopi - lonper
        else:
            lonper=None
        # -------- find true longitude - circular equatorial ------
        if (magr > small) and (str(typeorbit) == 'ce'):
            temp = r[0] / magr
            if (np.abs(temp) > 1.0):
                temp = np.sign(temp)
            truelon = np.arccos(temp)
            if (r[1] < 0.0):
                truelon = twopi - truelon
            if (incl > halfpi):
                truelon = twopi - truelon
            m = truelon
        else:
            truelon=None
        # ------------ find mean anomaly for all orbits -----------
        if typeorbit.startswith('e'):
            e, m = smu.newtonnu(ecc, nu)
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
           nu: float, arglat: float, truelon: float, lonper: float):
    """This function finds the position and velocity vectors in geocentric
    equatorial (ijk) system given the classical orbit elements.

    Parameters
    ----------
    p: float
        semilatus rectum
    ecc: float
        eccentricity
    incl: float
        inclination
    omega: float
        longitude of ascending node
    argp: float
        argument of pedigree
    nu: float
        true anomaly
    arglat: float
        argument of latitude
    truelon: float
        true longitude
    lonper: float
        longitude of periapsis

    Returns
    -------
        ijk position vector, ijk velocity vector
    """

    # -------------------------------------------------------------
    #       determine what type of orbit is involved and set up the
    #       set up angles for the special cases.
    # -------------------------------------------------------------
    if (ecc < small):
        # ----------------  circular equatorial  ------------------
        if (incl<small) or (abs(incl-math.pi)< small):
            argp = 0.0
            omega = 0.0
            nu = truelon
        else:
            # --------------  circular inclined  ------------------
            argp = 0.0
            nu = arglat
    else:
            # ---------------  elliptical equatorial  -----------------
        if ((incl<small) or (abs(incl-math.pi)<small)):
            argp = lonper
            omega = 0.0

    #print("here, argp =%.3f, omega =%.3f, nu =%.3f" % (argp, omega, nu))
    # ----------  form pqw position and velocity vectors ----------
    cosnu = math.cos(nu)
    sinnu = math.sin(nu)
    temp = p / (1.0  + ecc*cosnu)
    rpqw = np.zeros((3))
    rpqw[0] = temp*cosnu
    rpqw[1] = temp*sinnu
    rpqw[2] = 0.0
    if (abs(p) < 0.0001):
        p = 0.0001
    vpqw = np.zeros((3))
    vpqw[0] = -sinnu*np.sqrt(mu)  / np.sqrt(p)
    vpqw[1] = (ecc + cosnu)*np.sqrt(mu) / np.sqrt(p)
    vpqw[2] = 0.0

    # ----------------  perform transformation to ijk  ------------
    tempvec = smu.rot3(rpqw   , -argp)
    tempvec = smu.rot1(tempvec, -incl)
    r = smu.rot3(tempvec, -omega)

    tempvec = smu.rot3(vpqw   , -argp)
    tempvec = smu.rot1(tempvec, -incl)
    v = smu.rot3(tempvec, -omega)

    r = r.T #r = r'
    v = v.T #v = v'
    return r, v



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


def ecef2eci  (recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps):
    # ---- find matrices
    prec, psia, wa, ea, xa = obu.precess (ttt, '80')

    deltapsi, trueeps, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)

    st, stdot = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

    pm = smu.polarm(xp, yp, ttt, '80')

    # ---- perform transformations
    thetasa = earthrot * (1.0  - lod/86400.0)
    #        omegaearth = np.array([[0], [0], [thetasa]])
    #omegaearth = np.array([[0.0], [0.0], [thetasa]])
    omegaearth = np.array([0.0, 0.0, thetasa]).T
    print("in ecef2eci:")

    tmpmat = np.matmul(np.matmul(prec, nut), st)
    rpef = np.matmul(pm, recef)
    reci = np.matmul(tmpmat, rpef)
    #prec*nut*st*pm
    vpef = np.matmul(pm, vecef)
    oresh = np.reshape(omegaearth, (3, 1))
    print(oresh)
    temp = np.cross(omegaearth, rpef.T) #turn from column to row vectors for cross product
    veci = np.dot(tmpmat, (vpef + temp.T))

    # veci1 = prec*nut * (stdot*recef + st*pm*vecef)  % alt approach using sidereal rate

    # two additional terms not needed if satellite is not on surface
    # of the Earth
    c1 = np.cross(omegaearth, temp)
    c2 = 2.0*np.cross(omegaearth.T, vpef.T)
    aeci = np.matmul(tmpmat, np.matmul(pm, aecef))
    #print("aeci: ", aeci)

    aeci = aeci + c1.T + c2.T
    return reci, veci, aeci



#  this function converts range, azimuth, and elevation values with slant
#    range and velocity vectors for a satellite from a radar site in the
#    topocentric horizon (sez) system.
#
#  author        : david vallado                  719-573-2600   10 jun 2002
#
#  inputs          description                    range / units
#    rho         - satellite range from site      km
#    az          - azimuth                        0.0 to 2pi rad
#    el          - elevation                      -pi/2 to pi/2 rad
#
#  outputs       :
#    rhovec      - sez satellite range vector     km



def raz2rv(rho, az_deg, el_deg):
  az = az_deg/360.0*math.tau
  el = el_deg/90.0*math.pi/2.0

  # ----------------------- initialize values -------------------
  sinel = math.sin(el)
  cosel = math.cos(el)
  sinaz = math.sin(az)
  cosaz = math.cos(az)

  # ------------------- form sez range vector -------------------
  rhosez = [0.0, 0.0, 0.0]
  rhosez[0] = -rho*cosel*cosaz
  rhosez[1] = rho*cosel*sinaz
  rhosez[2] = rho*sinel
  return rhosez



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
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
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


def ecef2mod  (recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps):
        # ---- find matrices
        deltapsi, trueeps, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)

        st, stdot = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

        pm = smu.polarm(xp, yp, ttt, '80')

        # ---- perform transformations
        thetasa = earthrot * (1.0  - lod/86400.0)

        omegaearth = np.array([0.0, 0.0, thetasa])

#trueeps-meaneps
#deltapsi
#nut
        rpef = pm@recef
        rmod = nut@st@rpef

        vpef = pm@vecef
        vmod = nut@st@(vpef + np.cross(omegaearth, rpef.T).T)

        temp = np.cross(omegaearth, rpef.T)
        amod = nut@st@(pm@aecef + np.cross(omegaearth, temp).T \
               + 2.0*np.cross(omegaearth, vpef.T).T)

        return rmod, vmod, amod


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


def ecef2tod(recef, vecef, aecef, ttt, jdut1,
             lod, xp, yp, eqeterms, ddpsi, ddeps):
        # ---- find matrices - note nut is only needed for st argument inputs
        deltapsi, trueeps, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)

        st, stdot = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

        pm = smu.polarm(xp, yp, ttt, '80')

        # ---- perform transformations
        thetasa = earthrot * (1.0  - lod/86400.0)
        omegaearth = np.array([0.0, 0.0, thetasa])

        rpef = pm@recef
        rtod = st@rpef

        vpef = pm@vecef
        vtod = st@(vpef + np.cross(omegaearth, rpef.T).T)

        temp = np.cross(omegaearth, rpef.T)
        atod = st@(pm@aecef + np.cross(omegaearth, temp).T \
               + 2.0*np.cross(omegaearth, vpef.T).T)

        return rtod, vtod, atod



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


def ecef2teme(recef, vecef, aecef, ttt, jdut1, lod, xp, yp, eqeterms):

        # ------------------------ find gmst --------------------------
        gmst = stu.gstime(jdut1)

        # find omeage from nutation theory
        omega = 125.04452222  + (-6962890.5390 *ttt + \
                7.455 *ttt*ttt + 0.008 *ttt*ttt*ttt)  / 3600.0
        omega = np.fmod(omega, 360.0) * deg2rad

        # teme does not include the geometric terms here
        # after 1997, kinematic terms apply
        if (jdut1 > 2450449.5) and (eqeterms > 0):
            gmstg = gmst \
                   + 0.00264*math.pi /(3600*180)*math.sin(omega) \
                   + 0.000063*math.pi /(3600*180)*math.sin(2.0 *omega)
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

        ateme = st@(pm@aecef + tt1 \
               + 2.0*tt2)

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


def ecef2pef  (recef, vecef, aecef, opt, xp, yp, ttt):

        pm = smu.polarm(xp, yp, ttt, opt)

        rpef = pm@recef

        vpef = pm@vecef

        apef = pm@aecef

        return rpef, vpef, apef



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
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
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


def eci2ecef  (reci, veci, aeci, ttt, jdut1, lod, xp, yp, eqeterms, ddpsi, ddeps):
        prec, psia, wa, ea, xa = obu.precess (ttt, '80')

        deltapsi, trueeps, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)

        st, stdot = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)

        pm = smu.polarm(xp, yp, ttt, '80')

        thetasa = earthrot * (1.0  - lod/86400.0)
#        omegaearth = np.array([[0], [0], [thetasa]])

        omegaearth = np.array([0.0, 0.0, thetasa])

        tmpmat = st.T@nut.T@prec.T

        rpef = tmpmat@reci
        recef = pm.T@rpef

        temp = np.cross(omegaearth, rpef.T) #turn from column to row vectors for cross product
        vpef = tmpmat@veci - temp.T
        vecef = pm.T@vpef

        # two additional terms not needed if satellite is not on surface
        # of the Earth
        c1 = np.cross(omegaearth, temp)
        c2 = 2.0*np.cross(omegaearth, vpef.T)
        aecef = np.matmul(pm.T, np.matmul(tmpmat, aeci)) - c1.T - c2.T
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
#    reci      - position vector eci          km
#    veci      - velocity vector eci          km/s
#    aeci      - acceleration vector eci      km/s2
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


def eci2mod  (reci, veci, aeci, ttt):

        prec, psia, wa, ea, xa = obu.precess (ttt, '80')

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
#    ttt         - julian centuries of tt         centuries
#    ddpsi       - correction for iau2000         rad
#    ddeps       - correction for iau2000         rad
#    ddx
#    ddy
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


def eci2tod  (reci, veci, aeci, opt, ttt, ddpsi, ddeps, ddx, ddy):


    prec, psia, wa, ea, xa = obu.precess (ttt, opt)

    if opt == '80':
        deltapsi, trueeps, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)
    else:
        # ---- ceo based, iau2006
        if opt == '6c':
            x, y, s, pnb = obu.iau06xys (ttt, ddx, ddy)

        # ---- class equinox based, 2000a
        if opt == '6a':
             deltapsi, pnb, prec, nut, l, l1, f, d, omega, \
                lonmer, lonven, lonear, lonmar, lonjup, lonsat, \
                    lonurn, lonnep, precrate \
                            = obu.iau06pna (ttt)

        # ---- class equinox based, 2000b
        if opt == '6b':
            deltapsi, pnb, prec, nut, l, l1, f, d, omega, \
                lonmer, lonven, lonear, lonmar, lonjup, lonsat, \
                    lonurn, lonnep, precrate \
                        = obu.iau06pnb (ttt)
        prec = np.eye(3)
        nut = pnb

    if sh.show == 'y':
        conv = math.pi / (180.0*3600.0)
        print('dpsi %11.7f trueeps %11.7f mean eps %11.7f deltaeps %11.7f \n'
              %(deltapsi/conv, trueeps/conv, meaneps/conv,
                (trueeps-meaneps)/conv))
        print('nut iau 76 \n')
        print('%20.14f %20.14f %20.14f \n'% (nut[0], nut[1], nut[2]))

    rtod = nut.T@prec.T@reci

    vtod = nut.T@prec.T@veci

    atod = nut.T@prec.T@aeci

    if sh.show == 'y':
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


def eci2teme  (reci, veci, aeci, ttt, ddpsi, ddeps):

    prec, psia, wa, ea, xa = obu.precess (ttt, '80')

    deltapsi, trueeps, meaneps, omega, nut = obu.nutation  (ttt, ddpsi, ddeps)

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


def gd2gc (latgd):
    # -------------------------  implementation   -----------------
    latgc = math.atan((1.0  - eccearthsqrd)*math.tan(latgd))
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


def gc2gd (latgc):
        # -------------------------  implementation   -----------------
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

def dms2rad(deg, min, sec):
        # ------------------------  implementation   ------------------
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

def rad2dms(dms):
        # ------------------------  implementation   ------------------
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

def radec2rv(rr, rtasc, decl, drr, drtasc, ddecl):

        # -------------------------  implementation   -----------------
        r = np.zeros((3))
        r[0] = rr*math.cos(decl)*math.cos(rtasc)
        r[1] = rr*math.cos(decl)*math.sin(rtasc)
        r[2] = rr*math.sin(decl)
        r = r.T

        v = np.zeros([3])
        v[0] = drr*math.cos(decl)*math.cos(rtasc) \
            - rr*math.sin(decl)*math.cos(rtasc)*ddecl \
            - rr*math.cos(decl)*math.sin(rtasc)*drtasc
        v[1] = drr*math.cos(decl)*math.sin(rtasc) \
            - rr*math.sin(decl)*math.sin(rtasc)*ddecl \
            + rr*math.cos(decl)*math.cos(rtasc)*drtasc
        v[2] = drr*math.sin(decl) + rr*math.cos(decl)*ddecl
        v = v.T

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


def rv2radec(r, v):
        # -------------------------  implementation   -------------------------
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


def rv2tradec (reci, veci, latgd, lon, alt, ttt, jdut1, lod, xp, yp,
               terms, ddpsi, ddeps):

    # --------------------- implementation ------------------------
    # ----------------- get site vector in ecef -------------------
    rsecef, vsecef = obu.site (latgd, lon, alt)

    #rs
    #vs
    # -------------------- convert ecef to eci --------------------
    a = np.zeros((3))
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
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
#    terms       - number of terms for ast calculation 0, 2
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


def rv2razel (reci, veci, latgd, lon, alt, ttt, jdut1, lod, xp, yp, terms,
              ddpsi, ddeps):
    # --------------------- implementation ------------------------
    # ----------------- get site vector in ecef -------------------
    rsecef, vsecef = obu.site(latgd, lon, alt)
    #print('rsecef    %14.7f %14.7f %14.7f \n', rsecef)

    # -------------------- convert eci to ecef --------------------
    a = np.array([[0],[0],[0]])
    recef, vecef, aecef = eci2ecef(reci, veci, a, ttt, jdut1, lod, xp, yp, terms,
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
    tempvec = smu.rot3(rhoecef, lon)
    rhosez = smu.rot2(tempvec, halfpi-latgd)

    tempvec = smu.rot3(drhoecef, lon)
    drhosez = smu.rot2(tempvec, halfpi-latgd)

    # ------------- calculate azimuth and elevation ---------------
    temp = np.sqrt(rhosez[0]*rhosez[0] + rhosez[1]*rhosez[1])
    if (temp < small):           # directly over the north pole
        el = np.sign(rhosez[2])*halfpi   # +- 90 deg
    else:
        magrhosez = smu.mag(rhosez)
        el = math.asin(rhosez[2] / magrhosez)

    if (temp < small):
        az = math.atan2(drhosez[1], -drhosez[0])
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
        delx = (drhosez[2] - drho*np.sin(el)) / temp
    else:
        delx = 0.0
    return rho, az, el, drho, daz, delx



#
# position and velocity to ecliptic latitude longitude
# dav 28 mar 04
#
# [rr, ecllon, ecllat, drr, decllon, decllat] = rv2ell (rijk, vijk)


def rv2ell (rijk, vijk):

        # --------------------  implementation   ----------------------
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



#
# ecliptic latitude longitude to position and velocity
# dav 28 mar 04
#
# [rijk, vijk] = ell2rv (rr, ecllon, ecllat, drr, decllon, decllat)


def ell2rv (rr, ecllon, ecllat, drr, decllon, decllat):

        # --------------------  implementation   ----------------------
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

        rijk = smu.rot1 (r, -obliquity)
        vijk = smu.rot1 (v, -obliquity)

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


def mod2eci(rmod, vmod, amod, ttt):

    prec, psia, wa, ea, xa = obu.precess(ttt, '80')

    reci = prec@rmod
    veci = prec@vmod
    aeci = prec@amod

    return reci.T, veci.T, aeci.T


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

def tod2eci (rtod, vtod, atod, ttt, ddpsi, ddeps):

    prec, psia, wa, ea, xa = obu.precess(ttt, '80')

    deltapsi, trueeps, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)

    reci = prec@nut@rtod

    veci = prec@nut@vtod

    aeci = prec@nut@atod

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
#    rmag        - eci position vector magnitude  km
#    vmag        - eci velocity vector magnitude  km/sec
#    latgc       - geocentric latitude            rad
#    lon         - longitude                      rad
#    fpa         - sat flight path angle          rad
#    az          - sat flight path az             rad
#    ttt         - julian centuries of tt         centuries
#    jdut1       - julian date of ut1             days from 4713 bc
#    lod         - excess length of day           sec
#    xp          - polar motion coefficient       arc sec
#    yp          - polar motion coefficient       arc sec
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
# [reci, veci] = flt2rv (rmag, vmag, latgc, lon, fpa, az, ttt, jdut1, lod, xp, yp, terms, ddpsi, ddeps)
# ----------------------------------------------------------------------------

def flt2rv(rmag=None, vmag=None, latgc=None, lon=None, fpa=None,
           az=None, ttt=None, jdut1=None, lod=None, xp=None,
           yp=None, terms=None, ddpsi=None, ddeps=None):
    small = 1e-08
    # -------- form position vector
    recef = np.zeros(3)
    recef[0] = rmag * np.cos(latgc) * np.cos(lon)
    recef[1] = rmag * np.cos(latgc) * np.sin(lon)
    recef[2] = rmag * np.sin(latgc)
    recef = np.transpose(recef)
    # -------- convert r to eci
    vecef = np.array([[0], [0], [0]])

    aecef = np.array([[0], [0], [0]])
    reci, veci, aeci = ecef2eci(recef, vecef, aecef, ttt,
                              jdut1, lod, xp, yp, terms, ddpsi, ddeps)
    # ------------- calculate rtasc and decl ------------------
    temp = np.sqrt(reci[0] * reci[0] + reci[1] * reci[1])
    if (temp < small):
        # v needs to be defined herexxxxxxxxx
        rtasc = math.atan2(veci[1], veci[0])
    else:
        rtasc = math.atan2(reci[1], reci[0])

    decl = np.arcsin(reci[2] / rmag)
    # -------- form velocity vector
    fpav = np.pi * 0.5 - fpa
    veci = np.zeros(3)
    veci[0] = vmag * (- np.cos(rtasc) * np.sin(decl)
                      * (np.cos(az) * np.cos(fpav) - np.sin(rtasc)
                         * np.sin(az) * np.cos(fpav))
                      + np.cos(rtasc) * np.sin(decl) * np.sin(fpav))
    veci[1] = vmag * (- np.sin(rtasc) * np.sin(decl)
                      * (np.cos(az) * np.cos(fpav) + np.cos(rtasc)
                         * np.sin(az) * np.cos(fpav))
                      + np.sin(rtasc) * np.cos(decl) * np.sin(fpav))
    veci[2] = vmag * (np.sin(decl) * np.sin(fpav)
                      + np.cos(decl) * np.cos(az) * np.cos(fpav))
    reci = np.transpose(reci)
    veci = np.transpose(veci)
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

def rvs2raz(rhosez=None, drhosez=None):
    # -------------------------  implementation   -----------------
    small = 1e-08
    # ------------- calculate azimuth and elevation ---------------
    temp = np.sqrt(rhosez[0] * rhosez[0] + rhosez[1] * rhosez[1])
    if (np.abs(rhosez[1]) < small):
        if (temp < small):
            az = math.atan2(drhosez[1], - drhosez[0])
        else:
            if (rhosez[0] > 0.0):
                az = np.pi
            else:
                az = 0.0
    else:
        az = math.atan2(rhosez[1], - rhosez[0])

    rho = smu.mag(rhosez)


    if ((temp < small)):
        el = np.sign(rhosez[2]) * halfpi
    else:
        el = np.arcsin(rhosez[2] / rho)

    # -------  calculate range, azimuth and elevation rates -------
    drho = np.dot(rhosez, drhosez.T) / rho
    if (np.abs(temp * temp) > small):
        daz = (drhosez[0] * rhosez[1] - drhosez[1] * rhosez[0]) / (temp * temp)
    else:
        daz = 0.0

    if (np.abs(temp) > small):
        delv = (drhosez[2] - drho * np.sin(el)) / temp
    else:
        delv = 0.0

    return rho, az, el, drho, daz, delv


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
    #    fprintf(' jd #16.8f lon #11.5f  incl #11.5f raan #11.5f argp #11.5f \n', jdut1, lon * rad2deg, incl * rad2deg, raan * rad2deg, argp * rad2deg)
#    need to use their GMST calculation
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

def pef2eci(rpef=None, vpef=None, apef=None, ttt=None, jdut1=None,
            lod=None, eqeterms=None, ddpsi=None, ddeps=None):
    prec, psia, wa, ea, xa = obu.precess(ttt, '80')
    deltapsi, trueeps, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)
    st, stdot = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)
    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    reci = prec @ nut @ st @ rpef
    veci = prec @ nut @ st @ (vpef + np.cross(omegaearth, rpef.T))
    temp = np.cross(omegaearth, rpef.T)
    aeci = prec @ nut @ st @ (apef + np.cross(omegaearth, temp) + 2.0 * np.cross(omegaearth, vpef.T))
    return reci, veci, aeci

# ------------------------------------------------------------------------------
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

def tradec2rv(rho=None, trtasc=None, tdecl=None, drho=None, dtrtasc=None,
              dtdecl=None, rseci=None, vseci=None, lod=None):
    # --------------------- implementation ------------------------
    latgc = np.arcsin(rseci(3) / smu.mag(rseci))
    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    np.cross(omegaearth, rseci, vseci)
    # --------  calculate topocentric slant range vectors ------------------
    rhov = np.zeros(3)
    drhov = np.zeros(3)
    rhov[0] = rho * np.cos(tdecl) * np.cos(trtasc)
    rhov[1] = rho * np.cos(tdecl) * np.sin(trtasc)
    rhov[2] = rho * np.sin(tdecl)
    drhov[0] = (drho * np.cos(tdecl) * np.cos(trtasc) - rho * np.sin(tdecl)
                * np.cos(trtasc) * dtdecl
                - rho * np.cos(tdecl) * np.sin(trtasc) * dtrtasc)
    drhov[1] = (drho * np.cos(tdecl) * np.sin(trtasc)
                - rho * np.sin(tdecl) * np.sin(trtasc) * dtdecl
                + rho * np.cos(tdecl) * np.cos(trtasc) * dtrtasc)
    drhov[2] = drho * np.sin(tdecl) + rho * np.cos(tdecl) * dtdecl
    # ------ find eci range vector from site to satellite ------
    reci = rhov + rseci
    veci = drhov + np.cos(latgc) * vseci
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

def teme2eci(rteme=None, vteme=None, ateme=None,
             ttt=None, ddpsi=None, ddeps=None):
    prec, psia, wa, ea, xa = obu.precess(ttt, '80')
    deltapsi, trueeps, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)
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

def teme2ecef(rteme=None, vteme=None, ateme=None, ttt=None, jdut1=None,
              lod=None, xp=None, yp=None, eqeterms=None):
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
        gmstg = (gmst + 0.00264 * np.pi / (3600 * 180) * np.sin(omega)
                 + 6.3e-05 * np.pi / (3600 * 180) * np.sin(2.0 * omega))
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

def eci2pef(reci=None, veci=None, aeci=None, opt=None, ttt=None,
            jdut1=None, lod=None, eqeterms=None, ddpsi=None,
            ddeps=None, ddx=None, ddy=None):
    prec, psia, wa, ea, xa = obu.precess(ttt, opt)

#need to fix this if else to default to '80' when opt is not specified

    if opt == '80' or not opt:
        deltapsi, trueeps, meaneps, omega, nut = obu.nutation(ttt, ddpsi, ddeps)
        st, stdot = stu.sidereal(jdut1, deltapsi, meaneps, omega, lod, eqeterms)
    else:
        # ---- ceo based, iau2006
        if opt == '6c':
            x, y, s, pnb = obu.iau06xys(ttt, ddx, ddy)
            st = obu.iau06era(jdut1)
        # ---- class equinox based, 2000a
        if opt == '6a':
            deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
                lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pna(ttt)
            gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')
        # ---- class equinox based, 2000b
        if opt == '6b':
            deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, lonear, \
                lonmar, lonjup, lonsat, lonurn, lonnep, precrate = obu.iau06pnb(ttt)
            gst, st = obu.iau06gst(jdut1, ttt, deltapsi, '06')
        prec = np.eye(3)
        nut = pnb
        st = st

    thetasa = earthrot * (1.0 - lod / 86400.0)
    omegaearth = np.array([0.0, 0.0, thetasa])
    rpef = st.T@nut.T@prec.T@reci
    vpef = st.T@nut.T@prec.T@veci - np.cross(omegaearth, rpef.T).T

    temp = np.cross(omegaearth, rpef.T)
    apef = st.T@nut.T@prec.T@aeci - np.cross(omegaearth, temp).T \
        - 2.0*np.cross(omegaearth, vpef.T).T

    return rpef, vpef, apef

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

def cirs2ecefiau06(rcirs=None, vcirs=None, acirs=None, ttt=None,
                   jdut1=None, lod=None, xp=None, yp=None, option=None):

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

    rpef = st.T@rcirs
    recef = pm.T@rpef

    vpef = st.T@vcirs - np.cross(omegaearth, rpef.T)
    vecef = pm.T@vpef

    temp = np.cross(omegaearth, rpef.T)
    aecef = pm.T@(st.T@acirs - np.cross(omegaearth, temp).T - 2.0*np.cross(omegaearth, vpef.T))

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

def cirs2eciiau06(rcirs=None, vcirs=None, acirs=None, ttt=None,
                  option=None, ddx=None, ddy=None):
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

def eci2cirsiau06(reci=None, veci=None, aeci=None, ttt=None,
                  option=None, ddx=None, ddy=None):

    # ---- ceo based, iau2000
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

def eci2ecefiau06(reci=None, veci=None, aeci=None, ttt=None,
                  jdut1=None, lod=None, xp=None, yp=None,
                  option=None, ddx=None, ddy=None):
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

    rpef = st.T@pnb.T@reci
    recef = pm.T@rpef

    print("spv")
    print(st.T@pnb.T@veci)
    print(np.cross(omegaearth, rpef.T))
    print("pmt")
    print(pm.T)
    vpef = st.T@pnb.T@veci - np.cross(omegaearth, rpef.T).T
    print("vpef")
    print(vpef)
    vecef = pm.T@vpef

    temp = np.cross(omegaearth, rpef.T)
    aecef = pm.T@(st.T@pnb.T@aeci - np.cross(omegaearth, temp).T
                  - 2.0*np.cross(omegaearth, vpef.T).T)

    #        if iauhelp == 'y'
    rtirs = pnb.T @ reci
    vtirs = pnb.T @ veci

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
        print(rtirs)

    if option == 'c':
        print('CIRS          IAU-2006 CIO ', (rtirs))

    print(' v ' , (vtirs))
    print('TIRS          IAU-2006 %c   ' % (option))
    print(rpef)
    print(' v ' , (vpef))
    #          end

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

def ecef2cirsiau06(recef=None, vecef=None, aecef=None, ttt=None,
                   jdut1=None, lod=None, xp=None, yp=None, option=None):
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

def ecef2eciiau06(recef=None, vecef=None, aecef=None, ttt=None,
                  jdut1=None, lod=None, xp=None, yp=None,
                  option=None, ddx=None, ddy=None):
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

def rad2hms(hms):
        # ------------------------  implementation   ------------------
        temp = 15.0 * math.pi/180.0

        temp = hms   / temp
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

def hms2rad(hr, min, sec):

        # ------------------------  implementation   ------------------
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

def twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst):

    # sgp4fix no longer needed, put in satrec
    # global tumin radiusearthkm xke j2 j3 j4 j3oj2

    # Needs a constant in space_constants?
    xpdotp   =  1440.0 / (2.0*math.pi)   # 229.1831180523293  # [rev/day]/[rad/min]

    year   = 0
    satrec = {}
    satrec['error'] = 0


    longstr1 = list(longstr1)
    longstr2 = list(longstr2)

   #     // set the implied decimal points since doing a formated read
    #     // fixes for bad input data values (missing, ...)
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


    #     // ---- find no, ndot, nddot ----
    satrec['no_kozai'] = satrec['no_kozai'] / xpdotp #//* rad/min
    satrec['nddot'] = satrec['nddot'] * 10.0**nexp
    # note the implied decimal is set when adjusting longstr1 above
    satrec['bstar'] = satrec['bstar'] * 10.0**ibexp

    #     // ---- convert to sgp4 units ----
    #    satrec['a']    = (satrec['no']*tumin)^(-2/3)                # [er]
    satrec['ndot'] = satrec['ndot'] / (xpdotp * 1440.0)          # [rad/min^2]
    satrec['nddot'] = satrec['nddot'] / (xpdotp * 1440.0 * 1440)     # [rad/min^3]

    #     // ---- find standard orbital elements ----
    satrec['inclo'] = satrec['inclo'] * deg2rad
    satrec['nodeo'] = satrec['nodeo'] * deg2rad
    satrec['argpo'] = satrec['argpo'] * deg2rad
    satrec['mo'] = satrec['mo'] * deg2rad

    #       // sgp4fix not needed here
    #    satrec['alta'] = satrec['a']*(1.0 + satrec['ecco']) - 1.0
    #    satrec['altp'] = satrec['a']*(1.0 - satrec['ecco']) - 1.0

    #     // ----------------------------------------------------------------
    #     // find sgp4epoch time of element set
    #     // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
    #     // and minutes from the epoch (time)
    #     // --------------------------------------------------------------

    #     // ------------- temp fix for years from 1957-2056 ----------------
    #     // ------ correct fix will occur when year is 4-digit in 2le ------
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

    #     // input start stop times manually
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



    #     // ------------- initialize the orbit at sgp4epoch --------------
    sgp4epoch = satrec['jdsatepoch'] + satrec['jdsatepochf'] - 2433281.5 # days since 0 Jan 1950
    satrec = obu.sgp4init(whichconst, opsmode, satrec,
                        satrec['jdsatepoch'] + satrec['jdsatepochf'] - 2433281.5,
                        satrec['bstar'], satrec['ndot'], satrec['nddot'],
                        satrec['ecco'], satrec['argpo'], satrec['inclo'],
                        satrec['mo'], satrec['no_kozai'],satrec['nodeo'])
    return startmfe, stopmfe, deltamin, satrec

##############################################################################################################
##############################################################################################################
##############################################################################################################

if __name__ == '__main__':


    #i dont think these are accurate!!!-jmb
    latgc = math.pi*0.75
    latgd = gc2gd (latgc)
    print(latgd, latgc)
    latgc = gd2gc (latgd)
    print(latgd, latgc)


    # LEO test
    recef = np.array([[-1033.4793830],  [7901.2952754],  [6380.3565958]])
    vecef = np.array([[-3.225636520],  [-2.872451450],   [5.531924446]])
    aecef = np.array([[0.001], [0.002], [0.003]])

    conv = math.pi / (180.0*3600.0)

    year = 2004
    mon = 4
    day = 6
    hr = 7
    min = 51
    sec = 28.386009

    dut1 = -0.4399619  # sec
    dat = 32         # sec
    xp = -0.140682 * conv  # " to rad
    yp = 0.333309 * conv
    lod = 0.0015563
    ddpsi = -0.052195 * conv  # " to rad
    ddeps = -0.003875 * conv
    ddx = -0.000205 * conv  # " to rad
    ddy = -0.000136 * conv
    order = 106
    eqeterms = 2
    timezone =0
    opt = 'c' # specify the iau00 cio approach

    print('input data \n\n')
    print(' year %5i '% year)
    print(' mon %4i '% mon)
    print(' day %3i '% day)
    print(' %3i:%2i:%8.6f\n'% (hr, min, sec))
    print(' dut1 %8.6f s'% dut1)
    print(' dat %3i s'% dat)
    print(' xp %8.6f "'% (xp / conv))
    print(' yp %8.6f "'% (yp / conv))
    print(' lod %8.6f s\n'% lod)
    print(' ddpsi %8.6f " ddeps  %8.6f\n'% (ddpsi/conv, ddeps/conv))
    print(' ddx   %8.6f " ddy    %8.6f\n'% (ddx/conv, ddy/conv))
    print(' order %3i  eqeterms %3i  opt %3s \n'% (order, eqeterms, opt))
    print('units are km and km/s and km/s2\n')

    # -------- convtime    - convert time from utc to all the others
    print('convtime results\n')
    ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, \
        ttdb, jdtdb, jdtdbfrac \
            = stu.convtime (year, mon, day, hr, min, sec, timezone, dut1, dat)
    print('ut1 %8.6f tut1 %16.12f jdut1 %18.11f '% (ut1, tut1, jdut1+jdut1frac))

    print('input vectors:')
    pp(recef)
    pp(vecef)
    pp(aecef)


#jdut1 2453101.82740678312


    nu = lon2nu (jdut1, 7.020438698 * deg2rad, 0.070273056 * deg2rad, 19.90450011 * deg2rad,
                 352.5056022 * deg2rad)
    print("nu is ", nu)
    lon = nu2lon(jdut1, nu, 0.070273056 * deg2rad, 19.90450011 * deg2rad, 352.5056022 * deg2rad)
    print("lon is ", lon)

    #--------------------------------ecef2--------------------------------------------------------------
    rpef, vpef, apef = ecef2pef(recef, vecef, aecef, '80', xp, yp, ttt)
    print('ecef2pef returned: ')
    pp(rpef)
    pp(vpef)
    pp(apef)

    reci, veci, aeci = pef2eci(rpef, vpef, apef, ttt, jdut1, lod, eqeterms, ddpsi, ddeps)
    print('pef2eci returned: ')
    pp(reci)
    pp(veci)
    pp(aeci)


    rpef, vpef, apef = eci2pef  (reci, veci, aeci, "6a", ttt, jdut1, lod,
                               eqeterms, ddpsi, ddeps)
    print('eci2pef returned: ')
    pp(rpef)
    pp(vpef)
    pp(apef)

    rpef, vpef, apef = eci2pef  (reci, veci, aeci, "6b", ttt, jdut1, lod,
                               eqeterms, ddpsi, ddeps)
    print('eci2pef returned: ')
    pp(rpef)
    pp(vpef)
    pp(apef)

    rpef, vpef, apef = eci2pef  (reci, veci, aeci, "6c", ttt, jdut1, lod,
                               eqeterms, ddpsi, ddeps, ddx, ddy)
    print('eci2pef returned: ')
    pp(rpef)
    pp(vpef)
    pp(apef)


    #this is not complete
    #rcirs, vcirs, acirs = ecef2cirsiau06(recef, vecef, aecef, ttt, jdut1, lod, xp, yp, 'b')
    #print('ecef2cirsiau06 returned: ')
    #pp(rcirs)
    #pp(vcirs)
    #pp(acirs)

    rcirs, vcirs, acirs = eci2cirsiau06  (reci, veci, aeci, ttt, 'a', ddx, ddy)
    print('eci2cirsiau06 a returned: ')
    pp(rcirs)
    pp(vcirs)
    pp(acirs)

    rcirs, vcirs, acirs = eci2cirsiau06  (reci, veci, aeci, ttt, 'b', ddx, ddy)
    print('eci2cirsiau06 b returned: ')
    pp(rcirs)
    pp(vcirs)
    pp(acirs)

    rcirs, vcirs, acirs = eci2cirsiau06  (reci, veci, aeci, ttt, 'c', ddx, ddy)
    print('eci2cirsiau06 c returned: ')
    pp(rcirs)
    pp(vcirs)
    pp(acirs)


    #this is not complete
    #recef, vecef, aecef = cirs2ecefiau06(rcirs, vcirs, acirs, ttt, jdut1, lod, xp, yp, 'a')
    #print('cirs2ecefiau06 returned: ')
    #pp(recef)
    #pp(vecef)
    #pp(aecef)

    #recef, vecef, aecef = cirs2ecefiau06(rcirs, vcirs, acirs, ttt, jdut1, lod, xp, yp, 'b')
    #print('cirs2ecefiau06 returned: ')
    #pp(recef)
    #pp(vecef)
    #pp(aecef)

    reci, veci, aeci = cirs2eciiau06 (rcirs, vcirs, acirs, ttt, 'a', ddx, ddy)
    print('cirs2eciiau06 returned: ')
    pp(reci)
    pp(veci)
    pp(aeci)

    reci, veci, aeci = cirs2eciiau06 (rcirs, vcirs, acirs, ttt, 'b', ddx, ddy)
    print('cirs2eciiau06 returned: ')
    pp(reci)
    pp(veci)
    pp(aeci)

    reci, veci, aeci = cirs2eciiau06 (rcirs, vcirs, acirs, ttt, 'c', ddx, ddy)
    print('cirs2eciiau06 returned: ')
    pp(reci)
    pp(veci)
    pp(aeci)

    reci, veci, aeci = ecef2eciiau06(recef, vecef, aecef, ttt, jdut1,
                                   lod, xp, yp, 'a', ddx, ddy)
    print('ecef2eciiau06 a returned: ')
    pp(reci)
    pp(veci)
    pp(aeci)

    reci, veci, aeci = ecef2eciiau06(recef, vecef, aecef, ttt, jdut1,
                                   lod, xp, yp, 'b', ddx, ddy)
    print('ecef2eciiau06 b returned: ')
    pp(reci)
    pp(veci)
    pp(aeci)

    reci, veci, aeci = ecef2eciiau06(recef, vecef, aecef, ttt, jdut1,
                                   lod, xp, yp, 'c', ddx, ddy)
    print('ecef2eciiau06 c returned: ')
    pp(reci)
    pp(veci)
    pp(aeci)

    recef, vecef, aecef = eci2ecefiau06(reci, veci, aeci, ttt, jdut1,
                                      lod, xp, yp, 'a', ddx, ddy)
    print('eci2ecefiau06 a returned: ')
    pp(recef)
    pp(vecef)
    pp(aecef)

    recef, vecef, aecef = eci2ecefiau06(reci, veci, aeci, ttt, jdut1,
                                      lod, xp, yp, 'b', ddx, ddy)
    print('eci2ecefiau06 b returned: ')
    pp(recef)
    pp(vecef)
    pp(aecef)

    recef, vecef, aecef = eci2ecefiau06(reci, veci, aeci, ttt, jdut1,
                                      lod, xp, yp, 'c', ddx, ddy)
    print('eci2ecefiau06 c returned: ')
    pp(recef)
    pp(vecef)
    pp(aecef)

    rteme, vteme, ateme = ecef2teme(recef, vecef, aecef, ttt,
                                    jdut1+jdut1frac, lod, xp, yp, eqeterms)
    print('ecef2teme returned: ')
    pp(rteme)
    pp(vteme)
    pp(ateme)

    rtod, vtod, atod = ecef2tod(recef, vecef, aecef, ttt,
                                jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps)
    print('ecef2tod returned: ')
    pp(rtod)
    pp(vtod)
    pp(atod)

    rmod, vmod, amod = ecef2mod(recef, vecef, aecef, ttt,
                                jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps)
    print('ecef2mod returned: ')
    pp(rmod)
    pp(vmod)
    pp(amod)

    recig, vecig, aecig = ecef2eci(recef, vecef, aecef, ttt,
                                   jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps)
    print('ecef2eci returned: ')
    pp(recig)
    pp(vecig)
    pp(aecig)

    #--------------------------------eci2--------------------------------------------------------------
    rteme, vteme, ateme = eci2teme(recig, vecig, aecig, ttt, ddpsi, ddeps)
    print('eci2teme returned:')
    pp(rteme)
    pp(vteme)
    pp(ateme)

    rtod, vtod, atod = eci2tod(recig, vecig, aecig, '80',
                               ttt, ddpsi, ddeps, ddx, ddy)
    print('eci2tod returned:')
    pp(rtod)
    pp(vtod)
    pp(atod)

    rmod, vmod, amod = eci2mod(recig, vecig, aecig, ttt)
    print('eci2mod returned:')
    pp(rmod)
    pp(vmod)
    pp(amod)

    recef, vecef, aecef = eci2ecef(recig, vecig, aecig, ttt,
                                   jdut1+jdut1frac, lod, xp, yp, 2, ddpsi, ddeps)
    print('eci2ecef returned:')
    pp(recef)
    pp(vecef)
    pp(aecef)


    reci, veci, aeci = tod2eci(rtod, vtod, atod, ttt, ddpsi, ddeps)
    print('tod2eci returned: ')
    pp(reci)
    pp(veci)
    pp(aeci)
    print("----------------------------------------")




    a = 6860.7631
    ecc = 0.75
    p = a*(1.0 - ecc**2)
    incl = 15 * deg2rad #97.65184/rad
    omega = 79.54701 * deg2rad
    argp = 83.86041 * deg2rad
    nu = 65.21303 * deg2rad
    arglat = 0.0
    truelon = 0.0
    lonper = 0.0

    print('p km =%f  a km =%f  ecc =%f  incl deg =%f  raan deg =%f argp deg =%f  nu deg =%f'% \
            (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg))
    print('incl rad =%f  raan rad =%f  argp rad =%f  nu rad =%f'% \
            (incl, omega, argp, nu))
    reci, veci = coe2rv (p, ecc, incl, omega, argp, nu, arglat, truelon, lonper)
    print("coe2rv returned: ", reci, veci)



#    reci = np.array([[11074.95274], [40629.74421], [-32.1123199]])
#    veci = np.array([[-2.940822436], [0.9007122363], [0.002036330819]])
#    reci = np.array([11074.95274, 40629.74421, -32.1123199])
#    veci = np.array([-2.940822436, 0.9007122363, 0.002036330819])
#    reci = np.array([6524.834000000,  6862.875000000, 6448.296000000])
#    veci = np.array([4.9013270000,    5.5337560000,   -1.9763410000])
    p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe (reci, veci)
    print('p km =%f  a km =%f  ecc =%f  incl =%f  raan =%f  argp =%f  nu =%f  m =%f '% \
            (p, a, ecc, incl, omega, argp, nu, m))
    print('p km =%f  a km =%f  ecc =%f  incl deg =%f  raan deg =%f  argp deg =%f  nu deg =%f  m deg =%f '% \
            (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg))
    print('     arglat   truelon    lonper ', \
            arglat, truelon, lonper)
#            arglat * rad2deg, truelon * rad2deg, lonper * rad2deg)





























