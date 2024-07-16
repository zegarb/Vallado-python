import math
import sys
import numpy as np
from numpy.polynomial import Polynomial
from pprint import pprint as pp
from space_constants import *
import orbit_utils as obu

#some debug stuff
from space_constants import sethelp as sh

# ------------------------------------------------------------------------------
#
#                                  rot1mat
#
#  this function sets up a rotation matrix for an input angle about the first
#    axis.
#
#  author        : david vallado                  719-573-2600   10 jan 2003
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    xval        - angle of rotation              rad
#
#  outputs       :
#    outmat      - matrix result
#
#  locals        :
#    c           - cosine of the angle xval
#    s           - sine of the angle xval
#
#  coupling      :
#    none.
#
# [outmat] = rot1mat (xval)
# ----------------------------------------------------------------------------- }

def rot1mat(xval: float):
    """this function sets up a rotation matrix for an input angle about the first
    axis.

    Parameters
    ----------
    xval : float
        angle of rotation: rad

    Returns
    -------
    outmat: ndarray
        rotation matrix
    """
    c = math.cos(xval)
    s = math.sin(xval)

    outmat = np.zeros((3, 3))
    outmat[0, 0] = 1.0
    outmat[0, 1] = 0.0
    outmat[0, 2] = 0.0
    outmat[1, 0] = 0.0
    outmat[1, 1] = c
    outmat[1, 2] = s
    outmat[2, 0] = 0.0
    outmat[2, 1] = -s
    outmat[2, 2] = c
    return outmat


# ------------------------------------------------------------------------------
#
#                                  rot2mat
#
#  this function sets up a rotation matrix for an input angle about the second
#    axis.
#
#  author        : david vallado                  719-573-2600   10 jan 2003
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    xval        - angle of rotation              rad
#
#  outputs       :
#    outmat      - matrix result
#
#  locals        :
#    c           - cosine of the angle xval
#    s           - sine of the angle xval
#
#  coupling      :
#    none.
#
# [outmat] = rot2mat (xval)
# ----------------------------------------------------------------------------- }


def rot2mat(xval: float):
    """this function sets up a rotation matrix for an input angle about the second
    axis.

    Parameters
    ----------
    xval : float
        angle of rotation: rad

    Returns
    -------
    outmat: ndarray
        rotation matrix
    """

    c = math.cos(xval)
    s = math.sin(xval)

    outmat = np.zeros((3, 3))
    outmat[0, 0] = c
    outmat[0, 1] = 0.0
    outmat[0, 2] = -s
    outmat[1, 0] = 0.0
    outmat[1, 1] = 1.0
    outmat[1, 2] = 0.0
    outmat[2, 0] = s
    outmat[2, 1] = 0.0
    outmat[2, 2] = c
    return outmat



# ------------------------------------------------------------------------------
#
#                                  rot3mat
#
#  this function sets up a rotation matrix for an input angle about the third
#    axis.
#
#  author        : david vallado                  719-573-2600   10 jan 2003
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    xval        - angle of rotation              rad
#
#  outputs       :
#    outmat      - matrix result
#
#  locals        :
#    c           - cosine of the angle xval
#    s           - sine of the angle xval
#
#  coupling      :
#    none.
#
# [outmat] = rot3mat (xval)
# -----------------------------------------------------------------------------


def rot3mat(xval: float):
    """this function sets up a rotation matrix for an input angle about the third
    axis.

    Parameters
    ----------
    xval : float
        angle of rotation: rad

    Returns
    -------
    outmat: ndarray
        rotation matrix
    """

    c = math.cos(xval)
    s = math.sin(xval)

    outmat = np.zeros((3, 3))
    outmat[0, 0] = c
    outmat[0, 1] = s
    outmat[0, 2] = 0.0
    outmat[1, 0] = -s
    outmat[1, 1] = c
    outmat[1, 2] = 0.0
    outmat[2, 0] = 0.0
    outmat[2, 1] = 0.0
    outmat[2, 2] = 1.0
    return outmat




# ------------------------------------------------------------------------------
#
#                           function rngaz
#
#  this function calculates the range and azimuth between two specified
#    ground points on a spherical earth.  notice the range will always be
#    within the range of values listed since you for not know the direction of
#    firing, long or short.  the function will calculate rotating earth ranges
#    if the tof is passed in other than 0.0 . range is calulated in rad and
#    converted to er by s = ro, but the radius of the earth = 1 er, so it's
#    s = o.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                - dav 31 may 06 add elliptical model
#
#  inputs          description                    range / units
#    llat        - start geocentric latitude      -pi/2 to  pi/2 rad
#    llon        - start longitude (west -)       0.0  to 2pi rad
#    tlat        - end geocentric latitude        -pi/2 to  pi/2 rad
#    tlon        - end longitude (west -)         0.0  to 2pi rad
#    tof         - time of flight if icbm, or 0.0 min
#
#  outputs       :
#    range       - range between points           km
#    az          - azimuth                        0.0  to 2pi rad
#
#  locals        :
#    none.
#
#  coupling      :
#    site, rot3, cross, atan2, dot, unit
#
#  references    :
#    vallado       2001, 774-775, eq 11-3, eq 11-4, eq 11-5
#
# [range, az] = rngaz (llat, llon, tlat, tlon, tof)
# rad = 180/pi
# [range, az] = rngaz(60.0/rad, -80.0/rad , 41.0/rad, 74.0/rad, 0)
# ------------------------------------------------------------------------------

def rngaz(llat: float, llon: float, tlat: float, tlon: float, tof: float):
    """this function calculates the range and azimuth between two specified
    ground points on a spherical earth.  notice the range will always be
    within the range of values listed since you for not know the direction of
    firing, long or short.  the function will calculate rotating earth ranges
    if the tof is passed in other than 0.0. range is calulated in rad and
    converted to er by s = ro, but the radius of the earth = 1 er, so it's
    s = o.

    Parameters
    ----------
    llat : float
        start geocentric latitude: -pi/2 to pi/2 rad
    llon : float
        start geocentric longitude (west -): -2pi to 2pi rad
    tlat : float
        end geocentric latitude: -pi/2 to pi/2 rad
    tlon : float
        end geocentric longitude (west -): -2pi to 2pi rad
    tof : float
        time of flight if icbm, or 0.0 min

    Returns
    -------
    range : float
        range between 2 points
    az: float
        azimuth: 0 to 2pi rad
    """
    small = smalle8
    omegaearth = omegaearthradptu
    # fix units on tof and omegaearth

    range_ = math.acos(math.sin(llat) * math.sin(tlat)
                       + math.cos(llat) * math.cos(tlat)
                       * math.cos(tlon - llon + omegaearth * tof))
    # ------ check if the range is 0 or half the earth  ---------
    if (abs(math.sin(range_) * math.cos(llat)) < small):
        if (abs(range_ - math.pi) < small):
            az = math.pi
        else:
            az = 0.0
    else:
        az = math.acos((math.sin(tlat) - math.cos(range_) * math.sin(llat))
                       / (math.sin(range_) * math.cos(llat)))

    # ------ check if the azimuth is grt than pi (180deg) -------
    if (math.sin(tlon - llon + omegaearth * tof) < 0.0):
        az = twopi - az

    ###looks like this was test code -jmb


    print('spehrical range %11.7f km az %11.7f \n'
          % (range_ * 6378.1363, az * 180 / math.pi))
    # test ellipsoidal approach
    alt = 0.0
    rlch, vlch = obu.site(llat, llon, alt)
    rtgt, vtgt = obu.site(tlat, tlon, alt)
    print('U rlch %11.7f %11.7f  %11.7f \n' % (rlch[0], rlch[1], rlch[2]))
    print('  rtgt %11.7f %11.7f  %11.7f \n' % (rtgt[0], rtgt[1], rtgt[2]))
    rlu = unit(rlch)

    rtu = unit(rtgt)
    print('u rlu %11.7f %11.7f  %11.7f \n' % (rlu[0], rlu[1], rlu[2]))
    print('  rtu %11.7f %11.7f  %11.7f \n' % (rtu[0], rtu[1], rtu[2]))
    rp = 6356.0

    delta = re * re / (rp ** 2) - 1.0
    eps = delta * (rlu[2] ** 2 + rtu[2] ** 2)
    w = np.cross(rlu, rtu)
    print('w UxV %11.7f %11.7f  %11.7f \n' % (w[0], w[1], w[2]))
    if (w[2] < 0.0):
        w = np.cross(rtu, rlu)
        print('w UxV %11.7f %11.7f  %11.7f \n' % (w[0], w[1], w[2]))

    v = np.cross(rlu, w)
    print('v uxw %11.7f %11.7f  %11.7f \n' % (v[0], v[1], v[2]))
    phi = math.pi - 0.5 * math.atan2(- 2 * v[2] * rlu[2], v[2] ** 2
                                     - rlu[2] ** 2)

    print('phi %11.7f %11.7f \n' % (phi, phi * rad2deg))
    # phi = 0.5*math.atan2(-2*.47024*.86603, .47024^2-.86603^2)
    # fprintf(1, 'phi #11.7f #11.7f \n', phi, phi * rad2deg)

    temp = np.array([np.dot(rlch, rlu), np.dot(rlch, rtu), 0.0])
    uprime = rot3(temp, phi) / 6378.137

    print('uprime %11.7f %11.7f  %11.7f \n'
          % (uprime[0], uprime[1], uprime[2]))
    phi1 = 2 * math.pi + math.atan2(uprime[1] * math.sqrt(1 + eps), uprime[0])
    print('phi1 %11.7f %11.7f \n' % (phi1, phi1 * rad2deg))
    temp = np.array([np.dot(rtgt, rlu), np.dot(rtgt, rtu), 0.0])
    vprime = rot3(temp, phi) / 6378.137

    print('vprime %11.7f %11.7f  %11.7f \n'
          % (vprime[0], vprime[1], vprime[2]))
    phi2 = 2 * math.pi + math.atan2(vprime[1] * math.sqrt(1 + eps), vprime[0])
    print('phi1 %11.7f %11.7f \n' % (phi2, phi2 * rad2deg))
    e = 0.08181922
    # do each half of integral evaluation
    phi = phi2
    m = 1
    r = 0
    s1 = ((math.factorial(2 * m) * (math.factorial(r)) ** 2)
          / (2 ** (2 * m - 2 * r) * math.factorial(2 * r + 1)
             * (math.factorial(m)) ** 2)
          * (math.cos(phi) ** (2 * r + 1)))
    f1 = math.comb(m, 2 * m) * phi / 2 ** (2 * m) + math.sin(phi) * s1
    r = 0
    m = 2
    s1 = ((math.factorial(2 * m) * (math.factorial(r)) ** 2)
          / (2 ** (2 * m - 2 * r) * math.factorial(2 * r + 1)
             * (math.factorial(m)) ** 2)
          * (math.cos(phi) ** (2 * r + 1)))
    r = 1
    s2 = ((math.factorial(2 * m) * (math.factorial(r)) ** 2)
          / (2 ** (2 * m - 2 * r) * math.factorial(2 * r + 1)
             * (math.factorial(m)) ** 2)
          * (math.cos(phi) ** (2 * r + 1)))
    f2 = math.comb(m, 2 * m) * phi / 2 ** (2 * m) + math.sin(phi) * (s1 + s2)
    funct2 = (1.0 - e ** 2 / 2 * f1 - math.comb(m - 2, 2 * m - 3) * (e ** 2 * f2)
              / (m * 2 ** (2 * m - 2)))
    phi = phi1
    m = 1
    r = 0
    s1 = ((math.factorial(2 * m) * (math.factorial(r)) ** 2)
          / (2 ** (2 * m - 2 * r) * math.factorial(2 * r + 1)
             * (math.factorial(m)) ** 2)
          * (math.cos(phi) ** (2 * r + 1)))
    f1 = math.comb(m, 2 * m) * phi / 2 ** (2 * m) + math.sin(phi) * s1
    r = 0
    m = 2
    s1 = ((math.factorial(2 * m) * (math.factorial(r)) ** 2)
          / (2 ** (2 * m - 2 * r) * math.factorial(2 * r + 1)
             * (math.factorial(m)) ** 2)
          * (math.cos(phi) ** (2 * r + 1)))
    r = 1
    s2 = ((math.factorial(2 * m) * (math.factorial(r)) ** 2)
          / (2 ** (2 * m - 2 * r) * math.factorial(2 * r + 1)
             * (math.factorial(m)) ** 2)
          * (math.cos(phi) ** (2 * r + 1)))
    f2 = math.comb(m, 2 * m) * phi / 2 ** (2 * m) + math.sin(phi) * (s1 + s2)
    funct1 = (1.0 - e ** 2 / 2 * f1
              - math.comb(m - 2, 2 * m - 3) * (e ** 2 * f2)
              / (m * 2 ** (2 * m - 2)))
    range_ = (funct2 - funct1) * re


    return range_, az

# ------------------------------------------------------------------------------
#
#                           function recovqt
#
#  this function recovers the time and function values in quartic blending routines.
#
#  author        : david vallado                  719-573-2600    11 dec 2002
#
#  revisions
#                - misc fixes                                      2 feb 2004
#
#  inputs          description                    range / units
#    p1, p2, p3, p4, p5, p6
#                - function values used for blending
#    root        - root used as variable
#
#  outputs       :
#    funvalue    - function value
#
#  locals        :
#    none
#
#  coupling      :
#    none
#
#  references    :
#    vallado       2001, 900
#
# [funvalue] = recovqt (p1, p2, p3, p4, p5, p6, root)
# ------------------------------------------------------------------------------

def recovqt(p1: float, p2: float, p3: float, p4: float, p5: float, p6: float,
            root: float):
    """this function recovers the time and function values in quartic blending
    routines.

    Parameters
    ----------
    p1, p2, p3, p4, p5, p6 : float
        function values used for blending
    root:  float
        root used as variable

    Returns
    -------
    funvalue: float
        function value
    """

    # ------ set up function from C-45 --------
    # aqit5*x**5 + aqit4*x**4 + etc
    temp = 1.0 / 24.0
    aqit0 = p3
    aqit1 = (2 * p1 - 16 * p2 + 16 * p4 - 2 * p5) * temp
    aqit2 = (- 1 * p1 + 16 * p2 - 30 * p3 + 16 * p4 - p5) * temp
    aqit3 = (- 9 * p1 + 39 * p2 - 70 * p3 + 66 * p4 - 33 * p5 + 7 * p6) * temp
    aqit4 = (13 * p1 - 64 * p2 + 126 * p3 - 124 * p4 + 61 * p5 - 12 * p6) * temp
    aqit5 = (- 5 * p1 + 25 * p2 - 50 * p3 + 50 * p4 - 25 * p5 + 5 * p6) * temp
    # ----- recover the variable value
    funvalue = (aqit5 * root ** 5 + aqit4 * root ** 4 + aqit3 * root ** 3
                + aqit2 * root ** 2 + aqit1 * root + aqit0)
    return funvalue

# ------------------------------------------------------------------------------
#
#                           function recovqd
#
#  this function recovers the time and function values in parabolic blending routines.
#
#  author        : david vallado                  719-573-2600     3 jan 2003
#
#  revisions
#                - misc fixes                                      2 feb 2004
#
#  inputs          description                    range / units
#    p1, p2, p3    - function values used for blending
#    root        - root used as variable
#
#  outputs       :
#    funvalue    - function value
#
#  locals        :
#    none
#
#  coupling      :
#    none
#
#  references    :
#    vallado       2013, 1032
#
# [funvalue] = recovqd (p1, p2, p3, root)
# ------------------------------------------------------------------------------


def recovqd(p1: float, p2: float, p3: float, root: float):
    """this function recovers the time and function values in parabolic
    blending routines.

    Parameters
    ----------
    p1, p2, p3 : float
        function values used for blending
    root : float
        root used as variable

    Returns
    -------
    funvalue: float
        function value
    """
    # ------ set up function from C-39 --------
    aqd0 = p1
    aqd1 = (- 3.0 * p1 + 4.0 * p2 - p3) * 0.5
    aqd2 = (p1 - 2.0 * p2 + p3) * 0.5
    # ----- recover the variable value
    funvalue = aqd2 * root ** 2 + aqd1 * root + aqd0
    return funvalue

# ------------------------------------------------------------------------------
#
#                           function recovpar
#
#  this function recovers the time and function values in parabolic blending routines.
#
#  author        : david vallado                  719-573-2600    11 dec 2002
#
#  revisions
#                - fix for single parabbln                        19 sep 2003
#                - misc fixes                                      2 feb 2004
#
#  inputs          description                    range / units
#    p1, p2, p3, p4 - function values used for blending
#    root        - root used as variable
#
#  outputs       :
#    funvalue    - function value
#
#  locals        :
#    none
#
#  coupling      :
#    none
#
#  references    :
#    vallado       2001, 897
#
# [funvalue] = recovpar (p1, p2, p3, p4, root)
# ------------------------------------------------------------------------------

def recovpar(p1: float, p2: float, p3: float, p4: float, root: float):
    """this function recovers the time and function values in parabolic
    blending routines.

    Parameters
    ----------
    p1, p2, p3, p4 : float
        function values used for blending
    root : float
        root used as variable

    Returns
    -------
    funvalue: float
        function value
    """
    # ------ set up function from C-39 -------
    #  acut3*x**3 + acut2*x**2 + etc
    acut0 = p2
    acut1 = (- p1 + p3) * 0.5
    acut2 = p1 - 2.5 * p2 + 2.0 * p3 - 0.5 * p4
    acut3 = - 0.5 * p1 + 1.5 * p2 - 1.5 * p3 + 0.5 * p4
    # ----- recover the variable value
    funvalue = acut3 * root ** 3 + acut2 * root ** 2 + acut1 * root + acut0
    return funvalue


# ------------------------------------------------------------------------------
#
#                           function parabbln
#
#  this function performs parabolic blending of an input zero crossing
#  function in order to find event times.
#
#  author        : david vallado                  719-573-2600    18 dec 2002
#
#  revisions
#                - fix eqt ref                                     3 jan 2003
#                - misc fixes                                      2 feb 2004
#
#  inputs          description                    range / units
#    p1, p2, p3    - function values used for blending
#
#  outputs       :
#    minfound    - test of success
#    rootf       - root for the function
#    funrate     - function rate
#
#  locals        :
#
#  coupling      :
#    quadric     - find roots of a quadric
#
#  references    :
#    vallado       2007, 979
#
# [minfound, rootf, funrate] = parabbln(p1, p2, p3)
# ------------------------------------------------------------------------------

### definitely not finished; minfound never changed from 'n' -zeg
def parabbln(p1: float, p2: float, p3: float):
    """this function performs parabolic blending of an input zero crossing
    function in order to find event times.

    Parameters
    ----------
    p1, p2, p3: float
        function values used for blending

    Returns
    -------
    minfound : str
        test of success: always 'n'?
    rootf : float
        root of the function
    funrate : float
        function rate
    """

    rootf = 0.0
    funrate = 0.0
    minfound = 'n'
    # ------ set up function from C-37 --------
    aqd0 = p1
    aqd1 = (- 3.0 * p1 + 4.0 * p2 - p3) * 0.5
    aqd2 = (p1 - 2.0 * p2 + p3) * 0.5
    # --------------- solve roots of this function -------------
    opt = 'U'
    r1r, r1i, r2r, r2i = quadric(aqd2, aqd1, aqd0, opt)
    # ---------- search through roots to locate answers --------
    for indx2 in range(3):
        if (indx2 == 0):
            root = r1r
        if (indx2 == 1):
            root = r2r
        if ((root >= 0.0) and (root <= 2.0)):
            # [time] = recovqd(t1, t2, t3, root) # should be 0.0!!!!!!
            ans = recovqd(p1, p2, p3, root)
            # ----- recover the function value derivative
            funrate = 2.0 * aqd2 * root + aqd1

    return minfound, rootf, funrate

# ------------------------------------------------------------------------------
#
#                           function quartbln
#
#  this function performs quartic blending of an input zero crossing
#  function in order to find event times.
#
#  author        : david vallado                  719-573-2600   18 dec 2002
#
#  revisions
#                - misc fixes                                      2 feb 2004
#
#  inputs          description                    range / units
#    p1, p2, p3, p4, p5, p6
#                - function values used for blending
#
#  outputs       :
#    minfound    - test of success
#    rootf       - root for the function
#    funrate     - function rate
#
#  locals        :
#
#  coupling      :
#    quintic     - find roots of a quintic
#
#  references    :
#    vallado       2001, 899-901
#
# [minfound, rootf, funrate] = quartbln (p1, p2, p3, p4, p5, p6)
# ------------------------------------------------------------------------------

### np.roots is old, switch to numpy Polynomial class instead? -zeg
def quartbln(p1: float, p2: float, p3: float, p4: float, p5: float,
             p6: float):
    """this function performs quartic blending of an input zero crossing
    function in order to find event times.

    Parameters
    ----------
    p1, p2, p3, p4, p5, p6: float
        function values used for blending

    Returns
    -------
    minfound: str
        test of success: 'y' or 'n'
    rootf: float
        root of the function
    funrate: float
        function rate
    """
    rootf = 0.0
    funrate = 0.0
    minfound = 'n'
    # ------ set up function from C-45 --------
#  aqit5*x**5 + aqit4*x**4 + etc
    temp = 1.0 / 24.0
    aqi0 = p3
    aqi1 = (2 * p1 - 16 * p2 + 16 * p4 - 2 * p5) * temp
    aqi2 = (- 1 * p1 + 16 * p2 - 30 * p3 + 16 * p4 - p5) * temp
    aqi3 = (- 9 * p1 + 39 * p2 - 70 * p3 + 66 * p4 - 33 * p5 + 7 * p6) * temp
    aqi4 = (13 * p1 - 64 * p2 + 126 * p3 - 124 * p4 + 61 * p5 - 12 * p6) * temp
    aqi5 = (- 5 * p1 + 25 * p2 - 50 * p3 + 50 * p4 - 25 * p5 + 5 * p6) * temp
    # --------------- solve roots of this function -------------
    opt = 'U'
    #       [r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i] = ...
#                   quintic(aqi5, aqi4, aqi3, aqi2, aqi1, aqi0, opt)

    # rt = np.roots(np.array([aqi5, aqi4, aqi3, aqi2, aqi1, aqi0]))
    rt = Polynomial([aqi0, aqi1, aqi2, aqi3, aqi4, aqi5]).roots()

    # ---------- search through roots to locate answers --------
    for indx2 in range(5):
        root = 99999.9
        if (indx2 == 0):
            root = rt[indx2]
        if (indx2 == 1):
            root = rt[indx2]
        if (indx2 == 2):
            root = rt[indx2]
        if (indx2 == 3):
            root = rt[indx2]
        if (indx2 == 4):
            root = rt[indx2]
        if ((root >= 0.0) and (root <= 1.0)):
            minfound = 'y'
            rootf = root
            #               [time] = recovqt(t1, t2, t3, t4, t5, t6, root)
            ans = recovqt(p1, p2, p3, p4, p5, p6, root)
            print("recovqt:". ans)
            # ----- recover the function value derivative
            funrate = (5.0 * aqi5 * root ** 4
                       + 4.0 * aqi4 * root ** 3
                       + 3.0 * aqi3 * root ** 2
                       + 2.0 * aqi2 * root + aqi1)

    return minfound, rootf, funrate

# ------------------------------------------------------------------------------
#
#                           function quartic
#
#  this function solves for the four roots of a quartic equation.  there are
#    no restrictions on the coefficients, and imaginary results are passed
#    out as separate values.  the general form is y = ax4 + bx3 + cx2 + dx + e.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  revisions
#    vallado     - convert to matlab              719-573-2600   18 dec 2002
#
#  inputs          description                    range / units
#    a           - coeficient of x fourth term
#    b           - coefficient of x cubed term
#    c           - coefficient of x squared term
#    d           - coefficient of x term
#    e           - constant
#    opt         - option for output              I all roots including imaginary
#                                                 R only real roots
#                                                 U only unique real roots (no repeated)
#
#  outputs       :
#    r1r         - real portion of root 1
#    r1i         - imaginary portion of root 1
#    r2r         - real portion of root 2
#    r2i         - imaginary portion of root 2
#    r3r         - real portion of root 3
#    r3i         - imaginary portion of root 3
#    r4r         - real portion of root 4
#    r4i         - imaginary portion of root 4
#
#  locals        :
#    temp1       - temporary value
#    temp2       - temporary value
#    s           - alternate variable
#    h           - temporary value
#    hsqr        - h squared
#    hcube       - h cubed
#    p           - term in auxillary equation
#    q           - term in auxillary equation
#    r           - term in auxillary equation
#    delta       - discriminator for use with cardans formula
#    e0          - angle holder for trigonometric solution
#    phi         - angle used in trigonometric solution
#    cosphi      - cosine of phi
#    sinphi      - sine of phi
#    rprime      - values of roots before final work
#    temp        - temporary variable in finding max rprime
#    eta         - constant coefficient in quadric solutions
#    beta        - constant coefficient in quadric solutions
#
#  coupling      :
#    quadric     find roots of a quadric polynomial
#
#  references    :
#    vallado       2007, 976
#
# [r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i] = quartic(a, b, c, d, e, opt)
# ------------------------------------------------------------------------------

### all of the polynomial functions could just be replaced with
### numpy's polynomial class... -zeg
def quartic(a=None, b=None, c=None, d=None, e=None, opt=None):
    # --------------------  implementation   ----------------------
    onethird = 1.0 / 3.0
    small = 1e-08
    r1r = 0.0
    r1i = 0.0
    r2r = 0.0
    r2i = 0.0
    r3r = 0.0
    r3i = 0.0
    r4r = 0.0
    r4i = 0.0
    root1 = 0.0
    root2 = 0.0
    root3 = 0.0
    if (np.abs(a) > small):
        # ----------- force coefficients into std form ----------------
        b = b / a
        c = c / a
        d = d / a
        e = e / a
        h = - b / 4
        hsqr = h ** 2
        hcube = hsqr * h
        p = 6.0 * hsqr + 3.0 * b * h + c
        q = 4.0 * hcube + 3.0 * b * hsqr + 2.0 * c * h + d
        r = h * hcube + b * hcube + c * hsqr + d * h + e
        a = onethird * (- p * p - 12.0 * r)
        b = (1.0 / 27.0) * (- 2.0 * p * p * p + 72.0 * p * r - 27.0 * q * q)
        s = - 2.0 * onethird * p
        delta = (a * a * a / 27.0) + (b * b * 0.25)
        if (np.abs(q) > small):
            # ------------------ use cardans formula ------------------
            if (delta > small):
                temp1 = (- b * 0.5) + np.sqrt(delta)
                temp2 = (- b * 0.5) - np.sqrt(delta)
                temp1 = np.sign(temp1) * np.abs(temp1) ** onethird
                temp2 = np.sign(temp2) * np.abs(temp2) ** onethird
                root1 = temp1 + temp2
                root2 = - 0.5 * (temp1 + temp2)
                r2i = - 0.5 * np.sqrt(3.0) * (temp1 - temp2)
                root3 = - 0.5 * (temp1 + temp2)
                r3i = - r2i
            else:
                # --------------- evaluate zero point -----------------
                if (np.abs(delta) < small):
                    root1 = - 2.0 * np.sign(b) * np.abs(b * 0.5) ** onethird
                    root2 = np.sign(b) * np.abs(b * 0.5) ** onethird
                    root3 = root2
                else:
                    # ------------ use trigonometric identities -------
                    e0 = 2.0 * np.sqrt(- a * onethird)
                    cosphi = (- b / (2.0 * np.sqrt(- a * a * a / 27.0)))
                    sinphi = np.sqrt(1.0 - cosphi * cosphi)
                    phi = math.atan2(sinphi, cosphi)
                    if (phi < 0.0):
                        phi = phi + 2 * np.pi
                    root1 = e0 * np.cos(phi * onethird)
                    root2 = e0 * np.cos(phi * onethird + 120.0 * deg2rad)
                    root3 = e0 * np.cos(phi * onethird + 240.0 * deg2rad)
            # --------------- find largest value of root -------------
            rprime = root1 + s
            if ((rprime < root2 + s) and (np.abs(r2i) < 0.0001)):
                rprime = root2 + s
            if ((rprime < root3 + s) and (np.abs(r3i) < 0.0001)):
                rprime = root3 + s
            # -- evaluate coefficients of two resulting quadratics ---
            if (rprime > small):
                eta = 0.5 * (p + rprime - q / np.sqrt(rprime))
                beta = 0.5 * (p + rprime + q / np.sqrt(rprime))
            else:
                eta = 0.5 * p
                beta = 0.5 * p
                if (rprime < 0.0):
                    rprime = - rprime
            r1r, r1i, r2r, r2i = quadric(1.0, np.sqrt(rprime), eta, opt)
            r3r, r3i, r4r, r4i = quadric(1.0, - np.sqrt(rprime), beta, opt)
        else:
            # ------ case where solution reduces to a quadratic ------
            r1r, r1i, r3r, r3i = quadric(1.0, p, r, opt)
            r = np.sqrt(r1r * r1r + r1i * r1i)
            phi = math.atan2(r1i, r1r)
            if (phi < 0.0):
                phi = phi + 2 * np.pi
            r1r = np.sqrt(r) * np.cos(phi * 0.5)
            r1i = np.sqrt(r) * np.sin(phi * 0.5)
            if (r1i > 0.0):
                r2r = r1r
            else:
                r2r = - r1r
            r2i = - r1i
            r = np.sqrt(r3r * r3r + r3i * r3i)
            phi = math.atan2(r3i, r3r)
            if (phi < 0.0):
                phi = phi + 2 * np.pi
            r3r = np.sqrt(r) * np.cos(phi * 0.5)
            r3i = np.sqrt(r) * np.sin(phi * 0.5)
            if (r3i > 0.0):
                r4r = r3r
            else:
                r4r = - r3r
            r4i = - r3i
        r1r = r1r + h
        r2r = r2r + h
        r3r = r3r + h
        r4r = r4r + h
        if (opt == 'R'):
            #            test
            pass
    else:
        r1r, r1i, r2r, r2i, r3r, r3i = cubic(b, c, d, e, opt)
        r4r = 99999.9
        r4i = 99999.9

    return r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i

# ------------------------------------------------------------------------------
#
#                           function quintic
#
#  this function solves for the five roots of a quintic equation.  there are
#    no restrictions on the coefficients, and imaginary results are passed
#    out as separate values.  the general form is y = ax5 + bx4 + cx3 + dx2 +
#    ex + f.
#
#     this routine uses a newton-raphson search to find
#     the real roots (x) between 0 and 1 of the polynomial
#        a*x**5+b*x**4+...
#     the resulting 4th order equation is solved by calling root4,
#     if needed. (repeated roots are ignored)
#     taken from "numerical analysis" by maron, copyright 1982, pg 80
#
#  author        : david vallado                  719-573-2600   16 dec 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    a           - coefficient of x quintic term
#    b           - coefficient of x quartic term
#    c           - coefficient of x cubic term
#    d           - coefficient of x squared term
#    e           - coefficient of x term
#    f           - constant
#    opt         - option for output              I all roots including imaginary
#                                                 R only real roots
#                                                 U only unique real roots (no repeated)
#
#  outputs       :
#    r1r         - real portion of root 1
#    r1i         - imaginary portion of root 1
#    r2r         - real portion of root 2
#    r2i         - imaginary portion of root 2
#    r3r         - real portion of root 3
#    r3i         - imaginary portion of root 3
#    r4r         - real portion of root 4
#    r4i         - imaginary portion of root 4
#    r5r         - real portion of root 5
#    r5i         - imaginary portion of root 5
#
#  locals        :
#    temp1       - temporary value
#    temp2       - temporary value
#    p           - coefficient of x squared term where x cubed term is 1.0
#    q           - coefficient of x term where x cubed term is 1.0
#    r           - coefficient of constant term where x cubed term is 1.0
#    delta       - discriminator for use with cardans formula
#    e0          - angle holder for trigonometric solution
#    phi         - angle used in trigonometric solution
#    cosphi      - cosine of phi
#    sinphi      - sine of phi
#     f        polynomial coefficients
#     aa(5)       intermediate polynomial coefficients
#     dz          change in intermediate variable
#     f           value of polynomial function at root z
#     flag        flag for printing error message
#     fp          slope of polynomial function at root z
#     i           counter
#     maxdz       max step size for finding z
#     n           polynomial order
#     numrts      integer number of real roots
#     olddz       previous value of dz
#     small         tolerance
#     p5        real roots
#     xx(5)       intermediate real roots
#     z           intermediate root
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2001,
#
# [r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i] = quintic (a, b, c, d, e, f, opt)
# ------------------------------------------------------------------------------

def quintic(a=None, b=None, c=None, d=None, e=None, f=None, opt=None):
    # --------------------  implementation   ----------------------
    opt = 'U'
    onethird = 1.0 / 3.0
    small = 1e-08
    temp = 1.0 / 24.0
    r1r = 0.0
    r1i = 0.0
    r2r = 0.0
    r2i = 0.0
    r3r = 0.0
    r3i = 0.0
    r4r = 0.0
    r4i = 0.0
    r5r = 0.0
    r5i = 0.0
    # -----   test to see if there cannot be a crossing on the interval 0 to 1
    aa = np.zeros(6)
    aa[4] = e
    aa[3] = aa[4] + d
    aa[2] = aa[3] + c
    aa[1] = aa[2] + b
    aa[0] = aa[1] + a
    if (((f > 0.0) and (np.amin(aa) > - f))):
        return r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i

    if (((f < 0.0) and (np.amax(aa) < - f))):
        return r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i

    # ---------- check to see if it really is fifth order ----------
    if (np.abs(a) < small):
        r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i = quartic(a, b, c, d, e, opt)
        return r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i

    if (np.abs(f) < small):
        r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i = quartic(b, c, d, e, f, opt)
        return r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i

    # ----- find a good first guess for a real root between 0 and 1
    # ----- by assigning 5 evenly spaced points on the interval
    p1 = f
    z = 0.25
    p2 = ((((a * z + b) * z + c) * z + d) * z + e) * z + f
    z = 0.5
    p3 = ((((a * z + b) * z + c) * z + d) * z + e) * z + f
    z = 0.75
    p4 = ((((a * z + b) * z + c) * z + d) * z + e) * z + f
    p5 = a + b + c + d + e + f
    # -----   find polynomial coefficients (assumes the interval 0 to 4)
    aa = p1
    ab = (- 50 * p1 + 96 * p2 - 72 * p3 + 32 * p4 - 6 * p5) * temp
    ac = (35 * p1 - 104 * p2 + 114 * p3 - 56 * p4 + 11 * p5) * temp
    ad = (- 10 * p1 + 36 * p2 - 48 * p3 + 28 * p4 - 6 * p5) * temp
    ae = (p1 - 4 * p2 + 6 * p3 - 4 * p4 + p5) * temp
    # -----   find the roots of the fourth order polynomial and bound them
    r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i = quartic(aa, ab, ac, ad, ae, opt)
    #       call bound4(xx, numrts, 0.0, 4.0)

    numrts = 0
    if ((np.abs(r1r) > small) and (np.abs(r1i) < small)):
        numrts = numrts + 1

    if ((np.abs(r2r) > small) and (np.abs(r2i) < small)):
        numrts = numrts + 1

    if ((np.abs(r3r) > small) and (np.abs(r3i) < small)):
        numrts = numrts + 1

    if ((np.abs(r4r) > small) and (np.abs(r4i) < small)):
        numrts = numrts + 1

    # -----   if no real roots found   return to main routine
    if (numrts < 1):
        return r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i

    # -----   initialize parameters using xx[0]/4 as a first guess for a real root
    # -----   (this rescales the root to the  interval 0 to 1)
    maxdz = 0.1
    olddz = maxdz
    z = r1r * 0.25
    for i in range(25):
        25
        f = ((((a * z + b) * z + c) * z + d) * z + e) * z + f
        if (np.abs(f) < small):
            # ----- find roots of resulting 4th order equation and return
            a = a
            b = b + z * aa[0]
            c = c + z * aa[1]
            d = d + z * aa[2]
            e = e + z * aa[3]
            r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i = quartic(a, b, c, d, e, opt)
            return r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i
        # ----- find new dz and z
        fp = (((5 * a * z + 4 * b) * z + 3 * c) * z + 2 * d) * z + e
        if (np.abs(fp) > small):
            dz = - f / fp
        if (np.abs(dz) > maxdz):
            dz = np.sign(maxdz) * dz
            if (olddz * dz < 0.0):
                maxdz = 0.8 * maxdz
        olddz = dz
        z = z + dz
        # ----- test for z out of limits, if so there are no roots
        if (((z > 1.0) or (z < 0.0))):
            numrts = 0
            return r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i


    # ---- once here, newton root search did not find an answer ----
    print(' warning: no convergence in quintic ' % ())
    return r1r, r1i, r2r, r2i, r3r, r3i, r4r, r4i, r5r, r5i


# ----------------------------------------------------------------------------
#
#                           function printdiff
#
#  this function prints a covariance matrix difference
#
#  author        : david vallado                  719-573-2600   23 may 2003
#
#  revisions
#
#  inputs          description                    range / units
#    strin       - title
#    mat1        - 6x6 input matrix
#    mat2        - 6x6 input matrix
#
#  outputs       :
#
#  locals        :
#
#  references    :
#    none
#
#  printdiff(strin, cov1, cov2)
# ----------------------------------------------------------------------------

def printdiff(strin: str, mat1: np.ndarray, mat2: np.ndarray):
    """this function prints a covariance matrix difference

    Parameters
    ----------
    strin : str
        title string
    mat1 : np.ndarray
        6x6 input matrix
    mat2 : np.ndarray
        6x6 input matrix
    """

    small = smalle18
    print('\ndiff %s' % (strin))
    print(mat1.T - mat2.T)
    print('\npctdiff %s pct over 1e-18  \n' % (strin))
    #    fprintf(1, '#14.4f#14.4f#14.4f#14.4f#14.4f#14.4f \n', 100.0*((mat1' - mat2')/mat1'))
    #    fprintf(1, 'Check consistency of both approaches tmct2cl-inv(tmcl2ct) diff pct over 1e-18 \n')
    #    fprintf(1, '-------- accuracy of tm comparing ct2cl and cl2ct --------- \n')
    tm1 = mat1.T
    tm2 = mat2.T
    diffmm = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            if (abs(tm1[i, j] - tm2[i, j]) < small) or (abs(tm1[i, j]) < small):
                diffmm[i, j] = 0.0
            else:
                diffmm[i, j] = 100.0 * ((tm1[i, j] - tm2[i, j]) / tm1[i, j])

    print('diffmm:\n', (diffmm))

# ------- two recursion algorithms needed by the lambertbattin routine
def seebatt(v: float):
    """a recursion algorithm used in the lambertbattin routine

    Parameters
    ----------
    v : float
        the value

    Returns
    -------
    seebat: float
        algorithm result
    """
    c = np.zeros(21)
    c[0] = 0.2
    c[1] = 9.0 / 35.0
    c[2] = 16.0 / 63.0
    c[3] = 25.0 / 99.0
    c[4] = 36.0 / 143.0
    c[5] = 49.0 / 195.0
    c[6] = 64.0 / 255.0
    c[7] = 81.0 / 323.0
    c[8] = 100.0 / 399.0
    c[9] = 121.0 / 483.0
    c[10] = 144.0 / 575.0
    c[11] = 169.0 / 675.0
    c[12] = 196.0 / 783.0
    c[13] = 225.0 / 899.0
    c[14] = 256.0 / 1023.0
    c[15] = 289.0 / 1155.0
    c[16] = 324.0 / 1295.0
    c[17] = 361.0 / 1443.0
    c[18] = 400.0 / 1599.0
    c[19] = 441.0 / 1763.0
    c[20] = 484.0 / 1935.0
    sqrtopv = np.sqrt(1.0 + v)
    eta = v / (1.0 + sqrtopv) ** 2
    # ------------------- process forwards ----------------------
    delold = 1.0
    termold = c[0]

    sum1 = termold
    i = 1
    while ((i <= 20) and (np.abs(termold) > 1e-08)):

        del_ = 1.0 / (1.0 + c[i] * eta * delold)
        term = termold * (del_ - 1.0)
        sum1 = sum1 + term
        i = i + 1
        delold = del_
        termold = term


    seebatt = 1.0 / ((1.0 / (8.0 * (1.0 + sqrtopv)))
                     * (3.0 + sum1 / (1.0 + eta * sum1)))
    c[0] = 9.0 / 7.0
    c[1] = 16.0 / 63.0
    c[2] = 25.0 / 99.0
    c[3] = 36.0 / 143.0
    c[4] = 49.0 / 195.0
    c[5] = 64.0 / 255.0
    c[6] = 81.0 / 323.0
    c[7] = 100.0 / 399.0
    c[8] = 121.0 / 483.0
    c[9] = 144.0 / 575.0
    c[10] = 169.0 / 675.0
    c[11] = 196.0 / 783.0
    c[12] = 225.0 / 899.0
    c[13] = 256.0 / 1023.0
    c[14] = 289.0 / 1155.0
    c[15] = 324.0 / 1295.0
    c[16] = 361.0 / 1443.0
    c[17] = 400.0 / 1599.0
    c[18] = 441.0 / 1763.0
    c[19] = 484.0 / 1935.0
    ktr = 20
    sum2 = 0.0
    term2 = 1.0 + c[ktr-1] * eta
    for i in range(ktr - 3):
        sum2 = c[ktr - i] * eta / term2
        term2 = 1.0 + sum2

    seebatt = (8.0 * (1.0 + sqrtopv)
               / (3.0 + (1.0 / (5.0 + eta + ((9.0 / 7.0) * eta / term2)))))
    return seebatt

# ------- two recursion algorithms needed by the lambertbattin routine
def kbat(v: float):
    """a recursion algorithm used by the lamertbattin routine

    Parameters
    ----------
    v : float
        _description_
    """
    d = np.zeros(21)
    d[0] = 1.0 / 3.0
    d[1] = 4.0 / 27.0
    d[2] = 8.0 / 27.0
    d[3] = 2.0 / 9.0
    d[4] = 22.0 / 81.0
    d[5] = 208.0 / 891.0
    d[6] = 340.0 / 1287.0
    d[7] = 418.0 / 1755.0
    d[8] = 598.0 / 2295.0
    d[9] = 700.0 / 2907.0
    d[10] = 928.0 / 3591.0
    d[11] = 1054.0 / 4347.0
    d[12] = 1330.0 / 5175.0
    d[13] = 1480.0 / 6075.0
    d[14] = 1804.0 / 7047.0
    d[15] = 1978.0 / 8091.0
    d[16] = 2350.0 / 9207.0
    d[17] = 2548.0 / 10395.0
    d[18] = 2968.0 / 11655.0
    d[19] = 3190.0 / 12987.0
    d[20] = 3658.0 / 14391.0
    # ----------------- process forwards ------------------------
    sum1 = d[0]
    delold = 1.0
    termold = d[0]
    i = 2
    ktr = 21
    while ((i < ktr) and (np.abs(termold) > 1e-08)):

        del_ = 1.0 / (1.0 + d[i] * v * delold)
        term = termold * (del_ - 1.0)
        sum1 = sum1 + term
        i = i + 1
        delold = del_
        termold = term


    sum2 = 0.0
    term2 = 1.0 + d[ktr-1] * v
    for i in range(ktr - 3):
        sum2 = d[ktr - i - 1] * v / term2
        term2 = 1.0 + sum2

    kbatt = d[0] / term2
    #            test = d[0] / ...
#                   (1 + (d[1]*v / ...
#                         (1 + (d(3)*v / ...
#                               (1 + (d(4)*v / ...
#                                     (1 + (d(5)*v / ...
#                                           (1 + (d(6)*v / ...
#                                                 (1 + (d(7)*v / ...
#                                                       (1 + (d(8)*v / ...
#                                                             (1 + (d(9)*v / ...
#                                                                  (1 + (d(10)*v / ...
#                                                                        (1 + (d(11)*v)))))))) ...
#                                                                       ))))))))))))
#kbatt = test
    return kbatt

# ------------------------------------------------------------------------------
#
#                           function cubicspl1
#
#  this function performs cubic splining of an input zero crossing
#  function in order to find event times.
#
#  author        : david vallado                  719-573-2600    18 dec 2002
#
#  revisions
#                - fix for single cubicspl                         8 dec 2003
#                - misc fixes                                      2 feb 2004
#
#  inputs          description                    range / units
#    p1, p2, p3, p4 - function values used for blending
#
#  outputs       :
#    minfound    - test of success
#    rootf       - root for the function
#    funrate     - function rate
#
#  locals        : none
#
#  coupling      :
#    cubic       - find roots of a cubic
#
#  references    :
#    vallado       2007, 981
#
# [minfound, rootf, funrate] = cubicspl1(p1, p2, p3, p4)
# ------------------------------------------------------------------------------

def cubicspl1(p1: float, p2: float, p3: float, p4: float):
    """this function performs cubic splining of an input zero crossing
    function in order to find event times.

    Parameters
    ----------
    p1, p2, p3, p4: float
        function values used for blending

    Returns
    -------
    minfound: str
        test of success: 'y' or 'n'
    rootf: float
        root of the function
    funrate: float
        function rate
    """
    rootf = 0.0
    funrate = 0.0
    minfound = 'n'
    # ---- check for true condition on first two points
    # if (indx == 1) & (sign(p1) ~= sign(p2))
    #     [evt1ctr, evt1, evt2ctr, evt2] = cubicbln (ev1n, ev2n, timearr, funarr, indx)
    #     event1ctr = event1ctr + evt1ctr
    #     event1 = [event1;evt1]
    #     event2ctr = event2ctr + evt2ctr
    #     event2 = [event2;evt2]

    # ------ set up function from C-39 --------
    acu0 = p2
    acu1 = (- p1 + p3) * 0.5
    acu2 = p1 - 2.5 * p2 + 2.0 * p3 - 0.5 * p4
    acu3 = - 0.5 * p1 + 1.5 * p2 - 1.5 * p3 + 0.5 * p4
    # --------------- solve roots of this function -------------
    opt = 'U'
    r1r, r1i, r2r, r2i, r3r, r3i = cubic(acu3, acu2, acu1, acu0, opt)
    # ---------- search through roots to locate answers --------
    for indx2 in range(3):
        if (indx2 == 0):
            root = r1r
        if (indx2 == 1):
            root = r2r
        if (indx2 == 2):
            root = r3r
        if (root >= 0.0) and (root <= 1.0):
            minfound = 'y'
            rootf = root
            #               [time] = recovpar(t1, t2, t3, t4, root)
            ans = recovpar(p1, p2, p3, p4, root)
            # ----- recover the function value derivative
            funrate = 3.0 * acu3 * root ** 2 + 2.0 * acu2 * root + acu1

    # fprintf(1, ' roots #11.7f  #11.7f  #11.7f \n', r1r, r2r, r3r)
# rt = roots([acu3 acu2 acu1 acu0])

    # ---- check for true condition on last two points
#       if (sign(p3) ~= sign(p4))
#           [evt1ctr, evt1, evt2ctr, evt2] = cubicbln (ev1n, ev2n, timearr, funarr, indx)
#           event1ctr = event1ctr + evt1ctr
#           event1 = [event1;evt1]
#           event2ctr = event2ctr + evt2ctr
#           event2 = [event2;evt2]
#         end

    return minfound, rootf, funrate


#
# ------------------------------------------------------------------------------
#
#                           function cubicspl
#
#  this function performs cubic splining of an input zero crossing
#  function in order to find function values.
#
#  author        : david vallado                  719-573-2600     2 feb 2004
#
#  revisions
#                -
#  inputs          description                    range / units
#    p1, p2, p3, p4 - function values used for splining
#    t1, t2, t3, t4 - time values used for splining
#
#  outputs       :
#    acu0..acu3  - splined polynomial coefficients. acu3 t^3, etc
#
#  locals        : none
#
#  coupling      :
#    cubic       - find roots of a cubic
#
#  references    :
#    vallado       2013, 559, 1034
#
# [acu0, acu1, acu2, acu3] = cubicspl (p1, p2, p3, p4)
# ------------------------------------------------------------------------------


def cubicspl(p1: float, p2: float, p3: float, p4: float):
    """this function performs cubic splining of an input zero crossing
    function in order to find function values.

    Parameters
    ----------
    p1, p2, p3, p4: float
        function values used for splining

    Returns
    -------
    acu0, acu1, acu2, acu3: float
        splined polynomial coefficients
    """
    # ------ set up function from C-41 --------
#       det = t1^3*t2^2 + t1^2*t2 + t1*t2^3 - t1^3*t2 - t1^2*t2^3 - t1*t2^2

    #       acu0 = p1
#       acu1 = ((t2^3-t2^2)*(p2-p1) + (t1^2-t1^3)*(p3-p1) + (t1^3*t2^2-t1^2*t2^3)*(p4-p1)) / det
#       acu2 = ((t2-t2^3)*(p2-p1) + (t1^3-t1)*(p3-p1) + (t1*t2^3-t1^3*t2)*(p4-p1)) / det
#       acu3 = ((t2^2-t2)*(p2-p1) + (t1-t1^2)*(p3-p1) + (t1^2*t2-t1*t2^2)*(p4-p1)) / det

    acu0 = p2
    acu1 = -p1 / 3.0 - 0.5 * p2 + p3 - p4 / 6.0
    acu2 = 0.5 * p1 - p2 + 0.5 * p3
    acu3 = -p1 / 6.0 + 0.5 * p2 - 0.5 * p3 + p4 / 6.0
    return acu0, acu1, acu2, acu3


# -----------------------------------------------------------------------------
#
#                           function cubicinterp
#
#  this function performs a cubic spline. four points are needed.
#
#  author        : david vallado                  719-573-2600   1 dec  2005
#
#  revisions
#
#  inputs          description                    range / units
#    valuein     - kp
#
#  outputs       :
#    out         - ap
#
#  locals        :
#                -
#
#  coupling      :
#    cubicspl
#
#  references    :
#    vallado       2013, 1027
# --------------------------------------------------------------------------- */


def cubicinterp(p1a: float, p1b: float, p1c: float, p1d: float, p2a: float,
                p2b: float, p2c: float, p2d: float, valuein: float):

    # double kc0, kc1, kc2, kc3, ac0, ac1, ac2, ac3,
    #        r1r, r1i, r2r, r2i, r3r, r3i, value

    # -------- assign function points ---------
    ac0, ac1, ac2, ac3 = cubicspl(p1a, p1b, p1c, p1d)
    kc0, kc1, kc2, kc3 = cubicspl(p2a, p2b, p2c, p2d)
    # recover the original function values
    # use the normalized time first, but at an arbitrary interval
    r1r, r1i, r2r, r2i, r3r, r3i = cubic(kc3, kc2, kc1, kc0 - valuein, 'R')
    #fprintf(1, 'cubic #11.7f  #11.7f  #11.7f  #11.7f  #11.7f \n', ac0, ac1, kc0, r1r, r2r)
    if ((r1r >= - 1e-06) and (r1r <= 1.001)):
        value = r1r
    else:
        if ((r2r >= - 1e-06) and (r2r <= 1.001)):
            value = r2r
        else:
            if ((r3r >= - 1e-06) and (r3r <= 1.001)):
                value = r3r
            else:
                value = 0.0
                print('error in cubicinterp root {0} {1} {2} {3} \n'
                      % (valuein, r1r, r2r, r3r))

    answer = ac3 * value ** 3 + ac2 * value * value + ac1 * value + ac0
    # cubicinterp
    return answer

# ------------------------------------------------------------------------------
#
#                           function unit
#
#  this function calculates a unit vector given the original vector.  if a
#    zero vector is input, the vector is set to zero.
#
#  author        : david vallado                  719-573-2600    4 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    vec         - vector
#
#  outputs       :
#    outvec      - unit vector
#
#  locals        :
#    i           - index
#
#  coupling      :
#    mag           magnitude of a vector
#
# [outvec] = norm (vec)
# ------------------------------------------------------------------------------


def unit(vec: np.ndarray):
    """this function calculates a unit vector given the original vector.  if a
    zero vector is input, the vector is set to zero.

    Parameters
    ----------
    vec : ndarray
        vector

    Returns
    -------
    outved: ndarray
        unit vector
    """
    small = smalle6
    magv = mag(vec)
    outvec = np.zeros(3)
    if (magv > small):
        for i in range(3):
            outvec[i] = vec[i] / magv
    else:
        for i in range(3):
            outvec[i] = 0.0

    return outvec



#
#  this rountine accomplishes the iteration work for the double-r angles
#  only routine
#
#
# dav 12-23-03
#

def doubler(cc1=None, cc2=None, magrsite1=None, magrsite2=None,
            magr1in=None, magr2in=None, los1=None, los2=None,
            los3=None, rsite1=None, rsite2=None, rsite3=None,
            t1=None, t3=None, direct=None, mu=None):

    rho1 = (- cc1 + np.sqrt(cc1 ** 2 - 4.0
                            * (magrsite1 ** 2 - magr1in ** 2))) / 2.0
    rho2 = (- cc2 + np.sqrt(cc2 ** 2 - 4.0
                            * (magrsite2 ** 2 - magr2in ** 2))) / 2.0

    r1 = rho1 * los1 + rsite1
    r2 = rho2 * los2 + rsite2
    magr1 = mag(r1)
    magr2 = mag(r2)

    if sh.show:
        print('start of loop  %11.7f  %11.7f  \n' % (magr1in, magr2in))
        print("r1:")
        print(r1)
        print("r2:")
        print(r2)

    if direct == 'y':
        w = np.cross(r1, r2) / (magr1 * magr2)
    else:
        w = - np.cross(r1, r2) / (magr1 * magr2)

    # change to negative sign
    rho3 = - np.dot(rsite3, w) / np.dot(los3, w)
    #rho1
    #rho2
    #rho3
    #los1
    #los2
    #los3
    r3 = np.multiply(rho3, los3) + rsite3
    magr3 = mag(r3)
    if sh.show:
        print('r3')
        print(r3)
        print('after 1st mag  %11.7f  %11.7f  %11.7f \n' % (magr1, magr2, magr3))
    cosdv21 = np.dot(r2, r1) / (magr2 * magr1)
    if math.isclose(magr2 * magr1, 0.0):
      print("Dont know what to do when we divide by zero here -jmb")
      return None, None, None, None, None, None, None, None, None
    else:
      sindv21 = mag(np.cross(r2, r1)) / (magr2 * magr1)
    dv21 = math.atan2(sindv21, cosdv21)
    cosdv31 = np.dot(r3, r1) / (magr3 * magr1)
    #    sindv31 = mag(np.cross(r3, r1))/(magr3*magr1)
    sindv31 = np.sqrt(1.0 - cosdv31 ** 2)
    dv31 = math.atan2(sindv31, cosdv31)
    cosdv32 = np.dot(r3, r2) / (magr3 * magr2)
    sindv32 = mag(np.cross(r3, r2)) / (magr3 * magr2)
    dv32 = math.atan2(sindv32, cosdv32)
    if dv31 > np.pi:
        c1 = (magr2 * sindv32) / (magr1 * sindv31)
        c3 = (magr2 * sindv21) / (magr3 * sindv31)
        p = (c1 * magr1 + c3 * magr3 - magr2) / (c1 + c3 - 1)
    else:
        c1 = (magr1 * sindv31) / (magr2 * sindv32)
        c3 = (magr1 * sindv21) / (magr3 * sindv32)
        p = (c3 * magr3 - c1 * magr2 + magr1) / (- c1 + c3 + 1)

    ecosv1 = p / magr1 - 1
    ecosv2 = p / magr2 - 1
    ecosv3 = p / magr3 - 1
    if dv21 != np.pi:
        esinv2 = (- cosdv21 * ecosv2 + ecosv1) / sindv21
    else:
        esinv2 = (cosdv32 * ecosv2 - ecosv3) / sindv32

    e = np.sqrt(ecosv2 ** 2 + esinv2 ** 2)
    a = p / (1 - e ** 2)
    if e * e < 0.99:
        n = np.sqrt(mu / a ** 3)
        s = magr2 / p * np.sqrt(1 - e ** 2) * esinv2
        c = magr2 / p * (e ** 2 + ecosv2)
        sinde32 = magr3 / np.sqrt(a * p) * sindv32 - magr3 / p * (1 - cosdv32) * s
        cosde32 = 1 - magr2 * magr3 / (a * p) * (1 - cosdv32)
        deltae32 = math.atan2(sinde32, cosde32)
        sinde21 = magr1 / np.sqrt(a * p) * sindv21 + magr1 / p * (1 - cosdv21) * s
        cosde21 = 1 - magr2 * magr1 / (a * p) * (1 - cosdv21)
        deltae21 = math.atan2(sinde21, cosde21)
        deltam32 = (deltae32 + 2 * s * (np.sin(deltae32 / 2)) ** 2
                    - c * np.sin(deltae32))
        deltam12 = (- deltae21 + 2 * s * (np.sin(deltae21 / 2)) ** 2
                    + c * np.sin(deltae21))
    else:
        print('hyperbolic, e1 is greater than 0.99 %11.7f \n' % (e))
        n = np.sqrt(mu / - a ** 3)
        s = magr2 / p * np.sqrt(e ** 2 - 1) * esinv2
        c = magr2 / p * (e ** 2 + ecosv2)
        sindh32 = (magr3 / np.sqrt(- a * p) * sindv32
                   - magr3 / p * (1 - cosdv32) * s)
        sindh21 = (magr1 / np.sqrt(- a * p) * sindv21
                   + magr1 / p * (1 - cosdv21) * s)
        deltah32 = np.log(sindh32 + np.sqrt(sindh32 ** 2 + 1))
        deltah21 = np.log(sindh21 + np.sqrt(sindh21 ** 2 + 1)) #s (and esinv2) being negative ripples down and makes this negative
        deltam32 = (- deltah32 + 2 * s * (np.sinh(deltah32 / 2)) ** 2
                    + c * np.sinh(deltah32))
        deltam12 = (deltah21 + 2 * s * (np.sinh(deltah21 / 2)) ** 2
                    - c * np.sinh(deltah21))
        # what if ends on hperbolic solution.
        # how to pass back deltae32?
        deltae32 = deltah32

    if sh.show:
        print('dm32 %11.7f  dm12 %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n'
          % (deltam32, deltam12, c1, c3, p, a, e, s, c))
        print('%11.7f %11.7f %11.7f \n' % (dv21, dv31, dv32))
    f1 = t1 - deltam12 / n
    f2 = t3 - deltam32 / n
    q1 = np.sqrt(f1 ** 2 + f2 ** 2)
    return r2, r3, f1, f2, q1, magr1, magr2, a, deltae32



### returns the cotangent of an angle in radians
def cot(arg):
    return  1.0/math.tan(arg)


# ------------------------------------------------------------------------------
#
#                           function findc2c3
#
#  this function calculates the c2 and c3 functions for use in the universal
#    variable calculation of z.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    znew        - z variable                     rad2
#
#  outputs       :
#    c2new       - c2 function value
#    c3new       - c3 function value
#
#  locals        :
#    sqrtz       - square root of znew
#
#  coupling      :
#    sinh        - hyperbolic sine
#    cosh        - hyperbolic cosine
#
#  references    :
#    vallado       2001, 70-71, alg 1
#
# [c2new, c3new] = findc2c3 (znew)
# ------------------------------------------------------------------------------


def findc2c3(znew: float):
    """this function calculates the c2 and c3 functions for use in the universal
    variable calculation of z.

    Parameters
    ----------
    znew : float
        z variable

    Returns
    -------
    c2new: float
        c2 function value
    c3new: float
        c3 function value
    """

    if (znew > small):
        sqrtz = math.sqrt(znew)
        c2new = (1.0 -math.cos(sqrtz)) / znew
        c3new = (sqrtz-math.sin(sqrtz)) / (sqrtz**3)
    elif (znew < -small):
        sqrtz = math.sqrt(-znew)
        c2new = (1.0 -math.cosh(sqrtz)) / znew
        c3new = (math.sinh(sqrtz) - sqrtz) / (sqrtz**3)
    else:
        c2new = 0.5
        c3new = 1.0 /6.0
    return c2new, c3new




# ------------------------------------------------------------------------------
#
#                           function newtonm
#
#  this function performs the newton rhapson iteration to find the
#    eccentric anomaly given the mean anomaly.  the true anomaly is also
#    calculated.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    ecc         - eccentricity                   0.0  to
#    m           - mean anomaly                   -2pi to 2pi rad
#
#  outputs       :
#    e0          - eccentric anomaly              0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#
#  locals        :
#    e1          - eccentric anomaly, next value  rad
#    sinv        - sine of nu
#    cosv        - cosine of nu
#    ktr         - index
#    r1r         - cubic roots - 1 to 3
#    r1i         - imaginary component
#    r2r         -
#    r2i         -
#    r3r         -
#    r3i         -
#    s           - variables for parabolic solution
#    w           - variables for parabolic solution
#
#  coupling      :
#    cubic       - solves a cubic polynomial
#
#  references    :
#    vallado       2001, 72-75, alg 2, ex 2-1
#
# [e0, nu] = newtonm (ecc, m)
# ------------------------------------------------------------------------------


def newtonm (ecc: float, m: float):
    """this function performs the newton rhapson iteration to find the
    eccentric anomaly given the mean anomaly. the true anomaly is also
    calculated.

    Parameters
    ----------
    ecc : float
        eccentricity
    m : float
        mean anomaly: -2pi to 2pi rad

    Returns
    -------
    e0 : float
        eccentric anomaly: 0 to 2pi rad
    nu : float
        true anomaly: 0 to 2pi rad
    """

    numiter = 50

    # -------------------------- hyperbolic  ----------------------
    if ((ecc-1.0) > small):
        # -------------------  initial guess -----------------------
        if (ecc < 1.6):
            if (((m < 0.0) and (m > -math.pi)) or (m > math.pi)):
                e0 = m - ecc
            else:
                e0 = m + ecc
        else:
            if ((ecc < 3.6) and (abs(m) > math.pi)):
                e0= m - np.sign(m)*ecc
            else:
                e0= m/(ecc-1.0)
        ktr = 1
        e1 = e0 + ((m-ecc*math.sinh(e0)+e0) / (ecc*math.cosh(e0) - 1.0))
        while ((abs(e1-e0)>small) and (ktr<= numiter)):
            e0= e1
            e1 = e0 + ((m - ecc*math.sinh(e0) + e0) / (ecc*math.cosh(e0) - 1.0))
            ktr = ktr + 1
        # ----------------  find true anomaly  --------------------
        sinv = -(math.sqrt(ecc*ecc-1.0) * math.sinh(e1)) / (1.0 - ecc*math.cosh(e1))
        cosv = (math.cosh(e1) - ecc) / (1.0 - ecc*math.cosh(e1))
        nu = math.atan2(sinv, cosv)

    # --------------------- parabolic -------------------------
    elif (abs(ecc-1.0) < small):
        #                c = [ 1.0/3.0 0.0; 1.0; -m]
        #                [r1r] = roots (c)
        #                e0= r1r
        s = 0.5  * (halfpi - math.atan(1.5 * m))
        w = math.atan(math.tan(s)**(1.0 / 3.0))
        e0 = 2.0 *cot(2.0 * w)
        ktr = 1
        nu = 2.0 * math.atan(e0)

    # -------------------- elliptical ----------------------
    elif (ecc > small):
        # -----------  initial guess -------------
        if (((m < 0.0) and (m > -math.pi)) or (m > math.pi)):
            e0 = m - ecc
        else:
            e0 = m + ecc
        #e0
        ktr = 1
        e1 = e0 + (m - e0 + ecc*math.sin(e0)) / (1.0 - ecc*math.cos(e0))
        while ((abs(e1-e0) > small) and (ktr <= numiter)):
            ktr = ktr + 1
            e0 = e1
            e1 = e0 + (m - e0 + ecc*math.sin(e0)) / (1.0 - ecc*math.cos(e0))
        e0 = e1
        # -------------  find true anomaly  ---------------
        sinv = ((math.sqrt(1.0 -ecc*ecc) * math.sin(e1))
                / (1.0 -ecc*math.cos(e1)))
        cosv = (math.cos(e1)-ecc) / (1.0 - ecc*math.cos(e1))
        nu = math.atan2(sinv, cosv)

    # -------------------- circular -------------------
    else:
        ktr = 0
        nu = m
        e0= m

    return e0 , nu




# ------------------------------------------------------------------------------
#
#                           function cubic
#
#  this function solves for the three roots of a cubic equation.  there are
#    no restrictions on the coefficients, and imaginary results are passed
#    out as separate values.  the general form is y = ax3 + bx2 + cx + d0.  note
#    that r1i will always be zero since there is always at least one real root.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  revisions
#    vallado     - convert to matlab              719-573-2600   18 dec 2002
#
#  inputs          description                    range / units
#    a3          - coefficient of x cubed term
#    b2          - coefficient of x squared term
#    c1          - coefficient of x term
#    d0          - constant
#    opt         - option for output              I all roots including imaginary
#                                                 R only real roots
#                                                 U only unique real roots (no repeated)
#
#  outputs       :
#    r1r         - real portion of root 1
#    r1i         - imaginary portion of root 1
#    r2r         - real portion of root 2
#    r2i         - imaginary portion of root 2
#    r3r         - real portion of root 3
#    r3i         - imaginary portion of root 3
#
#  locals        :
#    temp1       - temporary value
#    temp2       - temporary value
#    p           - coefficient of x squared term where x cubed term is 1.0
#    q           - coefficient of x term where x cubed term is 1.0
#    r           - coefficient of constant term where x cubed term is 1.0
#    delta       - discriminator for use with cardans formula
#    e0          - angle holder for trigonometric solution
#    phi         - angle used in trigonometric solution
#    cosphi      - cosine of phi
#    sinphi      - sine of phi
#
#  coupling      :
#    quadric     - roots of second order polynomial
#
#  references    :
#    vallado       2007, 975
#
# [r1r, r1i, r2r, r2i, r3r, r3i] = cubic (a3, b2, c1, d0, opt)
# ------------------------------------------------------------------------------

def cubic (a3, b2, c1, d0, opt):
    onethird = 1.0 / 3.0
    r1r = 0.0
    r1i = 0.0
    r2r = 0.0
    r2i = 0.0
    r3r = 0.0
    r3i = 0.0

    if (abs(a3) > small):
        # ----------- force coefficients into std form ----------------
        p = b2/a3
        q = c1/a3
        r = d0/a3

        a3 = onethird*(3.0 *q - p*p)
        b2 = (1.0 /27.0)*(2.0 *p*p*p - 9.0 *p*q + 27.0 *r)

        delta = (a3*a3*a3/27.0) + (b2*b2*0.25)

        # ------------------ use cardans formula ----------------------
        if (delta > small):
            temp1 = (-b2*0.5)+math.sqrt(delta)
            temp2 = (-b2*0.5)-math.sqrt(delta)
            temp1 = np.sign(temp1)*abs(temp1)**onethird
            temp2 = np.sign(temp2)*abs(temp2)**onethird
            r1r = temp1 + temp2 - p*onethird

            if (opt =='I'):
                r2r = -0.5 *(temp1 + temp2) - p*onethird
                r2i = -0.5 *math.sqrt(3.0)*(temp1 - temp2)
                r3r = -0.5 *(temp1 + temp2) - p*onethird
                r3i = -r2i
            else:
                r2r = 99999.9
                r3r = 99999.9
        else:
            # --------------- evaluate zero point ---------------------
            if (abs(delta) < small):
                r1r = -2.0*np.sign(b2)*abs(b2*0.5)**onethird - p*onethird
                r2r = np.sign(b2)*abs(b2*0.5)**onethird - p*onethird
                # if (opt =='U')
                    # r3r = 99999.9
                # else
                r3r = r2r

            else:
                # ------------ use trigonometric identities -----------
                e0 = 2.0 *math.sqrt(-a3*onethird)
                cosphi = (-b2/(2.0 *math.sqrt(-a3*a3*a3/27.0)))
                sinphi = math.sqrt(1.0 -cosphi*cosphi)
                phi = math.atan2(sinphi, cosphi)
                if (phi < 0.0):
                    phi = phi + 2.0*math.pi
                r1r = e0*math.cos(phi*onethird) - p*onethird
                r2r = e0*math.cos(phi*onethird + 120.0 * deg2rad) - p*onethird
                r3r = e0*math.cos(phi*onethird + 240.0 * deg2rad) - p*onethird
    else:
        r1r, r1i, r2r, r2i = quadric(b2, c1, d0, opt)
        r3r = 99999.9
        r3i = 99999.9
    return r1r, r1i, r2r, r2i, r3r, r3i




# ------------------------------------------------------------------------------
#
#                           function quadric
#
#  this function solves for the two roots of a quadric equation.  there are
#    no restrictions on the coefficients, and imaginary results are passed
#    out as separate values.  the general form is y = ax2 + bx + c.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  revisions
#    vallado     - convert to matlab              719-573-2600    3 dec 2002
#
#  inputs          description                    range / units
#    a           - coefficient of x squared term
#    b           - coefficient of x term
#    c           - constant
#    opt         - option for output              I all roots including imaginary
#                                                 R only real roots
#                                                 U only unique real roots (no repeated)
#
#  outputs       :
#    r1r         - real portion of root 1
#    r1i         - imaginary portion of root 1
#    r2r         - real portion of root 2
#    r2i         - imaginary portion of root 2
#
#  locals        :
#    discrim     - discriminate b2 - 4ac
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 974
#
# [r1r, r1i, r2r, r2i] = quadric   (a, b, c, opt)
# ------------------------------------------------------------------------------


def quadric(a, b, c, opt):

    small = 0.00000001
    r1r = 0.0
    r1i = 0.0
    r2r = 0.0
    r2i = 0.0

    discrim = b*b - 4.0 *a*c
    # ---------------------  real roots  --------------------------
    if (abs(discrim) < small):
        r1r = -b / (2.0 *a)
        r2r = r1r
        # if (opt =='U')
            # r2r = 99999.9

    elif abs(a) < small:
        r1r = -c/b
    elif (discrim > 0.0):
        r1r = (-b + math.sqrt(discrim)) / (2.0 *a)
        r2r = (-b - math.sqrt(discrim)) / (2.0 *a)
    else:
        # ------------------ complex roots --------------------
        if (opt =='I'):
            r1r = -b / (2.0 *a)
            r2r = r1r
            r1i = math.sqrt(-discrim) / (2.0 *a)
            r2i = -math.sqrt(-discrim) / (2.0 *a)
        else:
            r1r = 99999.9
            r2r = 99999.9
    return r1r, r1i, r2r, r2i


# ------------------------------------------------------------------------------
#
#                           function newtonnu
#
#  this function solves keplers equation when the true anomaly is known.
#    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
#    the parabolic limit at 168 is arbitrary. the hyperbolic anomaly is also
#    limited. the hyperbolic sine is used because it's not double valued.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#    vallado     - fix small                                     24 sep 2002
#
#  inputs          description                    range / units
#    ecc         - eccentricity                   0.0  to
#    nu          - true anomaly                   -2pi to 2pi rad
#
#  outputs       :
#    e0          - eccentric anomaly              0.0  to 2pi rad       153.02 deg
#    m           - mean anomaly                   0.0  to 2pi rad       151.7425 deg
#
#  locals        :
#    e1          - eccentric anomaly, next value  rad
#    sine        - sine of e
#    cose        - cosine of e
#    ktr         - index
#
#  coupling      :
#    arcsinh     - arc hyperbolic sine
#    sinh        - hyperbolic sine
#
#  references    :
#    vallado       2007, 85, alg 5
#
# [e0, m] = newtonnu (ecc, nu)
# ------------------------------------------------------------------------------


def newtonnu (ecc: float, nu: float):
    """this function solves keplers equation when the true anomaly is known.
    the mean and eccentric, parabolic, or hyperbolic anomaly is also found.
    the parabolic limit at 168 is arbitrary. the hyperbolic anomaly is also
    limited. the hyperbolic sine is used because it's not double valued.

    Parameters
    ----------
    ecc : float
        eccentricity
    nu : float
        true anomaly: -2pi to 2pi rad

    Returns
    -------
    e0 : float
        eccentric anomaly: 0 to 2pi rad
    m : float
        mean anomaly: 0 to 2pi rad
    """

    e0 = undefined
    m = undefined
    small = 0.00000001
    # --------------------------- circular ------------------------
    if (abs(ecc) < small):
        m = -0.0
        e0 = -0.0
        #print("m is ", m)

    # ---------------------- elliptical -----------------------
    elif (ecc < 1.0 - small):
        sine = ((math.sqrt(1.0 -ecc*ecc) * math.sin(nu))
                / (1.0 +ecc*math.cos(nu)))
        cose = (ecc + math.cos(nu)) / (1.0  + ecc*math.cos(nu))
        e0 = math.atan2(sine, cose)
        m = e0 - ecc*math.sin(e0)
        #print("m here is ", m)

    # -------------------- hyperbolic  --------------------
    elif (ecc > 1.0 + small):
        nu_check = nu
        if nu_check > math.pi:
            nu_check = nu_check - twopi
        if (abs(nu_check)+0.00001 < math.pi-math.acos(1.0 /ecc)):
            sine = ((math.sqrt(ecc*ecc-1.0) * math.sin(nu))
                    / (1.0  + ecc*math.cos(nu)))
            e0 = math.asinh(sine)
            m = ecc*math.sinh(e0) - e0
            #print("m now is ", m)

    # ----------------- parabolic ---------------------
    else:
        nu_check = nu
        if nu_check > math.pi:
            nu_check = nu_check - twopi
        if (abs(nu_check) < 168.0*math.pi/180.0):
            e0= math.tan(nu*0.5)
            m = e0 + (e0*e0*e0)/3.0
            #print("m is ", m)

    if (ecc < 1.0):
        m = np.fmod(m, 2.0 * math.pi)
        if (m < 0.0):
            m = m + 2.0 * math.pi
        e0 = np.fmod(e0, 4.0 *math.pi)

    return e0, m

# ------------------------------------------------------------------------------
#
#                           function newtone
#
#  this function solves keplers equation when the eccentric, paraboic, or
#    hyperbolic anomalies are known. the mean anomaly and true anomaly are
#    calculated.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    ecc         - eccentricity                   0.0  to
#    e0          - eccentric anomaly              -2pi to 2pi rad
#
#  outputs       :
#    m           - mean anomaly                   0.0  to 2pi rad
#    nu          - true anomaly                   0.0  to 2pi rad
#
#  locals        :
#    sinv        - sine of nu
#    cosv        - cosine of nu
#
#  coupling      :
#    sinh        - hyperbolic sine
#    cosh        - hyperbolic cosine
#
#  references    :
#    vallado       2001, 85, alg 6
#
# [m, nu] = newtone (ecc, e0)
# ------------------------------------------------------------------------------


def newtone (ecc: float, e0: float):
    """this function solves keplers equation when the eccentric, paraboic, or
    hyperbolic anomalies are known. the mean anomaly and true anomaly are
    calculated.

    Parameters
    ----------
    ecc : float
        eccentricity
    e0 : float
        eccentric anomaly: -2pi to 2pi rad

    Returns
    -------
    m : float
        mean anomaly: 0 to 2pi rad
    nu : float
        true anomaly: 0 to 2pi rad
    """

    small = 0.00000001

    # ------------------------- circular --------------------------
    if (abs(ecc) < small):
        m = e0
        nu = e0

    # ----------------------- elliptical ----------------------
    elif (ecc < 0.999):
        m = e0 - ecc*math.sin(e0)
        sinv = ((math.sqrt(1.0 -ecc*ecc) * math.sin(e0))
                / (1.0 -ecc*math.cos(e0)))
        cosv = (math.cos(e0)-ecc) / (1.0  - ecc*math.cos(e0))
        nu = math.atan2(sinv, cosv)

    # ---------------------- hyperbolic  ------------------
    elif (ecc > 1.0001):
        m = ecc*math.sinh(e0) - e0
        sinv = ((math.sqrt(ecc*ecc-1.0) * math.sinh(e0))
                / (1.0  - ecc*math.cosh(e0)))
        cosv = (math.cosh(e0)-ecc) / (1.0  - ecc*math.cosh(e0))
        nu = math.atan2(sinv, cosv)

    # -------------------- parabolic ------------------
    else:
        m = e0 + (1.0 /3.0)*e0*e0*e0
        #nu = 2.0 *datan(e0)  ###what is this? -todo
        nu = 2.0 * math.atan(e0)  ###what is this? -todo

    return m, nu

# ------------------------------------------------------------------------------
#
#                            function angl
#
#  this function calculates the angle between two vectors.  the output is
#    set to nan to indicate an undefined value.  be sure to check for
#    this at the output phase.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#    vallado     - fix tolerances                                 5 sep 2002
#
#  inputs          description                    range / units
#    vec1        - vector number 1
#    vec2        - vector number 2
#
#  outputs       :
#    theta       - angle between the two vectors  -pi to pi
#
#  locals        :
#    temp        - temporary real variable
#
#  coupling      :
#
# [theta] = angl (vec1, vec2)
# ----------------------------------------------------------------------------- }


def angl(vec1: np.ndarray, vec2: np.ndarray):
    """this function calculates the angle between two vectors.

    Parameters
    ----------
    vec1 : ndarray
        vector 1
    vec2 : ndarray
        vector 2

    Returns
    -------
    theta : float
        angle between two vetors: -pi to pi, nan if undefined
    """

    small = smalle8
    magv1 = mag(vec1)
    magv2 = mag(vec2)

    if magv1*magv2 > small**2:
        temp = np.dot(vec1, vec2) / (magv1*magv2)
        if abs(temp) > 1.0:
            temp = np.sign(temp) * 1.0
        theta = math.acos(temp)
    else:
        theta=undefined
    return theta

# ------------------------------------------------------------------------------
#
#                                  rot1
#
#  this function performs a rotation about the 1st axis.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    vec         - input vector
#    xval        - angle of rotation              rad
#
#  outputs       :
#    outvec      - vector result
#
#  locals        :
#    c           - cosine of the angle xval
#    s           - sine of the angle xval
#    temp        - temporary extended value
#
#  coupling      :
#    none.
#
# [outvec] = rot1 (vec, xval)
# ----------------------------------------------------------------------------- }

def rot1 (vec: np.ndarray, xval: float):
    """this function performs a rotation about the 1st axis.

    Parameters
    ----------
    vec : ndarray
        input vector
    xval : float
        angle of rotation: rad

    Returns
    -------
    outvec: ndarray
        result vector
    """
    temp = vec[2]
    c = math.cos(xval)
    s = math.sin(xval)

    outvec = np.zeros(3)
    outvec[2] = c*vec[2] - s*vec[1]
    outvec[1] = c*vec[1] + s*temp
    outvec[0] = vec[0]
    return outvec



# ------------------------------------------------------------------------------
#
#                            function rot2
#
#  this function performs a rotation about the 2nd axis.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    vec         - input vector
#    xval        - angle of rotation              rad
#
#  outputs       :
#    outvec      - vector result
#
#  locals        :
#    c           - cosine of the angle xval
#    s           - sine of the angle xval
#    temp        - temporary extended value
#
#  coupling      :
#    none.
#
# [outvec] = rot2 (vec, xval)
# ----------------------------------------------------------------------------- }

def rot2(vec, xval):
    """this function performs a rotation about the 2nd axis.

    Parameters
    ----------
    vec : ndarray
        input vector
    xval : float
        angle of rotation: rad

    Returns
    -------
    outvec: ndarray
        result vector
    """
    temp = vec[2]
    c = math.cos(xval)
    s = math.sin(xval)

    outvec = np.zeros(3)
    outvec[2] = c*vec[2] + s*vec[0]
    outvec[0] = c*vec[0] - s*temp
    outvec[1] = vec[1]
    return outvec

# ------------------------------------------------------------------------------
#
#                            function rot3
#
#  this function performs a rotation about the 3rd axis.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    vec         - input vector
#    xval        - angle of rotation              rad
#
#  outputs       :
#    outvec      - vector result
#
#  locals        :
#    c           - cosine of the angle xval
#    s           - sine of the angle xval
#    temp        - temporary extended value
#
#  coupling      :
#    none.
#
# [outvec] = rot3 (vec, xval)
# ----------------------------------------------------------------------------- }

def rot3 (vec, xval):
    """this function performs a rotation about the 3rd axis.

    Parameters
    ----------
    vec : ndarray
        input vector
    xval : float
        angle of rotation: rad

    Returns
    -------
    outvec: ndarray
        result vector
    """

    temp = vec[1]
    c = math.cos(xval)
    s = math.sin(xval)

    outvec = np.zeros(3)
    outvec[1] = c*vec[1] - s*vec[0]
    outvec[0] = c*vec[0] + s*temp
    outvec[2] = vec[2]
    return outvec

# ------------------------------------------------------------------------------
#
#                            function mag
#
#  this function finds the magnitude of a vector.  the tolerance is set to
#    0.000001, thus the 1.0e-12 for the squared test of underflows.
#
#  author        : david vallado                  719-573-2600   30 may 2002
#
#  revisions
#    vallado     - fix tolerance to match coe, eq, etc            3 sep 2002
#
#  inputs          description                    range / units
#    vec         - vector
#
#  outputs       :
#    mag         - magnitude
#
#  locals        :
#    none.
#
#  coupling      :
#    none.
#
# mag = (vec)
# ----------------------------------------------------------------------------- }

# Note: numpy has a "norm" function to take the magnitude of a vector
# mag = np.linalg.norm(x)

def mag(vec: np.ndarray):
    """this function finds the magnitude of a vector. the tolerance is set to
    0.000001, thus the 1.0e-12 for the squared test of underflows.

    Parameters
    ----------
    vec : ndarray
        input vector

    Returns
    -------
    mag : float
        magnitude
    """

    temp = vec[0]**2 + vec[1]**2 + vec[2]**2

    if abs(temp) >= smalle12:
        mag = math.sqrt(temp)
    else:
        mag = 0.0
    return mag

# ----------------------------------------------------------------------------
#
#                           function polarm
#
#  this function calulates the transformation matrix that accounts for polar
#    motion. both the 1980 and 2000 theories are handled. note that the rotation
#    order is different between 1980 and 2000 .
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    xp          - polar motion coefficient       rad
#    yp          - polar motion coefficient       rad
#    ttt         - julian centuries of tt (00 theory only)
#    opt         - method option                  '01', '02', '80'
#
#  outputs       :
#    pm          - transformation matrix for ecef - pef
#
#  locals        :
#    arcsec2rad  - conversion from arcsec to rad
#    sp          - s prime value
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2004, 207-209, 211, 223-224
#
# [pm] = polarm (xp, yp, ttt, opt)
# ----------------------------------------------------------------------------

### should swap ttt and opt so that ttt isn't necessary for '80' theory -zeg
def polarm(xp: float, yp: float, ttt: float, opt: str):
    """this function calulates the transformation matrix that accounts for
    polar motion. both the 1980 and 2000 theories are handled. note that the
    rotation order is different between 1980 and 2000.

    Parameters
    ----------
    xp : float
        polar motion coefficient: rad
    yp : float
        polar motion coefficient: rad
    ttt : float, optional
        julian centuries of tt (00 theory only): centuries
    opt : str
        theory: '80' for 1980, else 2000

    Returns
    -------
    pm : ndarray
        transformation matrix for ecef - pef/tirs
    """

    cosxp = math.cos(xp)
    sinxp = math.sin(xp)
    cosyp = math.cos(yp)
    sinyp = math.sin(yp)

    pm = np.zeros((3, 3))
    if (opt == "80"):
        pm[0, 0] = cosxp
        pm[0, 1] = 0.0
        pm[0, 2] = -sinxp
        pm[1, 0] = sinxp * sinyp
        pm[1, 1] = cosyp
        pm[1, 2] = cosxp * sinyp
        pm[2, 0] = sinxp * cosyp
        pm[2, 1] = -sinyp
        pm[2, 2] = cosxp * cosyp

        # a1 = rot2mat(xp)
        # a2 = rot1mat(yp)
        # pm = a2*a1
        # Approximate matrix using small angle approximations
        #pm(1, 1) = 1.0
        #pm(2, 1) = 0.0
        #pm(3, 1) = xp
        #pm(1, 2) = 0.0
        #pm(2, 2) = 1.0
        #pm(3, 2) = -yp
        #pm(1, 3) = -xp
        #pm(2, 3) = yp
        #pm(3, 3) = 1.0
    else:
        # approximate sp value in rad
        sp = -47.0e-6 * ttt * arcsec2rad
        print('xp %16.14f, %16.14f sp %16.14g \n' % (xp, yp, sp))
        cossp = math.cos(sp)
        sinsp = math.sin(sp)

        #fprintf(1, ' sp  %14.11f mas \n', sp/arcsec2rad)

        # form the matrix
        pm[0, 0] = cosxp * cossp
        pm[0, 1] = -cosyp * sinsp + sinyp * sinxp * cossp
        pm[0, 2] = -sinyp * sinsp - cosyp * sinxp * cossp
        pm[1, 0] = cosxp * sinsp
        pm[1, 1] = cosyp * cossp + sinyp * sinxp * sinsp
        pm[1, 2] = sinyp * cossp - cosyp * sinxp * sinsp
        pm[2, 0] = sinxp
        pm[2, 1] = -sinyp * cosxp
        pm[2, 2] = cosyp * cosxp

        # a1 = rot1mat(yp)
        # a2 = rot2mat(xp)
        # a3 = rot3mat(-sp)
        # pm = a3*a2*a1
    return pm

# ----------------------------------------------------------------------------
#
#                           function fundarg
#
#  this function calulates the delauany variables and planetary values for
#  several theories.
#
#  author        : david vallado                  719-573-2600   16 jul 2004
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    ttt         - julian centuries of tt
#    opt         - method option                  '06', '02', '96', '80'
#
#  outputs       :
#    l           - delaunay element               rad
#    ll          - delaunay element               rad
#    f           - delaunay element               rad
#    d           - delaunay element               rad
#    omega       - delaunay element               rad
#    planetary longitudes                         rad
#
#  locals        :
#    ttt2, ttt3,  - powers of ttt
#
#  coupling      :
#    none        -
#
#  references    :
#    vallado       2004, 212-214
#
# [ l, l1, f, d, omega, ...
#   lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
# ] = fundarg(ttt, opt)
# ----------------------------------------------------------------------------


def fundarg(ttt: float, opt: str):
    """this function calulates the delauany variables and planetary values for
    several theories.

    Parameters
    ----------
    ttt : float
        julian centuries of tt: centuries
    opt : str
        theory: '06' (2006), '02' (2000b), '96' (1996), '80' (1980)

    Returns
    -------
    l : float
        delauney element
    l1 : float
        delauney element
    f : float
        delauney element
    d : float
        delauney element
    omega : float
        delauney element
    lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep : float
        planetary longitudes
    precrate : float

    """

    # ---- determine coefficients for iau 2000 nutation theory ----
    ttt2 = ttt*ttt
    ttt3 = ttt2*ttt
    ttt4 = ttt2*ttt2

    # ---- iau 2006 theory
    if opt == '06':
        # ------ form the delaunay fundamental arguments in deg
        l = 134.96340251  + (1717915923.2178 *ttt + \
                    31.8792 *ttt2 + 0.051635 *ttt3 - 0.00024470 *ttt4) / 3600.0
        l1 = 357.52910918  + (129596581.0481 *ttt - \
                    0.5532 *ttt2 + 0.000136 *ttt3 - 0.00001149*ttt4)  / 3600.0
        f = 93.27209062  + (1739527262.8478 *ttt - \
                    12.7512 *ttt2 - 0.001037 *ttt3 + 0.00000417*ttt4)  / 3600.0
        d = 297.85019547  + (1602961601.2090 *ttt - \
                    6.3706 *ttt2 + 0.006593 *ttt3 - 0.00003169*ttt4)  / 3600.0
        omega = 125.04455501  + (-6962890.5431 *ttt + \
                    7.4722 *ttt2 + 0.007702 *ttt3 - 0.00005939*ttt4)  / 3600.0

        # ------ form the planetary arguments in deg
        lonmer = 252.250905494  + 149472.6746358  *ttt
        lonven = 181.979800853  +  58517.8156748  *ttt
        lonear = 100.466448494  +  35999.3728521  *ttt
        lonmar = 355.433274605  +  19140.299314   *ttt
        lonjup = 34.351483900  +   3034.90567464 *ttt
        lonsat = 50.0774713998 +   1222.11379404 *ttt
        lonurn = 314.055005137  +    428.466998313*ttt
        lonnep = 304.348665499  +    218.486200208*ttt
        precrate = 1.39697137214*ttt + 0.0003086*ttt2

    # ---- iau 2000b theory
    elif opt == '02':
        # ------ form the delaunay fundamental arguments in deg
        l = 134.96340251  + (1717915923.2178 *ttt) / 3600.0
        l1 = 357.52910918  + (129596581.0481 *ttt) / 3600.0
        f = 93.27209062  + (1739527262.8478 *ttt) / 3600.0
        d = 297.85019547  + (1602961601.2090 *ttt) / 3600.0
        omega = 125.04455501  + ( -6962890.5431 *ttt) / 3600.0

        # ------ form the planetary arguments in deg
        lonmer = 0.0
        lonven = 0.0
        lonear = 0.0
        lonmar = 0.0
        lonjup = 0.0
        lonsat = 0.0
        lonurn = 0.0
        lonnep = 0.0
        precrate = 0.0

    # ---- iau 1996 theory
    elif opt == '96':
        l = 134.96340251  + (1717915923.2178 *ttt + \
                    31.8792 *ttt2 + 0.051635 *ttt3 - 0.00024470 *ttt4) / 3600.0
        l1 = 357.52910918  + (129596581.0481 *ttt - \
                    0.5532 *ttt2 - 0.000136 *ttt3 - 0.00001149*ttt4)  / 3600.0
        f = 93.27209062  + (1739527262.8478 *ttt - \
                    12.7512 *ttt2 + 0.001037 *ttt3 + 0.00000417*ttt4)  / 3600.0
        d = 297.85019547  + (1602961601.2090 *ttt - \
                    6.3706 *ttt2 + 0.006593 *ttt3 - 0.00003169*ttt4)  / 3600.0
        omega = 125.04455501  + ( -6962890.2665 *ttt + \
                    7.4722 *ttt2 + 0.007702 *ttt3 - 0.00005939*ttt4)  / 3600.0
        # ------ form the planetary arguments in deg
        lonmer = 0.0
        lonven = 181.979800853  +  58517.8156748  *ttt   # deg
        lonear = 100.466448494  +  35999.3728521  *ttt
        lonmar = 355.433274605  +  19140.299314   *ttt
        lonjup = 34.351483900  +   3034.90567464 *ttt
        lonsat = 50.0774713998 +   1222.11379404 *ttt
        lonurn = 0.0
        lonnep = 0.0
        precrate = 1.39697137214*ttt + 0.0003086*ttt2

        # ---- iau 1980 theory
    elif opt == '80':
        l = ((((0.064) * ttt + 31.310) * ttt
                + 1717915922.6330) * ttt) / 3600.0 + 134.96298139
        l1 = ((((-0.012) * ttt - 0.577) * ttt
                + 129596581.2240) * ttt) / 3600.0 + 357.52772333
        f = ((((0.011) * ttt - 13.257) * ttt
                + 1739527263.1370) * ttt) / 3600.0 + 93.27191028
        d = ((((0.019) * ttt - 6.891) * ttt
                + 1602961601.3280) * ttt) / 3600.0 + 297.85036306
        omega = ((((0.008) * ttt + 7.455) * ttt
                    - 6962890.5390) * ttt) / 3600.0 + 125.04452222

        # ------ form the planetary arguments in deg
        lonmer = 252.3 + 149472.0 *ttt
        lonven = 179.9 +  58517.8 *ttt   # deg
        lonear = 98.4 +  35999.4 *ttt
        lonmar = 353.3 +  19140.3 *ttt
        lonjup = 32.3 +   3034.9 *ttt
        lonsat = 48.0 +   1222.1 *ttt
        lonurn = 0.0
        lonnep = 0.0
        precrate = 0.0

    # ---- convert units to rad
    l = np.fmod(l, 360.0)     * deg2rad # rad
    l1 = np.fmod(l1, 360.0)    * deg2rad
    f = np.fmod(f, 360.0)     * deg2rad
    d = np.fmod(d, 360.0)     * deg2rad
    omega = np.fmod(omega, 360.0) * deg2rad

    lonmer = np.fmod(lonmer, 360.0) * deg2rad  # rad
    lonven = np.fmod(lonven, 360.0) * deg2rad
    lonear = np.fmod(lonear, 360.0) * deg2rad
    lonmar = np.fmod(lonmar, 360.0) * deg2rad
    lonjup = np.fmod(lonjup, 360.0) * deg2rad
    lonsat = np.fmod(lonsat, 360.0) * deg2rad
    lonurn = np.fmod(lonurn, 360.0) * deg2rad
    lonnep = np.fmod(lonnep, 360.0) * deg2rad
    precrate = np.fmod(precrate, 360.0) * deg2rad

    if sh.iauhelp:
        print('fa %11.7f  %11.7f  %11.7f  %11.7f  %11.7f deg \n'
                % (l*180/math.pi, l1*180/math.pi, f*180/math.pi,
                    d*180/math.pi, omega*180/math.pi))
        print('fa %11.7f  %11.7f  %11.7f  %11.7f deg \n'
                % (lonmer*180/math.pi, lonven*180/math.pi,
                    lonear*180/math.pi, lonmar*180/math.pi))
        print('fa %11.7f  %11.7f  %11.7f  %11.7f deg \n'
                % (lonjup*180/math.pi, lonsat*180/math.pi,
                    lonurn*180/math.pi, lonnep*180/math.pi))
        print('fa %11.7f  \n' % (precrate*180/math.pi))


    # test if they are equivalent
    # most around 1e-10, but some at 1e-6

#         oo3600 = 1.0 / 3600.0
#         deg2rad = math.pi / 180.0
#         ttt = 0.34698738576
#         twopi = 2.0 * math.pi
#         lonmer = mod((908103.259872 + 538101628.688982 * ttt) * oo3600, 360)*deg2rad
#         lonven = mod((655127.283060 + 210664136.433548 * ttt) * oo3600, 360)*deg2rad
#         lonear = mod((361679.244588 + 129597742.283429 * ttt) * oo3600, 360)*deg2rad
#         lonmar = mod((1279558.798488 + 68905077.493988 * ttt) * oo3600, 360)*deg2rad
#         lonjup = mod((123665.467464 + 10925660.377991 * ttt) * oo3600, 360)*deg2rad
#         lonsat = mod((180278.799480 + 4399609.855732 * ttt) * oo3600, 360)*deg2rad
#         lonurn = mod((1130598.018396 + 1542481.193933 * ttt) * oo3600, 360)*deg2rad
#         lonnep = mod((1095655.195728 + 786550.320744 * ttt) * oo3600, 360)*deg2rad
#         precrate = ((1.112022 * ttt + 5028.8200) * ttt) * oo3600*deg2rad
#
#         lonmer1 = mod (4.402608842 + 2608.7903141574 * ttt , twopi)
#         lonven1 = mod (3.176146697 + 1021.3285546211 * ttt , twopi)
#         lonear1 = mod(1.753470314 + 628.3075849991 * ttt , twopi)
#         lonmar1 = mod(6.203480913 + 334.0612426700 * ttt , twopi)
#         lonjup1 = mod(0.599546497 + 52.9690962641 * ttt , twopi)
#         lonsat1 = mod(0.874016757 + 21.3299104960 * ttt , twopi)
#         lonurn1 = mod(5.481293872 + 7.4781598567 * ttt , twopi)
#         lonnep1 = mod(5.311886287 + 3.8133035638 * ttt , twopi)
#         precrate1 = (0.024381750 + 0.00000538691 * ttt) *ttt
#
#         lonmer-lonmer1
#         lonven-lonven1
#         lonear-lonear1
#         lonmar-lonmar1
#         lonjup-lonjup1
#         lonsat-lonsat1
#         lonurn-lonurn1
#         lonnep-lonnep1
#         precrate-precrate1
#
#


    return l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, \
        lonsat, lonurn, lonnep, precrate

# -----------------------------------------------------------------------------
#
#                           procedure finitediff
#
# this procedure perturbs the components of the state vector for processing
# with the finite differencing for the a matrix.
#
#  author        : david vallado                  719-573-2600   15 jan 2008
#
#  inputs          description                    range / units
#    whichconst  - parameter for sgp4 constants   wgs72, wgs721, wgs84
#    pertelem    - which element to perturb
#    percentchg  - amount to modify the vectors   0.001
#                  by in finite differencing
#    deltaamtchg - tolerance for small value in
#                  finite differencing            0.0000001
#    statetype   - type of elements (equinoctial, etc)  'e', 't'
#    xnom        - state vector                   varied
#    scalef      - scale factor for state         all set to 1.0 now
#
#  outputs       :
#    deltaamt    - amount each elemnt is perturbed
#    satrec      - satellite record
#
#  locals        :
#    jj          - index
#
#  coupling      :
#    getgravc    - get the constants for a gravity field for sgp4
#    state2satrec- conversion between state and satellite structure
#    sgp4init    - intiialize sgp4 constants
#
#  references    :
#    vallado       2007, 753-765
# --------------------------------------------------------------------------- */

def finitediff(pertelem: int, percentchg: float, deltaamtchg: float,
               xnom: np.ndarray):
    """this procedure perturbs the components of the state vector for processing
    with the finite differencing for the a matrix.

    Parameters
    ----------
    pertelem : int
        which element to perturb
    percentchg : float
        amount to modify vectors by in finite differencing
    deltaamtchg : float
        tolerance for small matching in finite differencing
    xnom : ndarray
        matrix

    Returns
    -------
    deltaamt: float
        amount each element is perturbed
    xnomp: ndarray
        altered matrix
    """


    # chk if perturbing amt is too small. if so, up the percentchg and try again
# this will execute 5 times, but leaves percentchg the same after each run
    jj = 1
    deltaamt = 0.0
    xnomp = np.copy(xnom)
    while ((np.abs(deltaamt) < deltaamtchg) and (jj < 5)):

        if (jj > 1):
            print('too large\n')
        deltaamt = xnom[pertelem] * percentchg
        xnomp[pertelem, 0] = xnom[pertelem, 0] + deltaamt
        #          state2satrec(xnom, scalef, statetype, statesize, eTo, satrec)
        if (np.abs(deltaamt) < deltaamtchg):
            percentchg = percentchg * 1.4
            print(' %i percentchg chgd %11.8f \n' % (jj, percentchg))
        jj = jj + 1


    # printf(" \n")
    return deltaamt, xnomp

# A quick function to shorten all of the sine and cosine calls made. -zeg
def getsincos(*args):
    """a function that takes any number of angles in radians and returns
    the sine and cosine values.

    Parameters
    ----------
    arg1, arg2, arg3...
        any number of angles in radians

    Returns
    -------
    sin_arg1, cos_arg1, sin_arg2, cos_arg2, sin_arg3, cos_arg3...
        the sine and cosine values of the arguments
    """
    results = []
    for arg in args:
        results.append(math.sin(arg))
        results.append(math.cos(arg))
    return results

def legPoly(x=None, i=None):
    # LEGPOLY calculates the legendre polynomial of i-th order in x
# x can be a scalar, a vector or a matrix

    # pn : nth legendre polynomial
# pn_1 : n-1 legendre polynomial
# pn_p1 : n+1 legendre polynomial

    # Konstantinos G.
# Modified version of the C code in "Numerical recipes in C"

    pn_p1 = 1
    pn = np.ones((x.shape[1-1], x.shape[2-1]))
    if (i > 0):
        pn_1 = np.multiply(x, pn)
        if (i == 1):
            pn = pn_1
        else:
            for i in range(1, i+1):
                pn_p1 = (np.multiply(np.multiply(x, (2 * i - 1)), pn_1)
                         - (i - 1) * pn) / i
                pn = pn_1
                pn_1 = pn_p1
            pn = pn_p1

    return pn


# ELLIPTIC12 evaluates the value of the Incomplete Elliptic Integrals
# of the First, Second Kind and Jacobi's Zeta Function.

#   [F, E, Z] = ELLIPTIC12(U, M, TOL) where U is a phase in radians, 0<M<1 is
#   the module and TOL is the tolerance (optional). Default value for
#   the tolerance is eps = 2.220e-16.

#   ELLIPTIC12 uses the method of the Arithmetic-Geometric Mean
#   and Descending Landen Transformation described in [1] Ch. 17.6,
#   to determine the value of the Incomplete Elliptic Integrals
#   of the First, Second Kind and Jacobi's Zeta Function [1], [2].

#       F(phi, m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);
#       E(phi, m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
#       Z(phi, m) = E(u, m) - E(m)/K(m)*F(phi, m).

#   Tables generating code ([1], pp. 613-621):
#       [phi, alpha] = meshgrid(0:5:90, 0:2:90);                  # modulus and phase in degrees
#       [F, E, Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  # values of integrals

#   See also ELLIPKE, ELLIPJ, ELLIPTIC12I, ELLIPTIC3, THETA, AGM.

#   References:
#   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
#       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
#   [2] D. F. Lawden, "Elliptic Functions and Applications"
#       Springer-Verlag, vol. 80, 1989

# Copyright Elliptic Project 2011
# For support, please reply to
#     moiseev.igor[at]gmail.com
#     Moiseev Igor,

# The code is optimized for ordered inputs produced by the functions
# meshgrid, ndgrid. To obtain maximum performace (up to 30#) for singleton,
# 1-dimensional and random arrays remark call of the function unique(.)
# and edit further code.

def elliptic12(u: np.ndarray, m: np.ndarray, tol: float = None):
    if not isinstance(u, np.ndarray):
        u = np.array([u])

    if not isinstance(m, np.ndarray):
        m = np.array([m])

    if tol == None:
        tol = sys.float_info.epsilon

    if not np.isreal(u).all() or not np.isreal(m).all():
        raise Exception('Input arguments must be real. Use ELLIPTIC12i for complex arguments.')

    # if len(m) == 1:
    #     m = m(np.ones(u.shape))

    # if len(u) == 1:
    #     u = u(np.ones(m.shape))

    # if not m.shape == u.shape:
    #     raise Exception('U and M must be the same size.')

    F = np.zeros(u.shape)
    E = F.copy()
    Z = E.copy()

    m = m.T
    u = u.T

    if np.any(m < 0) or np.any(m > 1):
        raise Exception('M must be in the range 0 <= M <= 1.')

    # cdav change for small eccentricities
    # for i in range(len(m)):
    #     if abs(m[i]) < 1e-07:
    #         m[i] = 1e-07

    I = np.nonzero((m != 1) & (m != 0))[0]
    if len(I):
        # mu, J, K = unique(m(I))
        mu, K = np.unique(m[I], False, True)
        # K = uint32(K)
        mumax = len(mu)
        signU = np.sign(u[I])
        # pre-allocate space and augment if needed
        chunk = 7
        a = np.zeros((chunk, mumax))
        c = a.copy()
        b = a.copy()
        a[0, :] = np.ones(mumax)
        c[0, :] = np.sqrt(mu)
        b[0, :] = np.sqrt(1 - mu)
        n = np.zeros(mumax)
        i = 0
        while (abs(c[i, :]) > tol).any(): # Arithmetic-Geometric Mean of A, B and C

            i = i + 1
            if i > a.shape[0]:
                np.append(a, np.zeros((1, mumax)), axis=0)
                np.append(b, np.zeros((1, mumax)), axis=0)
                np.append(c, np.zeros((1, mumax)), axis=0)

            a[i, :] = 0.5 * (a[i - 1, :] + b[i - 1, :])
            b[i, :] = np.sqrt(a[i - 1, :] * b[i - 1, :])
            c[i, :] = 0.5 * (a[i - 1, :] - b[i - 1, :])
            in_ = np.nonzero((abs(c[i, :]) <= tol)
                             & (abs(c[i - 1, :]) > tol))[0]
            if len(in_):
                n[in_] = i

        mmax = len(I)
        mn = int(np.max(n))
        phin = np.zeros(mmax)
        C = np.zeros(mmax)
        Cp = C.copy()
        e = C.copy()
        phin = signU * u[I]
        i = 0
        c2 = c ** 2
        while i < mn: # Descending Landen Transformation

            in_ = np.nonzero(n[K] > i)[0]
            if len(in_):
                phin[in_] = (np.arctan(b[i, K[in_]] / a[i, K[in_]]
                                         * np.tan(phin[in_]))
                               + np.pi * np.ceil(phin[in_] / np.pi - 0.5)
                               + phin[in_])
                e[in_] = 2.0 ** i
                C[in_] = C[in_] + e[in_[0]] * c2[i, K[in_]]
                Cp[in_] = Cp[in_] + c[i + 1, K[in_]] * np.sin(phin[in_])
            i = i + 1

        Ff = phin / (a[mn, K] * e * 2)
        F[I] = Ff * signU
        Z[I] = Cp * signU
        E[I] = (Cp + ((1 - 1 / 2 * C) *  Ff)) * signU

    # Special cases: m == {0, 1}
    m0 = np.nonzero(m == 0)[0]
    if len(m0):
        F[m0] = u[m0]
        E[m0] = u[m0]
        Z[m0] = 0

    m1 = np.nonzero(m == 1)[0]
    um1 = abs(u[m1])
    if len(m1):
        N = np.floor((um1 + np.pi / 2) / np.pi)
        M = np.nonzero(um1 < np.pi / 2)[0]
        F[m1[M]] = np.log(np.tan(np.pi / 4 + u[m1[M]] / 2))
        F[m1[um1 >= np.pi / 2]] = math.inf *  np.sign(u[m1[um1 >= np.pi / 2]])
        E[m1] = (((-1) ** N * np.sin(um1)) + 2 * N) * np.sign(u[m1])
        Z[m1] = -1 ** N * np.sin(u[m1])

    return F, E, Z

# INVERSELLIPTIC2 evaluates the value of the INVERSE Incomplete Elliptic Integrals
# of the Second Kind.

# INVE = INVERSELLIPTIC2(E,M,TOL) where E is a value of the integral to
# be inverted, 0<M<1 is the module and TOL is the tolerance (optional).
# Default value for the tolerance is eps = 2.220e-16.

# INVERSELLIPTIC2 uses the method described by Boyd J. P.
# to determine the value of the inverse Incomplete Elliptic Integrals
# of the Second Kind using the Empirical initialization to
# the Newtons iteration method [1].

# NOTICE. Please pay attention to the definition of the elliptic functions
# which follows the Abramovitz et al [2], for more theory on elliptic
# functions please consult the Lawden book [3].

# Elliptic integral of the second kind:

# E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);

# Empirical initialization [1]:

# T0(z,m) = pi/2 + sqrt(r)/(theta ? pi/2)

# where
# z \in [?E(pi/2,m), E(pi/2,m)]x[0, 1], value of the entire parameter space
# r = sqrt((1-m)^2 + zeta^2)
# zeta = 1 - z/E(pi/2,m)
# theta = atan((1 - m)/zeta)


# Example:
# # modulus and phase in degrees
# [phi,alpha] = meshgrid(0:5:90, 0:2:90);
# # values of integrals
# [F,E] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);
# # values of inverse
# invE = inverselliptic2(E, sin(pi/180*alpha).^2);
# # the difference between phase phi and invE should close to zero
# phi - invE * 180/pi

# See also ELLIPKE, ELLIPTIC12.

# References:
# [1] J. P. Boyd, "Numerical, Perturbative and Chebyshev Inversion
# of the Incomplete Elliptic Integral of the Second Kind", Applied Mathematics and Computation (January 2012)
# [2] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
# Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
# [3] D. F. Lawden, "Elliptic Functions and Applications"
# Springer-Verlag, vol. 80, 1989

# Copyright (C) 2011 by Elliptic Project. All rights reserved.

# GNU GENERAL PUBLIC LICENSE Version 2, June 1991
# http://www.gnu.org/licenses/gpl.html
# Everyone is permitted to copy and distribute verbatim copies of this
# script under terms and conditions of GNU GENERAL PUBLIC LICENSE.

# For support, please reply to
# moiseev.igor[at]gmail.com
# Moiseev Igor

# ELLIPTIC PROJECT: http://elliptic.googlecode.com
# Group:

def inverselliptic2(E: np.ndarray, m: np.ndarray, tol: float = None):
    if not isinstance(E, np.ndarray):
        E = np.array([E])

    if not isinstance(m, np.ndarray):
        m = np.array([m])

    if tol == None:
        tol = sys.float_info.epsilon

    if not np.isreal(E).all() or not np.isreal(m).all() :
        raise Exception('Input arguments must be real.')

    # if len(m) == 1:
    #     m = m(np.ones(E.shape))

    # if len(E) == 1:
    #     E = E(np.ones(m.shape))

    # if not m.shape == E.shape:
    #     raise Exception('E and M must be the same size.')

    invE = np.zeros(E.shape)
    # make a row vector
    m = m
    E = E

    if np.any(m < 0) or np.any(m > 1):
        raise Exception('M must be in the range 0 <= M <= 1.')

    # cdav change for small eccentricities
    # for i in range(len(m)):
    #     if abs(m[i]) < 1e-07:
    #         m[i] = 1e-07


    # inputs
    z = E
    mu = 1 - m
    # complete integral initialization
    # E1 = special.ellipe(m, tol)
    u = np.full(m.shape, np.pi / 2)
    __, E1, _ = elliptic12(u, m)
    zeta = 1 - z / E1
    r = np.sqrt(zeta * zeta + mu * mu)
    theta = np.arctan(mu / (z + sys.float_info.epsilon))
    # Empirical initialization [1]
    invE = np.pi / 2 + np.sqrt(r) * (theta - (np.pi / 2))
    for _ in range(4):
        __, E, _ = elliptic12(invE, m, tol)
        invE = invE -(E - z) / np.sqrt(1 - m * np.sin(invE) ** 2)

    return invE

#ARCLENGTH_ELLIPSE Calculates the arclength of ellipse.
#
#   ARCLENGTH_ELLIPSE(A, B, THETA0, THETA1) Calculates the arclength of ellipse
#   using the precise formulas based on the representation of
#   the arclength by the Elliptic integral of the second kind.
#
#   Ellipse parameters:
#       T - measured in radians from 0 in the positive direction,
#           Period: 2*Pi
#       A - major axis
#       B - minor axis
#
#   Parametric equations:
#       x(t) = a.cos(t)
#       y(t) = b.sin(t)
#
#   Cartesian equation:
#   x^2/a^2 + y^2/b^2 = 1
#
#   Eccentricity:
#       e = Sqrt(1 - (a/b)^2)
#
#   Focal parameter:
#       b^2/Sqrt(a^2 - b^2)
#
#   Foci:
#       (-Sqrt(a^2 - b^2), 0)   OR   (Sqrt(a^2 - b^2), 0)
#
#   Arclength:
#       b*EllipticE( t, 1 - (a/b)^2 )
#
#   Mathematica Test 1:
#       In:= b = 10; a = 5;
#            SetPrecision[b*EllipticE[2Pi, 1.0- a^2/b^2],20]
#      Out:= 48.442241102738385905
#
#   Mathematica Test 2:
#       In:= b = 10; a = 5;
#            SetPrecision[b*(EllipticE[Pi/2-Pi/10, 1.0- a^2/b^2]-EllipticE[Pi/10, 1.0- a^2/b^2]),20]
#      Out:= 7.3635807913930495516
#
#   MATLAB Test 1:
#       # full ellipse
#       arclength = arclength_ellipse(5,10)
#       arclength =
#           48.442241102738436
#
#   MATLAB Test 2:
#       # arclength ellipse
#       arclength = arclength_ellipse(5,10,pi/10,pi/2)
#       arclength =
#           7.363580791393055
#
#   References:
#   @see http://mathworld.wolfram.com/Ellipse.html
#   @see http://www.wolframalpha.com/input/?i=ellipse+arc+length&lk=1&a=ClashPrefs_*PlaneCurve.Ellipse.PlaneCurveProperty.ArcLength-
#

# Special thanks to for bug correction
#    drbitboy (Brian Carcich) https://github.com/drbitboy
# 2015-07-14 (New Horizons flyby of Pluto)
#
# 1) Old code returned values that were in error
# 1.1)  arclength_ellipse(1., .5, pi*.001, pi*.002) returned 0
# 1.2)  arclength_ellipse(1., .5, pi*.002, pi*.001) returned -.0003*pi instead of pi correct .0005*pi
# 1.3)  arclength_ellipse(1., .5, theta0, theta1) did not return the negative of the same call with the thetas reversed
# 2) Angles theta0 and theta1 were always interpreted as measured from the semi-minor axis
#
# 3) Corrected code:
# 3.1) Angle theta is measured from the positive a axis
# 3.2) The standard form of the b*E(phi,m) arc length integral has m = 1 - (a/b)^2
# 3.2.1) N.B. That only only works if b>a
# 3.3) If a>b, then an alternate formula is used:  a*E(PI/2 - phi, m') where m' = 1 - (b/a)^2
# 3.4) A few simple cases will show that the new code is correct
#        arclength_ellipse(1, .5, pi*.001, pi*.002) ~  pi*.0005
#        arclength_ellipse(1, .5, pi*.002, pi*.001) = -arclength(1, .5, pi*.001, pi*.002) ~ -pi*.0005
#        arclength_ellipse(1., 2., pi*.001, pi*.002) ~ pi*.002
#        arclength_ellipse(1, .5, pi/2 - pi*.002, pi/2 - pi*.001) ~ -pi*.001
#        arclength_ellipse(1, 2., pi/2 - pi*.002, pi/2 - pi*.001) ~ -pi*.001
#        etc.

# Copyright Elliptic Project 2011
# For support,
#     moiseev.igor[at]gmail.com
#     Moiseev Igor

def arclength_ellipse(a: np.ndarray, b: np.ndarray, theta0: np.ndarray = np.array([]),
                      theta1: np.ndarray = np.array([])):

    if not isinstance(a, np.ndarray):
        a = np.array([a])

    if not isinstance(b, np.ndarray):
        b = np.array([b])

    if theta0.size == 0 or theta1.size == 0:
        if theta0.size == 0 and theta1.size == 0:
            theta0 = np.full(a.shape, 0)
            theta1 = np.full(a.shape, 2*math.pi)
        else:
            print("Error: requires both theta0 and theta1 set or neither!")
            return

    if not isinstance(theta0, np.ndarray):
        theta0 = np.array([theta0])

    if not isinstance(theta1, np.ndarray):
        theta1 = np.array([theta1])

    arclength = a * (theta1 - theta0)


    theta0a = np.zeros(len(a))
    theta1a = np.zeros(len(a))
    ab = np.zeros(len(a))
    for i in range(len(a)):
        if a[i] <= b[i]:
            ab[i] = 1 - (a[i] / b[i]) ** 2
            theta0a[i] = theta0[i]
            theta1a[i] = theta1[i]
        elif a[i] > b[i]:
            ab[i] =  1 - (b[i] / a[i]) ** 2
            theta0a[i] = math.pi / 2 - theta0[i]
            theta1a[i] = math.pi / 2 - theta1[i]

    F1, E1, Z1 = elliptic12(theta1a, ab)
    F0, E0, Z0 = elliptic12(theta0a, ab)

    arclength = np.zeros(len(a))
    for i in range(len(a)):
        if a[i] <= b[i]:
            arclength[i] = b[i] * (E1[i] - E0[i])
        elif a[i] > b[i]:
            arclength[i] = a[i] * (E0[i] - E1[i])

    return arclength


# ------------------------------------------------------------------------------
#
#                                  findlos
#
#  finds the line of site vector given an array or arrays of right topocentric
#  ascension and declination angles in eci
#
#  author        : michael courville                            jul 11, 2024
#
#  revisions
#                -
#
#  inputs               description                             range / units
#    radec_angles        - array(s) of rtasc and decl angles    rad
#                        - np.array([rtasc,decl])
#
#  outputs       :
#    los                - line of site unit vector(s)           unitless
#
#
# [los1x,los1y,los1z], [los2x, los2y, los2z],... =
#           findlos(np.array([rtasc1, decl1]), np.array([rtasc2, decl2]), ...)
# -----------------------------------------------------------------------------
def findlos(*radec_angles: np.ndarray):
    """finds the line of site vector given an array or arrays
    of right topocentric ascension and declination angles in eci

    Parameters
    ----------
    radec_angles: ndarray([float, float)]
        array(s) of rtasc and decl angles: rad
        ex: np.array([rtasc,decl])
    Returns
    -------
    los: ndarray([float, float, float)]
        line of site unit vector(s)
    """

    los = np.zeros((len(radec_angles),3))

    i = 0
    for arg in radec_angles:
        rtasc = arg[1]
        decl = arg[1]

        los[i,0] = math.cos(decl) * math.cos(rtasc)
        los[i,1] = math.cos(decl) * math.sin(rtasc)
        los[i,2] = math.sin(decl)
        i = i + 1
    return los



##############################################################################################################
##############################################################################################################
##############################################################################################################

if __name__ == '__main__':


    c2new, c3new = findc2c3(math.pi/7.0)
    print("findc2c3 returned: ", c2new, c3new)


    m, nu = newtone(0.4, 334.566986 * deg2rad)
    print("newtone returned: ", m, nu)

    e0, nu = newtonm(0.1, math.pi/2.0)
    print("newtonm returned: ", e0, nu)

    e0, ma = newtonnu(0.1, math.pi/2.0)
    print("newtonnu returned: ", e0, ma)

    jd = 60206

    r1r, r1i, r2r, r2i, r3r, r3i = cubic(1, 2, 3, 4, 'I')
    print("cubic returned: ", r1r, r1i, r2r, r2i, r3r, r3i)

    r1r, r1i, r2r, r2i = quadric(1, 2, 3, 'I')
    print("quadric returned: ", r1r, r1i, r2r, r2i)


    ttt = (jd - 2451545.0)/ 36525.0  #0.34698738576
    print("ttt is ", ttt)
    opt = "80"

    l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, \
        lonsat, lonurn, lonnep, precrate = fundarg(ttt, opt)
    print("fundarg returned ", l, l1, f, d, omega, lonmer, lonven,
            lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate)

    xp = math.pi
    yp = math.pi

    pm = polarm (xp, yp, ttt, opt)
    print("polarm returned ", pm)


    vec = np.array([1, 2, 3])
    print("ORIG: ", vec)

    print("ANGL---")
    vecx = angl(vec, vec)
    print("ZERO: ", vecx)

    print("ROT1---")
    vecx = rot1(vec, 0.0)
    print("ZERO: ", vecx)
    vecx = rot1(vec, math.pi/2.0)
    print("PIO2: ", vecx)
    vecx = rot1(vec, math.pi)
    print("PI  : ", vecx)
    vecx = rot1(vec, 2.0*math.pi)
    print("2PI : ", vecx)

    print("ROT2---")
    vecx = rot2(vec, 0.0)
    print("ZERO: ", vecx)
    vecx = rot2(vec, math.pi/2.0)
    print("PIO2: ", vecx)
    vecx = rot2(vec, math.pi)
    print("PI  : ", vecx)
    vecx = rot2(vec, 2.0*math.pi)
    print("2PI : ", vecx)

    print("ROT3---")
    vecx = rot3(vec, 0.0)
    print("ZERO: ", vecx)
    vecx = rot3(vec, math.pi/2.0)
    print("PIO2: ", vecx)
    vecx = rot3(vec, math.pi)
    print("PI  : ", vecx)
    vecx = rot3(vec, 2.0*math.pi)
    print("2PI : ", vecx)


    print("MAG---")
    vec = np.array([math.sqrt(1.1e-16), math.sqrt(1.1e-16), math.sqrt(1.1e-16)])
    vecx = mag(vec)
    print("MAG : ", vecx)
    vec = np.array([math.sqrt(1.0e-16/3.0), math.sqrt(1.0e-16/3.0),
                    math.sqrt(1.0e-16/3.0)])
    vecx = mag(vec)
    print("MAG : ", vecx)

    #arclength_ellipse tests
    answers = np.array([.0005, -.0005, .002, .001, .001, 48.4422411 / math.pi,
                        7.36358/math.pi])
    a = np.array([1, 1, 1, 1, 1, 5, 5])
    b = np.array([.5, .5, 2, .5, 2, 10, 10])
    theta0 = np.array([.001, .002, .001, .5 - .002, .5 - .002, 0, 0.1])
    theta1 = np.array([.002, .001, .002, .5 - .001, .5 - .001, 2, 0.5])
    theta0 = theta0 * math.pi
    theta1 = theta1 * math.pi

    arclength = arclength_ellipse(a, b, theta0, theta1)

    print(f'{answers}')
    print(f'{arclength / math.pi}')
    print(f'{(arclength / math.pi) / answers}')






















