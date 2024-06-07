import math, os
import numpy as np
from pprint import pprint as pp
from space_constants import *
import spacetime_utils as stu
import sethelp as sh
import spacemath_utils as smu
import space_conversions as sc

# this could be done much easier with the pandas package, but i wanted to maintain
# less packages (fwf =fixed width file)
# lins is the fwf string, maxsp is the max space size between chars (round up)
def clean_fwf(lins, maxsp, isint=True, isfloat=False):
  ints = lins
  for ii in range(maxsp+1, 1, -1):
    ints = ints.replace(' '*ii, ' ')
  ints = ints.strip().replace(' ', ',')
  ints = ints.split(',')
  if isfloat:
    ints = [float(xx) for xx in ints]
  if isint:
    ints = [int(xx) for xx in ints]
  return ints



# ------------------------------------------------------------------------------
#
#                           function findtof
#
#  this function finds the time of flight given the initial position vectors,
#    semi-parameter, and the sine and cosine values for the change in true
#    anomaly.  the result uses p-iteration theory to analytically find the result.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#    vallado     - fix tolerances                                 5 sep 2002
#
#  inputs          description                    range / units
#    ro          - interceptor position vector    km
#    r           - target position vector         km
#    p           - semiparameter                  km
#
#  outputs       :
#    tof         - time for transfer              sec
#
#  locals        :
#    sindnu      - sine of change in nu           rad
#    cosdnu      - cosine of change in nu         rad
#    deltae      -
#    deltah      -
#    k           -
#    l           -
#    m           -
#    a           -
#    f           -
#    g           -
#    fdot        -
#    sindeltae   - sine value
#    cosdeltae   - cosine value
#    rcrossr     - cross product of two positions
#
#  coupling      :
#    cross       - cross product of two vectors
#    sinh        - hyperbolic sine
#    arccosh     - arc hyperbolic cosine
#
#  references    :
#    vallado       2007, 134-135, alg 11
#
# [tof] = findtof (ro, r, p)
# ------------------------------------------------------------------------------

def findtof (ro, r, p):
        small = 0.00000001

        # -------------------------  implementation   -----------------
        magr = smu.mag(r)
        magro = smu.mag(ro)
        cosdnu = np.dot(ro, r)/(magro*magr)
        rcrossr = np.cross(ro, r)
        sindnu = smu.mag(rcrossr)/(magro*magr)

        k = magro * magr*(1.0 -cosdnu)
        l = magro + magr
        m = magro * magr*(1.0 +cosdnu)
        a = (m*k*p) / ((2.0 *m-l*l)*p*p + 2.0 *k*l*p - k*k)

        # ------  use f and g series to find velocity vectors  --------
        f = 1.0  - (magr/p)*(1.0 -cosdnu)
        g = magro*magr*sindnu/math.sqrt(mu * p)
        alpha = 1.0 /a

        if (alpha > small):
            # ------------------------ elliptical ---------------------
            dnu = math.atan2(sindnu, cosdnu)
            fdot = math.sqrt(mu / p) * math.tan(dnu*0.5)* \
                    (((1.0 - cosdnu) / p)-(1.0 / magro)-(1.0 / magr))
            cosdeltae = 1.0 - (magro/a)*(1.0 -f)
            sindeltae = (-magro*magr*fdot)/math.sqrt(mu * a)
            deltae = math.atan2(sindeltae, cosdeltae)
            tof = g + math.sqrt(a*a*a/mu)*(deltae-sindeltae)
        else:
            # ------------------------ hyperbolic ---------------------
            if (alpha < -small):
                deltah = math.acosh(1.0 - (magro/a) * (1.0-f))
                tof = g + math.sqrt(-a*a*a/mu)*(math.sinh(deltah)-deltah)
            else:
                # -------------------- parabolic ----------------------
                dnu = math.atan2(sindnu, cosdnu)
                c = math.sqrt(magr*magr+magro*magro - \
                           2.0 *magr*magro*math.cos(dnu))
                s = (magro+magr+c) * 0.5
                tof = (2.0 /3.0) * math.sqrt(s*s*s*0.5/mu) * \
                           (1.0  -  ((s-c)/s)^1.5)
        return tof



# ------------------------------------------------------------------------------
#
#                           function makeorbitrv
#
#  this function propagtes an orbit around for one rev using either kepler
#  or pkepler (secular j2).
#
#  author        : david vallado                  719-573-2600    14 mar 2006
#
#  inputs          description                    range / units
#    rint        - initial position vector of int km
#    vint        - initial velocity vector of int km/s
#    rtgt        - initial position vector of tgt km
#    vtgt        - initial velocity vector of tgt km/s
#    dm          - direction of motion for gauss  'l', 's'
#    kind        - type of propagator             'k', 'p'
#    dtsec       - time of flight to the int      s
#
#  outputs       :
#    v1t         - initial transfer velocity vec  km/s
#    v2t         - final transfer velocity vec    km/s
#    dv1         - initial change velocity vec    km/s
#    dv2         - final change velocity vec      km/s
#    error       - error flag from gauss          'ok', ...
#
#  locals        :
#
#  coupling      :
#    kepler      - find r and v at future time

#
#  references    :
#    vallado       2004,
#
# [x, y, z] = makeorbitrv (2455545.0, 'j', [5003.400903511, -3817.812007872, 4720.200666830], [5.489294908, 3.005055561, -3.39013016])
# ------------------------------------------------------------------------------

def makeorbitrv(jd=None, kind=None, reci=None, veci=None):
    # -------------------------  implementation   -------------------------
    error = 'ok'
    if (kind == 'k'):
        outfile = open('makeorbrvKep.out', 'wt')
    else:
        if (kind == 'p'):
            outfile = open('makeorbrvJ2.out', 'wt')
        else:
            if (kind == 'j'):
                outfile = open('makeorbrvJ4.out', 'wt')

    # --------------- approximate ast with gst for this simple demo --------------
#    gmst = gstime(jd)
#    st(1, 1) = cos(gmst)
#    st(1, 2) = -sin(gmst)
#    st(1, 3) = 0.0
#    st(2, 1) = sin(gmst)
#    st(2, 2) = cos(gmst)
#    st(2, 3) = 0.0
#    st(3, 1) = 0.0
#    st(3, 2) = 0.0
#    st(3, 3) = 1.0
#    recef = st'*reci

    #    fprintf(1, 'st gmst #11.7f ast #11.7f ome  #11.7f \n', gmst*180/pi, ast*180/pi, omegaearth*180/pi)

    x = np.zeros((3, 3))
    y = np.zeros((3, 3))
    z = np.zeros((3, 3))
    x[0, 0] = reci[0]
    y[0, 0] = reci[1]
    z[0, 0] = reci[2]
    ndot = 0.0
    nddot = 0.0
    #    period = 2.0*pi*sqrt(a*a*a/398600.4418)
    jdut1 = jd
    outfile.write(' 1601   xx \n' % ())
    outfile.write('Time (UTCG)                x (km)             y (km)             z (km)         vx (km/sec)     vy (km/sec)     vz (km/sec)\n' % ())
    outfile.write('-------------------------    ---------------    ---------------    ---------------    ------------    ------------    ------------\n' % ())
    for i in range(101):
        #        dtsec = (i-1)*period/120.0
        dtsec = i * 960.0
        jdut1 = jd + dtsec / 86400.0
        # ----------- propogate satellite forward by time ----------------
        if (kind == 'k'):
            reci1, veci1 = kepler(reci, veci, dtsec)
        if (kind == 'p'):
            reci1, veci1 = pkepler(reci, veci, dtsec, ndot, nddot)
        if (kind == 'j'):
            reci1, veci1 = pkeplerj4(reci, veci, dtsec, ndot, nddot)
        #        fprintf(1, 'reci1 #4i x #11.7f  #11.7f  #11.7f  #11.7f  #11.7f  #11.7f \n', i, reci1, veci1)
#        gmst = gstime(jd+dtsec/86400.0)
#        st(1, 1) = cos(gmst)
#        st(1, 2) = -sin(gmst)
#        st(1, 3) = 0.0
#        st(2, 1) = sin(gmst)
#        st(2, 2) = cos(gmst)
#        st(2, 3) = 0.0
#        st(3, 1) = 0.0
#        st(3, 2) = 0.0
#        st(3, 3) = 1.0
#        recef = st'*reci1'
        #         [latgc, latgd, lon, hellp] = ijk2ll (reci1, jdut1)
#         if lon < 0.0
#             lon =lon + 2*pi
#         end
#        fprintf(outfile, '#11.7f  #11.7f   1.0 \n', latgd * rad2deg, lon * rad2deg)
        h, m, s = stu.sec2hms(dtsec)
        outfile.write('1 Jun 2018 %2i:%2i:%7.4f  %16.9f  %16.9f %16.9f  %16.9f %16.9f  %16.9f \n' % \
                    (h, m, s, reci1[0], reci1[1], reci1[2], veci1[0], veci1[1], veci1[2]))
        #        x(i, 1) = reci1(1)
#        y(i, 1) = reci1(2)
#        z(i, 1) = reci1(3)


    return x, y, z


#
# ----------------------------------------------------------------------------
#
#                           function iau06in
#
#  this function initializes the matricies needed for iau 2006 reduction
#    calculations. the routine uses the files listed as inputs, but they are
#    are not input to the routine as they are static files.
#
#  author        : david vallado                  719-573-2600   16 jul 2004
#
#  revisions
#    dav 14 apr 11  update for iau2006 conventions
#
#  inputs          description                    range / units
#    none
#    iau00x.dat  - file for x coefficient
#    iau00y.dat  - file for y coefficient
#    iau00s.dat  - file for s coefficient
#    iau00n.dat  - file for nutation coefficients
#    iau00pl.dat notused - file for planetary nutation coefficients
#    iau00gs.dat - file for gmst coefficients
#
#  outputs       :
#    axs0        - real coefficients for x        rad
#    a0xi        - integer coefficients for x
#    ays0        - real coefficients for y        rad
#    a0yi        - integer coefficients for y
#    ass0        - real coefficients for s        rad
#    a0si        - integer coefficients for s
#    apn         - real coefficients for nutation rad
#    apni        - integer coefficients for nutation
#    ape         - real coefficients for obliquity rad
#    apei        - integer coefficients for obliquity
#    agst        - real coefficients for gst      rad
#    agsti       - integer coefficients for gst
#
#  locals        :
#    convrt      - conversion factor to radians
#    i           - index
#
#  coupling      :
#    none        -
#
#  references    :
#    vallado     2004, pg 205-219, 910-912
#
# [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti] = iau06in
# -----------------------------------------------------------------------------

def iau06in():
    #function [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, ape, apei, agst, agsti] = iau06in


    # ------------------------  implementation   -------------------
# " to rad
    convrtu = (1e-06 * np.pi) / (180.0 * 3600.0)

    convrtm = (0.001 * np.pi) / (180.0 * 3600.0)

    # ------------------------------
#  note that since all these coefficients have only a single
#  decimal place, one could store them as integres, and then simply
#  divide by one additional power of ten. it woul dmake memeory
#  storage much smaller and potentially faster.
# ------------------------------

    # xys values
    fn = os.path.join(os.path.dirname(__file__), "data", "iau06xtab5.2.a.dat")
    filein = open(fn, 'r').readlines()
    axs0 = []
    a0xi = []

    for line in filein:
        if line == "" or line == "\n": continue
        if line.startswith("#"): continue
        f1 = float(line[5:20].strip())
        f2 = float(line[20:35].strip())
        axs0.append([f1*convrtu, f2*convrtu])
        ints = clean_fwf(line[35:], 4)
        #print(ints, len(ints))
        a0xi.append(ints)
    axs0 = np.array(axs0)
    a0xi = np.array(a0xi)
    pp(axs0.shape)
    pp(a0xi.shape)

    fn = os.path.join(os.path.dirname(__file__), "data", "iau06ytab5.2.b.dat")
    filein = open(fn, 'r').readlines()
    ays0 = []
    a0yi = []
    for line in filein:
        if line == "" or line == "\n": continue
        if line.startswith("#"): continue
        f1 = float(line[5:20].strip())
        f2 = float(line[20:35].strip())
        ays0.append([f1*convrtu, f2*convrtu])
        ints = clean_fwf(line[35:], 4)
        #print(ints, len(ints))
        a0yi.append(ints)
    ays0 = np.array(ays0)
    a0yi = np.array(a0yi)
    pp(ays0.shape)
    pp(a0yi.shape)

    fn = os.path.join(os.path.dirname(__file__), "data", "iau06stab5.2.d.dat")
    filein = open(fn, 'r').readlines()
    ass0 = []
    a0si = []
    for line in filein:
        if line == "" or line == "\n": continue
        if line.startswith("#"): continue
        f1 = float(line[5:20].strip())
        f2 = float(line[20:35].strip())
        ass0.append([f1*convrtu, f2*convrtu])
        ints = clean_fwf(line[35:], 4)
        #print(ints, len(ints))
        a0si.append(ints)
    ass0 = np.array(ass0)
    a0si = np.array(a0si)
    pp(ass0.shape)
    pp(a0si.shape)

    # nutation values old approach iau2003
    fn = os.path.join(os.path.dirname(__file__), "data", "iau03n.dat")
    filein = open(fn, 'r').readlines()
    apni = []
    apn = []
    for line in filein:
        if line == "" or line == "\n": continue
        ints = clean_fwf(line[0:16], 2)
        #print(ints, len(ints))
        apni.append(ints)
        flts = clean_fwf(line[29:], 6, isfloat = True, isint = False)
        flts = [fff*convrtm for fff in flts]
        #print(flts, len(flts))
        apn.append(flts)
    apn = np.array(apn)
    apni = np.array(apni)
    pp(apn.shape)
    pp(apni.shape)

    # planetary nutation values
    # original code extricated 5 int values but only used 4. i use all 5.
    fn = os.path.join(os.path.dirname(__file__), "data", "iau03pl.dat")
    filein = open(fn, 'r').readlines()
    appli = []
    appl = []
    for line in filein:
        if line == "" or line == "\n": continue
        ints = clean_fwf(line[4:60], 2)
        #print(ints, len(ints))
        appli.append(ints)
        flts = clean_fwf(line[72:], 6, isfloat = True, isint = False)
        flts = [fff*convrtm for fff in flts]
        #print(flts, len(flts))
        appl.append(flts)
    appl = np.array(appl)
    appli = np.array(appli)
    pp(appl.shape)
    pp(appli.shape)

    # nutation values planetary now included new iau2006
#       load iau00n.dat  # luni-solar
#       apn = iau00n(:, 2:3)
#       apni = iau00n(:, 4:17)
#       for i =1:size(apn)
#           apn[i, 0] = apn[i, 0] * convrtu
#           apn[i, 1] = apn[i, 1] * convrtu
#       end

    #       load iau00e.dat  # planetary
#       ape = iau00n(:, 2:3)
#       apei = iau00n(:, 4:17)
#       for i =1:size(ape)
#           ape[i, 0] = ape[i, 0] * convrtu
#           ape[i, 1] = ape[i, 1] * convrtu
#       end

    # gmst values
# note - these are very similar to the first 34 elements of iau00s.dat,
# but they are not the same.
    fn = os.path.join(os.path.dirname(__file__), "data", "iau06gsttab5.2.e.dat")
    filein = open(fn, 'r').readlines()
    agst = []
    agsti = []
    for line in filein:
        if line == "" or line == "\n": continue
        if line.startswith("#"): continue
        f1 = float(line[5:20].strip())
        f2 = float(line[20:34].strip())
        agst.append([f1*convrtu, f2*convrtu])
        ints = clean_fwf(line[34:], 4)
        #print(ints, len(ints))
        agsti.append(ints)
    agst = np.array(agst)
    agsti = np.array(agsti)
    pp(agst.shape)
    pp(agsti.shape)

    return axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti



# -----------------------------------------------------------------------------
#
#                           function getgravc
#
#  this function gets constants for the propagator. note that mu is identified to
#    facilitiate comparisons with newer models.
#
#  author        : david vallado                  719-573-2600   21 jul 2006
#
#  inputs        :
#    whichconst  - which set of constants to use  721, 72, 84
#
#  outputs       :
#    tumin       - minutes in one time unit
#    mu          - earth gravitational parameter
#    radiusearthkm - radius of the earth in km
#    xke         - reciprocal of tumin
#    j2, j3, j4  - un-normalized zonal harmonic values
#    j3oj2       - j3 divided by j2
#
#  locals        :
#
#  coupling      :
#
#  references    :
#    norad spacetrack report #3
#    vallado, crawford, hujsak, kelso  2006
# [tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(whichconst)
#  --------------------------------------------------------------------------- */

def getgravc(whichconst=None):
    #global tumin mu radiusearthkm xke j2 j3 j4 j3oj2
    if 721 == whichconst:
        # -- wgs-72 low precision str#3 constants --
        mu = 398600.79964
        radiusearthkm = 6378.135
        xke = 0.0743669161
        tumin = 1.0 / xke
        j2 = 0.001082616
        j3 = - 2.53881e-06
        j4 = - 1.65597e-06
        j3oj2 = j3 / j2
    else:
        if 72 == whichconst:
            # ------------ wgs-72 constants ------------
            mu = 398600.8
            radiusearthkm = 6378.135
            xke = 60.0 / np.sqrt(radiusearthkm * radiusearthkm
                                 * radiusearthkm / mu)
            tumin = 1.0 / xke
            j2 = 0.001082616
            j3 = - 2.53881e-06
            j4 = - 1.65597e-06
            j3oj2 = j3 / j2
        else:
            if 84 == whichconst:
                # ------------ wgs-84 constants ------------
                mu = 398600.5
                radiusearthkm = 6378.137
                xke = 60.0 / np.sqrt(radiusearthkm * radiusearthkm
                                     * radiusearthkm / mu)
                tumin = 1.0 / xke
                j2 = 0.00108262998905
                j3 = - 2.53215306e-06
                j4 = - 1.61098761e-06
                j3oj2 = j3 / j2
            else:
                print('unknown gravity option (%d)\n' % (whichconst))


    return tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2


# ------------------------------------------------------------------------------
#
#                           function pathm
#
#  this function determines the end position for a given range and azimuth
#    from a given point.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    llat        - start geocentric latitude      -pi/2 to  pi/2 rad
#    llon        - start longitude (west -)       0.0  to 2pi rad
#    range       - range between points           er
#    az          - azimuth                        0.0  to 2pi rad
#
#  outputs       :
#    tlat        - end geocentric latitude        -pi/2 to  pi/2 rad
#    tlon        - end longitude (west -)         0.0  to 2pi rad
#
#  locals        :
#    sindeltan   - sine of delta n                rad
#    cosdeltan   - cosine of delta n              rad
#    deltan      - angle between the two points   rad
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2001, 774-776, eq 11-6, eq 11-7
#
# [tlat, tlon] = pathm (llat, llon, range, az)
# ------------------------------------------------------------------------------

def pathm(llat=None, llon=None, range_=None, az=None):
    # -------------------------  implementation   -----------------
    small = 1e-08
    az = np.fmod(az, twopi)
    if (llon < 0.0):
        llon = twopi + llon

    if (range_ > twopi):
        range_ = np.fmod(range_, twopi)

    # ----------------- find geocentric latitude  -----------------
    tlat = np.arcsin(np.sin(llat) * np.cos(range_) + np.cos(llat)
                     * np.sin(range_) * np.cos(az))
    # ---- find delta n, the angle between the points -------------
    if (((np.abs(np.cos(tlat)) > small) and (np.abs(np.cos(llat)) > small))):
        sindn = np.sin(az) * np.sin(range_) / np.cos(tlat)
        cosdn = ((np.cos(range_) - np.sin(tlat) * np.sin(llat))
                 / (np.cos(tlat) * np.cos(llat)))
        deltan = math.atan2(sindn, cosdn)
    else:
        # ------ case where launch is within 3nm of a pole --------
        if (np.abs(np.cos(llat)) <= small):
            if (((range_ > np.pi) and (range_ < twopi))):
                deltan = az + np.pi
            else:
                deltan = az
        # ----- case where end point is within 3nm of a pole ------
        if (np.abs(np.cos(tlat)) <= small):
            deltan = 0.0

    tlon = llon + deltan
    if (np.abs(tlon) > twopi):
        tlon = np.fmod(tlon, twopi)

    if (tlon < 0.0):
        tlon = twopi + tlon

    return tlat, tlon


# ------------------------------------------------------------------------------
#
#                           function satfov
#
#  this function finds parameters relating to a satellite's fov.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    incl        - inclination                    rad
#    az          - azimuth                        rad
#    slatgd      - geodetic latitude of sat       rad
#    slon        - longitude of sat               rad
#    salt        - altitudeof satellite           er
#    tfov        - total field of view            rad
#    etactr      - ctr where sensor looks         rad
#
#  outputs       :
#    fovmax      - maximum field of view          rad
#    totalrng    -
#    rhomax      -                                km
#    rhomin      -                                km
#    tgtlat      -
#    tgtlon      -
#
#  locals        :
#    r           -
#    etahopriz   -
#    rhohoriz    -
#    gamma       -
#    rho         -
#    fovmin      -
#    lat         -
#    lon         -
#    maxlat      -
#    minlkat     -
#    i           - index
#
#  coupling      :
#    path        - finds tgt location given initial location, range, and az
#
#  references    :
#    vallado       2001, 776-781, eq 11-8 to eq 11-13, ex 11-1
#
# [totalrng, rhomax, rhomin, tgtlat, tgtlon] = ...
#  satfov (incl, az, slatgd, slon, salt, tfov, etactr, fovmax)
# ------------------------------------------------------------------------------

def satfov(incl=None, az=None, slatgd=None,
           slon=None, salt=None, tfov=None, etactr=None):
    # -------------------------  implementation   -----------------
# ------- find satellite parameters and limiting cases --------
    r = re + salt
    etahoriz = np.arcsin(re / r)
    rhohoriz = r * np.cos(etahoriz)

    print('etahoriz %11.7f rhohoriz %11.7f km \n'
          % (etahoriz * rad2deg, rhohoriz))
    # ---------------- find ground range angle --------------------
    lambda_ = np.arccos(re / r)
    print('lambda %11.7f  %11.7f km  \n' % (lambda_ * rad2deg, lambda_ * re))
    print('maximum locations \n' % ())
    # -------- for maximum, if the sensor looks off axis ----------
    fovmin = etactr + tfov * 0.5
    gamma = np.pi - np.arcsin(r * np.sin(fovmin) / re)

    rhomax = re * np.cos(gamma) + r * np.cos(fovmin)
    print('fovmin %11.7f gamma %11.7f gamma %11.7f rho %11.7f  \n'
          % (fovmin * rad2deg, gamma * rad2deg, (np.pi - gamma) * rad2deg, rhomax))
    # --------------------- slant range --------------------
    lambda_ = np.arcsin(rhomax * np.sin(fovmin) / re)
    rhomin = lambda_ * re
    print('lambda %11.7f rhomin %11.7f \n' % (lambda_ * rad2deg, rhomin))
    #         end

    # -------------- find location of center of fov ---------------
    if (np.abs(etactr) > 1e-05):
        lat, lon = pathm(slatgd, slon, lambda_, az)
    else:
        lat = slatgd
        lon = slon

    print('max NS lat %11.7f lon %11.7f \n'
          % ((lat + lambda_) * rad2deg, lon * rad2deg))
    print('min NS lat %11.7f lon %11.7f \n'
          % ((lat - lambda_) * rad2deg, lon * rad2deg))
    print('max EW lat %11.7f lon %11.7f \n'
          % (lat * rad2deg, (lon + lambda_) * rad2deg))
    print('min EW lat %11.7f lon %11.7f \n'
          % (lat * rad2deg, (lon - lambda_) * rad2deg))
    print('sat ctr of fov lat %11.7f lon %11.7f \n'
          % (lat * rad2deg, lon * rad2deg))
    return rhomin, rhomax



# ------------------------------------------------------------------------------
#
#                           function lambertu
#
#  this function solves the lambert problem for orbit determination and returns
#    the velocity vectors at each of two given position vectors.  the solution
#    uses universal variables for calculation and a bissection technique
#    updating psi.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                         range / units
#    r1          - ijk position vector 1               km
#    v1          - ijk velocity vector 1               needed for 180 deg transfer  km / s
#    r2          - ijk position vector 2               km
#    dm          - direction of motion                 'S', 'L'  short long period
#    df/de          - direction of flight                 'L', 'H'  low high energy
#    dtsec       - time between r1 and r2              s
#    nrev        - multiple revoluions                 0, 1, ...
#    tbi         - time of the bottom interval -       only needed for multi-rev cases
#                  this is a two-dimension array of psi and tof
#  outputs       :
#    v1          - ijk velocity vector                 km / s
#    v2          - ijk velocity vector                 km / s
#    error       - error flag                          'ok', ...
#
#  locals        :
#    vara        - variable of the iteration,
#                  not the semi-axis
#    y           - area between position vectors
#    upper       - upper bound for z
#    lower       - lower bound for z
#    cosdeltanu  - cosine of true anomaly change        rad
#    f           - f expression
#    g           - g expression
#    gdot        - g dot expression
#    x           - old universal variable x
#    xcubed      - x cubed
#    zold        - old value of z
#    znew        - new value of z
#    c2          - c2(z) function
#    c3          - c3(z) function
#    timenew     - new time                             s
#    small       - tolerance for roundoff errors
#    i, j        - index
#
#  coupling      :
#    mag         - magnitude of a vector
#    dot         - dot product of two vectors
#    findc2c3    - find c2 and c3 functions
#
#  references    :
#    vallado       2013, 489-493, alg 58, ex 7-5
#
# [vo, v, errorl] = lambertu (r1, v1, r2, dm, df, nrev, dtsec, tbi)
# ------------------------------------------------------------------------------


def lambertu (r1, v1, r2, dm, de, nrev, dtwait, dtsec, tbi, outfile):

    # -------------------------  implementation   -------------------------
    numiter = 20
    errorl = '      ok'
    v1dv = np.zeros((3))
    v2dv = np.zeros((3))

    dtold = 0.0

    # try canonical units for testing
    #constastro
    #mu = 1.0
    #r1 = r1/re
    #r2 = r2/re
    #dtsec = dtsec / tusec

    # ---- find parameters that are constant for the initial geometry
    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)

    # this value stays constant in all calcs, vara changes with df
    cosdeltanu = np.dot(r1, r2) / (magr1*magr2)
    if abs(cosdeltanu) > 1.0:
        cosdeltanu = np.sign(cosdeltanu) * 1.0
    if (dm == 'L'):  #dm == 'l'
        vara = -math.sqrt(magr1*magr2*(1.0 + cosdeltanu))
    else:
        vara = math.sqrt(magr1*magr2*(1.0 + cosdeltanu))

    # setup variables for speed
    oomu = 1.0 / math.sqrt(mu)

    # --------- set up initial bounds for the bissection ----------
    if (nrev == 0):
        lower = -16.0*math.pi*math.pi  # allow hyperbolic and parabolic solutions
        upper = 4.0*math.pi*math.pi  # could be negative infinity for all cases
    else:
        # set absolute limits for multi-rev cases
        lower = 4.0*nrev**2*math.pi*math.pi
        upper = 4.0*(nrev + 1.0)**2*math.pi*math.pi
        # adjust based on long or short way if dm == 'l'
        #if ((dm == 'l') and (df == 'd')) or ((dm == 's') and (df == 'r'))
        #if ((df == 'r') and (dm == 's')) or ((df == 'd') and (dm == 'l'))
        if (de == 'H'): #  and (dm == 'L')) or ((de == 'L') and (dm == 'L'))
            upper = tbi[nrev, 0]
        else:
            lower = tbi[nrev, 0]

    # ---------------  form initial guesses   ---------------------
    dtdpsi = 0.0
    xold = 0.0
    psinew = 0.0
    if (nrev == 0):
        # use log to get initial guess
        # empirical relation here from 10000 random draws
        # 10000 cases up to 85000 dtsec  0.11604050x + 9.69546575
        psiold = (math.log(dtsec) - 9.61202327)/0.10918231
        if psiold > upper:
            psiold = upper - math.pi
        # try arora-russel 2010 approach tbi(1 = psi, 2 = tof in sec)
        #         k1 = 1.0
        #         k2 = 5.712388980384687
        #         yx = magr1 + magr2 - k1*vara
        #         taux = math.sqrt(yx)/math.sqrt(mu) * (k2*yx + vara)
    else:
        #if (de == 'L')  # dm == 's' seemed better other way?????????????
        psiold = lower + (upper - lower)*0.5
        # else
        #            psiold = lower + (upper - lower)*0.6
        #         end
        # try arora-russel 2010 approach tbi(1 = psi, 2 = tof in sec)
        #         if dm == 'S'
        #             psiaux = 0.3 * tbi(nrev, 1) + 0.74 * nrev**2 * math.pi**2
        #         else
        #             psiaux = 0.3 * tbi(nrev, 1) + 0.74 * (nrev + 1)**2 * math.pi**2
        #         end
        #         [c2new, c3new] = smu.findc2c3(psiaux)
        #         y = magr1 + magr2 - (vara*(1.0-psiaux*c3new) / math.sqrt(c2new))
        #         taux = (xold**3*c3new + vara*math.sqrt(y)) * oomu
        #         a = (math.log(taux) - math.log(tbi(nrev, 2))) / (psiaux - tbi(nrev, 1))
        #         b = -2.0 * a * tbi(nrev, 1)
        #         c = math.log(tbi(nrev, 2)) + a * tbi(nrev, 1)**2
        #         a = a /math.log(dtsec)
        #         b = b /math.log(dtsec)
        #         c = c /math.log(dtsec)
        #         discrim = b*b - 4.0 *a*c
        #         if (discrim > 0.0)
        #             psiold1 = (-b + math.sqrt(discrim)) / (2.0 *a)
        #                psiold = (-b - math.sqrt(discrim)) / (2.0 *a)
        #psiold = lower + (upper - lower)*0.3

    c2new, c3new = smu.findc2c3(psiold)

    oosqrtmu = 1.0 / math.sqrt(mu)

    # find initial dtold from psiold
    print(magr1, magr2, vara, psiold, c3new)
    print(1.0 - psiold * c3new)
    print(vara * (1.0 - psiold * c3new))
    if (abs(c2new) > small):
        y = magr1 + magr2 - (vara * (1.0 - psiold * c3new) / math.sqrt(c2new))
    else:
        y = magr1 + magr2
    print(y, c2new)

    if ((vara > 0.0) and (y < 0.0)):  # (vara > 0.0) &
        ynegktr = 1
        while ((y < 0.0) and (ynegktr < 10)):
            psinew = 0.8*(1.0 / c3new)*(1.0 - (magr1 + magr2)
                                        *math.sqrt(c2new)/vara)
            # -------- find c2 and c3 functions -----------
            c2new, c3new = smu.findc2c3(psinew)
            psiold = psinew
            lower = psiold
            if (abs(c2new) > small):
                y = magr1 + magr2 - (vara*(1.0-psiold*c3new)
                                    / math.sqrt(c2new))
            else:
                y = magr1 + magr2
            # outfile
            print('y %11.7f lower %11.7f c2 %11.7f psinew %11.7f yneg %3i \n'
                    % (y, lower, c2new, psinew, ynegktr))
            ynegktr = ynegktr + 1

    if (abs(c2new) > small):
        xold = math.sqrt(y / c2new)
    else:
        xold = 0.0
    xoldcubed = xold * xold * xold
    dtold = (xoldcubed * c3new + vara * math.sqrt(y)) * oosqrtmu

    # -------  determine if  the orbit is possible at all ---------
    if (abs(vara) > 0.2):   # not exactly zero
        loops = 0
        ynegktr = 1  # y neg ktr
        dtnew = -10.0

        while ((abs(dtnew-dtsec) >= small)
               and (loops < numiter) and (ynegktr <= 10)):
            # print('%3i  dtnew-dtsec %11.7f yneg %3i \n', loops, dtnew-dtsec, ynegktr)
            if (abs(c2new) > small):
                y = magr1 + magr2 - (vara*(1.0 - psiold*c3new)/math.sqrt(c2new))
            else:
                y = magr1 + magr2
            # ----------- check for negative values of y ----------
            if ((vara > 0.0) and (y < 0.0)):  # (vara > 0.0) &
                ynegktr = 1
                while ((y < 0.0) and (ynegktr < 10)):
                    psinew = 0.8*(1.0 / c3new)*(1.0 - (magr1 + magr2)
                                                *math.sqrt(c2new)/vara)
                    # -------- find c2 and c3 functions -----------
                    c2new, c3new = smu.findc2c3(psinew)
                    psiold = psinew
                    lower = psiold
                    if (abs(c2new) > small):
                        y = magr1 + magr2 - (vara*(1.0-psiold*c3new)
                                            / math.sqrt(c2new))
                    else:
                        y = magr1 + magr2
                    # outfile
                    print('yneg %3i  y %11.7f lower %11.7f c2 %11.7f psinew %11.7f yneg %3i \n'
                          % (loops, y, lower, c2new, psinew, ynegktr))
                    ynegktr = ynegktr + 1
                #end # while
            #end  # if  y neg

            loops = loops + 1

            if (ynegktr < 10):
                if (abs(c2new) > small):
                    xold = math.sqrt(y / c2new)
                else:
                    xold = 0.0
                xcubed = xold**3
                dtnew = (xcubed*c3new + vara*math.sqrt(y)) * oomu

                # try newton rhapson iteration to update psi
                if (abs(psiold) > 1e-5):
                    c2dot = 0.5/psiold * (1.0 - psiold*c3new - 2.0*c2new)
                    c3dot = 0.5/psiold * (c2new - 3.0*c3new)
                else:  # case for parabolic orbit
                    c2dot = (-1.0/np.math.factorial(4) + 2.0*psiold/np.math.factorial(6)
                             - 3.0*psiold**2/np.math.factorial(8)
                             + 4.0*psiold**3/np.math.factorial(10)
                             - 5.0*psiold**4/np.math.factorial(12))
                    c3dot = (-1.0/np.math.factorial(5) + 2.0*psiold/np.math.factorial(7)
                             - 3.0*psiold**2/np.math.factorial(9)
                             + 4.0*psiold**3/np.math.factorial(11)
                             - 5.0*psiold**4/np.math.factorial(13))
                dtdpsi = (xcubed*(c3dot - 3.0*c3new*c2dot/(2.0*c2new))
                          + 0.125*vara
                          * (3.0*c3new*math.sqrt(y)/c2new + vara/xold)) * oomu
                # Newton iteration test to see if it keeps within the bounds
                psinew = psiold - (dtnew - dtsec)/dtdpsi

                # check if newton guess for psi is outside bounds (too steep a slope)
                if ((psinew > upper) or (psinew < lower)):
                    # --------  readjust upper and lower bounds -------
                    if ((de == 'L') or (nrev == 0)):
                        if (dtold < dtsec):
                            lower = psiold
                        else:
                            upper = psiold
                    else:
                        if (dtold < dtsec):
                            upper = psiold
                        else:
                            lower = psiold
                    psinew = (upper+lower) * 0.5
                    psilast = psinew

                # ------------- find c2 and c3 functions ----------
                c2new, c3new = smu.findc2c3(psinew)
                #if nrev > 0
                print('%3i  y %11.7f x %11.7f %11.7f dtnew %11.7f %11.7f %11.7f psinew %11.7f %11.7f \n' %
                    (loops, y, xold, dtsec, dtnew, lower, upper, psinew, dtdpsi)) #(dtnew - dtsec)/dtdpsi)( # c2dot, c3dot
                #end
                psilast = psiold  # keep previous iteration
                psiold = psinew
                dtold = dtnew

                # --- make sure the first guess isn't too close ---
                if ((abs(dtnew - dtsec) < small) and (loops == 1)):
                    dtnew = dtsec - 1.0
            #end  # if  ynegktr < 10

            #              print('#3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n', loops, y, x, dtnew, psinew)
            #              print('%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n', loops, y/re, x/math.sqrt(re), dtnew/tusec, psinew)
            #              print('%3i  y %11.7f x %11.7f dtnew %11.7f psinew %11.7f \n', loops, y/re, x/math.sqrt(re), dtnew/60.0, psinew)
            ynegktr = 1
        #end # while loop

        if ((loops >= numiter) or (ynegktr >= 10)):
            errorl = 'gnotconv' + str(abs(dtnew - dtsec))
            if (ynegktr >= 10):
                errorl = 'y negati'
        else:
            # --- use f and g series to find velocity vectors -----
            f = 1.0 - y/magr1
            gdot = 1.0 - y/magr2
            g = 1.0 / (vara*math.sqrt(y/mu))  # 1 over g
            #print('%11.7f  %11.7f  %11.7f \n', f, gdot, g)

            #  fdot = math.sqrt(mu*y)*(-magr2-magr1 + y)/(magr1*magr2*vara)
            #  f*gdot - fdot*g
            for i in range(3):
                v1dv[i] = (r2[i] - f*r1[i])*g
                v2dv[i] = (gdot*r2[i] - r1[i])*g

        #end   # if  the answer has converged
    else:
        errorl = 'impos180'

        #call Battin...

        # do hohmann but in 3-d...
        # can't do bissection because w series is not accurate
        mum = 3.986004418e5   # 14 m3/s2
        atx = (mum*(dtsec/(1.0*math.pi))**2)**(1.0/3.0)  # 1pi since half period
        v1tmag = math.sqrt(2.0*mum/(magr1) - mum/atx)
        v2tmag = math.sqrt(2.0*mum/(magr2) - mum/atx)
        wx = np.cross(r1, v1)
        wxu = smu.unit(wx)
        v1dir = np.cross(r1, wxu) # get retro direc
        v2dir = np.cross(r2, wxu) # get retro direc
        v1diru = smu.unit(v1dir)
        v2diru = smu.unit(v2dir)
        v1t = -v1tmag*v1diru
        v2t = -v2tmag*v2diru
        print('%11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n'
              % (r1[0], r1[1], r1[2], v1t[0], v1t[1], v1t[2]))
        print('%11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n'
              % (r2[0], r2[1], r2[2], v2t[0], v2t[1], v2t[2]))

        v1dv = v1t # / velkmps
        v2dv = v2t # / velkmps

        tof = dtsec

        #             # use JGCD 2011 v34 n6 1925 to solve 180 deg case
        #             p = 2.0*magr1*magr2 / (magr1 + magr2)
        #             ecc = math.sqrt(1.0 - 4.0*magr1*magr2 / ((magr1+magr2)**2))
        #             dt = math.sqrt(math.pi**2/mu * (p / (1.0 - ecc**2))**3)
        #             dnu = acos(cosdeltanu)
        #      print('hodo %11.6f   %11.6f  %11.6f  %11.6f %14.10f \n', p, ecc, dt, dtsec, vara)
        #             if abs(dt-dtsec) < 160
        #                 [v1dv, v2dv] = lambhodograph(r1, v1, r2, p, ecc, dnu, dtsec)
        #             end
    #end  # if  var a > 0.0

    if (errorl != '      ok'):
        print("\n\n-----Error found in lambertu: ", errorl)
        print("\n\n")
        p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = sc.rv2coe (r1, v1dv)
        print('%10s %3i %3i %2s %2s %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f case %11.7f %11.7f %11.7f %11.7f %11.7f ' % \
            (errorl, loops, nrev, dm, de, dtnew, y, xold, v1dv[0], v1dv[1], v1dv[2], v2dv[0], v2dv[1], v2dv[2], lower, upper, psinew, dtdpsi, ecc)) #(dtnew - dtsec)/dtdpsi, ecc))  # c2dot, c3dot
        print('C%3i %3i %2s %2s %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f  %11.7f dnu %11.7f \n' % \
            (loops, nrev, dm, de, dtnew, magr1, magr2, vara, y, xold, psinew, math.acos(cosdeltanu)*180.0/math.pi))
    else:
        #fprintf(outfile, '#s \n', errorl)
        print('#s \n', errorl)


    return v1dv, v2dv, errorl


# ------------------------------------------------------------------------------
#
#                           function lambertb
#
#  this function solves lambert's problem using battins method. the method is
#    developed in battin (1987) and explained by Thompson 2018. it uses continued
#    fractions to speed the solution and has several parameters that are defined
#    differently than the traditional gaussian technique.
#
#  author        : david vallado                  719-573-2600   12 feb 2018
#
#  inputs          description                    range / units
#    r1          - ijk position vector 1          km
#    v1          - ijk velocity vector 1 needed for 180 deg transfer  km / s
#    r2          - ijk position vector 2          km
#    dm          - dir of motion (long, short)        'l', 's'
#                  this is really a period discriminator
#    df          - dir of flight (direct, retrograde) 'd', 'r'
#                  this is the inclination discriminator
#    nrev        - number of revs to complete     0, 1, ...
#    dtsec       - time between r1 and r2         s
#
#  outputs       :
#    v1t         - ijk velocity vector            km / s
#    v2t         - ijk velocity vector            km / s
#    error       - error flag                     'ok', ...
#
#  locals        :
#    i           - index
#    loops       -
#    u           -
#    b           -
#    sinv        -
#    cosv        -
#    rp          -
#    x           -
#    xn          -
#    y           -
#    l           -
#    m           -
#    cosdeltanu  -
#    sindeltanu  -
#    dnu         -
#
#  coupling      :
#    mag         - magnitude of a vector
#    arcsinh     - inverse hyperbolic sine
#    arccosh     - inverse hyperbolic cosine
#    sinh        - hyperbolic sine
#
#  references    :
#    vallado       2013, 493-497, ex 7-5
#    thompson      2018
#
# [v1dv, v2dv, errorb] = lambertb (r1, v1, r2, dm, df, nrev, dtsec)
# ------------------------------------------------------------------------------

def lambertb(r1=None, v1=None, r2=None, dm=None, df=None, nrev=None,
             dtsec=None):
    errorb = '      ok'
    y = 0.0
    k2 = 0.0
    u = 0.0
    v1dv = np.array([1000, 1000, 1000])
    v2dv = np.array([1000, 1000, 1000])
    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)
    cosdeltanu = np.dot(r1, r2) / (magr1 * magr2)
    # make sure it's not more than 1.0
    if (np.abs(cosdeltanu) > 1.0):
        cosdeltanu = 1.0 * np.sign(cosdeltanu)

    rcrossr = np.cross(r1, r2)
    magrcrossr = smu.mag(rcrossr)
    if df == 'd':
        sindeltanu = magrcrossr / (magr1 * magr2)
    else:
        sindeltanu = - magrcrossr / (magr1 * magr2)

    dnu = math.atan2(sindeltanu, cosdeltanu)
    # the angle needs to be positive to work for the long way
    if dnu < 0.0:
        dnu = 2.0 * np.pi + dnu

    # these are the same
    chord = np.sqrt(magr1 * magr1 + magr2 * magr2
                    - 2.0 * magr1 * magr2 * cosdeltanu)
    # chord = smu.mag(r2-r1)

    s = (magr1 + magr2 + chord) * 0.5
    ror = magr2 / magr1
    eps = ror - 1.0
    lam = 1.0 / s * np.sqrt(magr1 * magr2) * np.cos(dnu * 0.5)
    L = ((1.0 - lam) / (1.0 + lam)) ** 2
    m = 8.0 * mu * dtsec * dtsec / (s ** 3 * (1.0 + lam) ** 6)
    #        tan2w = 0.25*eps*eps / (sqrt(ror) + ror * (2.0 + sqrt(ror)))
#        rp = sqrt(magr1*magr2)*((cos(dnu*0.25))^2 + tan2w)
#        if (dnu < pi)
#            L = ((sin(dnu*0.25))^2 + tan2w) / ((sin(dnu*0.25))^2 + tan2w + cos(dnu*0.5))
#        else
#            L = ((cos(dnu*0.25))^2 + tan2w - cos(dnu*0.5)) / ((cos(dnu*0.25))^2 + tan2w)
#        end
#        m = mu * dtsec*dtsec / (8.0*rp*rp*rp)
# initial guess
    if (nrev > 0):
        xn = 1.0 + 4.0 * L
    else:
        xn = L

    #    lim1 = sqrt(m/L)
# alt approach for high energy (long way, retro multi-rev) case
    if (dm == 'l') and (nrev > 0):
        xn = 1e-20
        x = 10.0
        loops = 1
        while ((np.abs(xn - x) >= small) and (loops <= 20)):

            x = xn
            temp = 1.0 / (2.0 * (L - x * x))
            temp1 = np.sqrt(x)
            temp2 = (nrev * np.pi * 0.5 + np.arctan(temp1)) / temp1
            h1 = temp * (L + x) * (1.0 + 2.0 * x + L)
            h2 = temp * m * temp1 * ((L - x * x) * temp2 - (L + x))
            b = 0.25 * 27.0 * h2 / ((temp1 * (1.0 + h1)) ** 3)
            if b < 0.0:
                f = 2.0 * np.cos(1.0 / 3.0 * np.arccos(np.sqrt(b + 1.0)))
            else:
                A = (np.sqrt(b) + np.sqrt(b + 1.0)) ** (1.0 / 3.0)
                f = A + 1.0 / A
            y = 2.0 / 3.0 * temp1 * (1.0 + h1) * (np.sqrt(b + 1.0) / f + 1.0)
            xn = 0.5 * ((m / (y * y) - (1.0 + L)) - np.sqrt((m / (y * y) - (1.0 + L)) ** 2 - 4.0 * L))
            print(' %3i yh %11.6f x %11.6f h1 %11.6f h2 %11.6f b %11.6f f %11.7f \n'
                  % (loops, y, x, h1, h2, b, f))
            loops = loops + 1

        print(' %3i yh %11.6f x %11.6f h1 %11.6f h2 %11.6f b %11.6f f %11.7f \n'
              % (loops, y, x, h1, h2, b, f))
        x = xn
        a = s * (1.0 + lam) ** 2 * (1.0 + x) * (L + x) / (8.0 * x)
        p = ((2.0 * magr1 * magr2 * (1.0 + x) * np.sin(dnu * 0.5) ** 2)
             / (s * (1 + lam) ** 2 * (L + x)))
        ecc = np.sqrt(1.0 - p / a)
        v1dv, v2dv = lambhodograph(r1, v1, r2, p, ecc, dnu, dtsec)
        print('high v1t %16.8f %16.8f %16.8f \n' % (v1dv))
    else:
        # standard processing
# note that the dr nrev =0 case is not represented
        loops = 1
        y1 = 0.0
        x = 10.0
        while ((np.abs(xn - x) >= small) and (loops <= 30)):

            if (nrev > 0):
                x = xn
                temp = 1.0 / ((1.0 + 2.0 * x + L) * (4.0 * x ** 2))
                temp1 = (nrev * np.pi * 0.5 + np.arctan(np.sqrt(x))) / np.sqrt(x)
                h1 = (temp * (L + x) ** 2
                      * (3.0 * (1.0 + x) ** 2 * temp1 - (3.0 + 5.0 * x)))
                h2 = (temp * m * ((x * x - x * (1.0 + L) - 3.0 * L)
                                  * temp1 + (3.0 * L + x)))
            else:
                x = xn
                tempx = smu.seebatt(x)
                denom = 1.0 / ((1.0 + 2.0 * x + L) * (4.0 * x + tempx * (3.0 + x)))
                h1 = (L + x) ** 2 * (1.0 + 3.0 * x + tempx) * denom
                h2 = m * (x - L + tempx) * denom
            # ----------------------- evaluate cubic ------------------
            b = 0.25 * 27.0 * h2 / ((1.0 + h1) ** 3)
            #        if b < -1.0 # reset the initial condition
#fprintf(1, 'xx #11.6f  #11.6f  #11.6f \n', L, xn, b)
#            xn = 1.0 - 2.0*L
#        end
#        else
#            if y1 > lim1
#                xn = xn * (lim1/y1)
#            end
#            else
            u = 0.5 * b / (1.0 + np.sqrt(1.0 + b))
            k2 = smu.kbat(u)
            y = (((1.0 + h1) / 3.0) * (2.0 + np.sqrt(1.0 + b)
                                       / (1.0 + 2.0 * u * k2 * k2)))
            xn = (np.sqrt(((1.0 - L) * 0.5) ** 2 + m / (y * y))
                  - (1.0 + L) * 0.5)
            #                    xn = sqrt(l*l + m/(y*y)) - (1.0 - l) alt, doesn't seem to work
#            end
#        end
            y1 = np.sqrt(m / ((L + x) * (1.0 + x)))
            loops = loops + 1
            print(' %3i yb %11.6f x %11.6f k2 %11.6f b %11.6f u %11.6f y1 %11.7f \n'
                  % (loops, y, x, k2, b, u, y1))

        print(' %3i yb %11.6f x %11.6f k2 %11.6f b %11.6f u %11.6f y1 %11.7f \n'
              % (loops, y, x, k2, b, u, y1))
        if (loops < 30):
            # blair approach use y from solution
#       lam = 1.0/s * sqrt(magr1*magr2) * cos(dnu*0.5)
#       m = 8.0*mu*dtsec*dtsec / (s^3*(1.0 + lam)^6)
#       L = ((1.0 - lam)/(1.0 + lam))^2
#a = s*(1.0 + lam)^2*(1.0 + x)*(lam + x) / (8.0*x)
# p = (2.0*magr1*magr2*(1.0 + x)*sin(dnu*0.5)^2)^2 / (s*(1 + lam)^2*(lam + x))  # loechler, not right?
            p = (2.0 * magr1 * magr2 * y * y * (1.0 + x) ** 2
                 * np.sin(dnu * 0.5) ** 2) / (m * s * (1 + lam) ** 2)
            ecc = np.sqrt((eps ** 2 + 4.0 * magr2 / magr1
                           * np.sin(dnu * 0.5) ** 2 * ((L - x) / (L + x)) ** 2)
                           / (eps ** 2 + 4.0 * magr2 / magr1
                              * np.sin(dnu * 0.5) ** 2))
            v1dv, v2dv = lambhodograph(r1, v1, r2, p, ecc, dnu, dtsec)
            #            fprintf(1, 'oldb v1t #16.8f #16.8f #16.8f #16.8f\n', v1dv, smu.mag(v1dv))
#         r_180 = 0.001  # 1 meter
#         [v1dvh, v2dvh] = lambert_vel(r1, v1, r2, dnu, p, ecc, mu, dtsec, r_180)
#         fprintf(1, 'newb v1t #16.8f #16.8f #16.8f #16.8f\n', v1dvh, smu.mag(v1dvh))
            # Battin solution to orbital parameters (and velocities)
# thompson 2011, loechler 1988
            if dnu > np.pi:
                lam = - np.sqrt((s - chord) / s)
            else:
                lam = np.sqrt((s - chord) / s)
            #      x = xn
            # loechler pg 21 seems correct!
            v1dvl = (1.0 / (lam * (1.0 + lam))
                     * np.sqrt(mu * (1.0 + x) / (2.0 * s ** 3 * (L + x)))
                     * ((r2 - r1) + s * (1.0 + lam) ** 2 * (L + x)
                        / (magr1 * (1.0 + x)) * r1))
            # added v2
            v2dvl = (1.0 / (lam * (1.0 + lam))
                     * np.sqrt(mu * (1.0 + x) / (2.0 * s ** 3 * (L + x)))
                     * ((r2 - r1) - s * (1.0 + lam) ** 2 * (L + x)
                        / (magr2 * (1.0 + x)) * r2))
            #fprintf(1, 'loe v1t #16.8f #16.8f #16.8f #16.8f\n', v1dvl, smu.mag(v1dvl))
#fprintf(1, 'loe v2t #16.8f #16.8f #16.8f #16.8f\n', v2dvl, smu.mag(v2dvl))

    return v1dv, v2dv, errorb

#    test case
#    ro = [ 2.500000    0.000000    0.000000]*6378.137
#    r = [ 1.9151111   1.6069690   0.000000]*6378.137
#
# lambert minimum energy, not min time!!
# min time is approximated by the parabolic case which is sort of a limit of
# what could be done
#
#  inputs          description                    range / units
#    r1          - ijk position vector 1          km
#    r2          - ijk position vector 2          km
#    dm          - direction of motion            'l', 's'
#    nrev        - multiple revoluions            0, 1, ...
#
#  outputs       :
#    v           - ijk velocity vector            km / s
#    aminenergy  - semimajor axis min energy      km / s
#    tminenergy  - time min energy                km / s
#    tminabs     - time min energy - parabolic    km / s
#
#

def lambertmin(r1=None, r2=None, dm=None, nrev=None):
    mu = 398600.4418

    # ---- find parameters that are constant for the initial geometry
    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)
    cosdeltanu = np.dot(r1, r2) / (magr1 * magr2)
    c = np.sqrt(magr1 ** 2 + magr2 ** 2 - 2.0 * magr1 * magr2 * cosdeltanu)
    s = 0.5 * (magr1 + magr2 + c)
    aminenergy = 0.5 * s
    alphae = np.pi
    betae = 2.0 * np.arcsin(np.sqrt((s - c) / s))
    if (dm == 'L'):
        tminenergy = (np.sqrt(aminenergy ** 3 / mu)
                      * (2.0 * nrev * np.pi + alphae - (betae - np.sin(betae))))
        #tminabs1 = sqrt(s^3/(8.0*mu))*(alphae - (betae-sin(betae)))
    else:
        tminenergy = (np.sqrt(aminenergy ** 3 / mu)
                      * (2.0 * nrev * np.pi + alphae + (betae - np.sin(betae))))
        #tminabs1 = sqrt(s^3/(8.0*mu))*(alphae + (betae-sin(betae)))

    # find parabolic tof - this will be the minimum limit for tof
# negative sign should be smallest
    tminabs = 1.0 / 3.0 * np.sqrt(2.0 / mu) * (s ** 1.5 - (s - c) ** 1.5)
    #tminabs = 1.0/3.0*sqrt(2.0/mu)*(s^1.5+(s-c)^1.5)

    # if calc min velocity
    rcrossr = np.cross(r1, r2)
    magrcrossr = smu.mag(rcrossr)
    pmin = magr1 * magr2 / c * (1.0 - cosdeltanu)
    if dm == 'S':
        sindeltanu = magrcrossr / (magr1 * magr2)
    else:
        sindeltanu = - magrcrossr / (magr1 * magr2)

    v = np.zeros(3)
    for i in range(3):
        v[i] = (np.sqrt(mu * pmin) / (magr1 * magr2 * sindeltanu)
                * (r2[i] - (1.0 - magr2 / pmin * (1.0 - cosdeltanu)) * r1[i]))

    return v, aminenergy, tminenergy, tminabs

    return v, aminenergy, tminenergy, tminabs


#------------------------------------------------------------------------------
#
#                           procedure lambertminT
#
#  this procedure solves lambert's problem and finds the miniumum time for
#  multi-revolution cases.
#
#  author        : david vallado                  719-573-2600   22 mar 2018
#
# inputs          description                    range / units
#    r1          - ijk position vector 1          km
#    r2          - ijk position vector 2          km
#    dm          - direction of motion            'l', 's'
#    df          - direction of flight            'd', 'r'
#    nrev        - number of revs to complete     0, 1, 2, 3, ...
#
#  outputs       :
#    tmin        - minimum time of flight         sec
#    tminp       - minimum parabolic tof          sec
#    tminenergy  - minimum energy tof             sec
#
#  locals        :
#    i           - index
#    loops       -
#    cosdeltanu  -
#    sindeltanu  -
#    dnu         -
#    chord       -
#    s           -
#
#  coupling      :
#    mag         - magnitude of a vector
#    dot         - dot product
#
#  references    :
#    vallado       2013, 494, Alg 59, ex 7-5
#    prussing      JAS 2000
#-----------------------------------------------------------------------------*/

def lambertminT(r1=None, r2=None, dm=None, de=None, nrev=None):
    mu = 398600.4415
    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)
    cosdeltanu = np.dot(r1, r2) / (magr1 * magr2)
    # make sure it's not more than 1.0
    if (np.abs(cosdeltanu) > 1.0):
        cosdeltanu = 1.0 * np.sign(cosdeltanu)

    rcrossr = np.cross(r1, r2)
    if (de == 'L'):
        sindeltanu = smu.mag(rcrossr) / (magr1 * magr2)
    else:
        sindeltanu = - smu.mag(rcrossr) / (magr1 * magr2)

    dnu = math.atan2(sindeltanu, cosdeltanu)
    # the angle needs to be positive to work for the long way
    if (dnu < 0.0):
        dnu = 2.0 * np.pi + dnu

    # ------------- try prussing test his numbers are wrong -------
#dnu = 75.0 / (180.0 / pi)
#cosdeltanu = cos(dnu)
#magr1 = 1.0
#magr2 = 1.524
#mu = 4.0 * pi * pi
#dnu = 90.0 / (180.0 / pi)
#cosdeltanu = cos(dnu)
#magr1 = 1.0  # er
#magr2 = 1.0
#mu = 1.0

    # these are the same
#    if (de == 'L')
    chord = np.sqrt(magr1 * magr1 + magr2 * magr2
                    - 2.0 * magr1 * magr2 * cosdeltanu)
    #    else
#    chord = -sqrt(magr1 * magr1 + magr2 * magr2 - 2.0 * magr1 * magr2 * cosdeltanu)
#    end
#chord = smu.mag(r2 - r1)

    s = (magr1 + magr2 + chord) * 0.5
    xi = 0.0
    eta = 0.0
    # ------------- calc tmin parabolic tof to see if the orbit is possible
# ----- no ellitpical orbits exist below this --------
    if (dm == 'S'):
        tminp = ((1.0 / 3.0) * np.sqrt(2.0 / mu)
                 * ((s ** 1.5) - (s - chord) ** 1.5))
    else:
        tminp = ((1.0 / 3.0) * np.sqrt(2.0 / mu)
                 * ((s ** 1.5) + (s - chord) ** 1.5))

    # could do this just for nrev cases, but you can also find these for any nrev if (nrev > 0)
# ------------- this is the min energy ellipse tof
    amin = 0.5 * s
    #alpha = pi
    beta = 2.0 * np.arcsin(np.sqrt((s - chord) / s))
    if (dm == 'S'):
        tminenergy = ((amin ** 1.5)
                      * ((2.0 * nrev + 1.0) * np.pi - beta + np.sin(beta))
                      / np.sqrt(mu))
    else:
        tminenergy = ((amin ** 1.5)
                      * ((2.0 * nrev + 1.0) * np.pi + beta - np.sin(beta))
                      / np.sqrt(mu))

    # -------------- calc min tof ellipse (prussing 1992 aas, 2000 jas, stern 1964 pg 230)
    an = 1.001 * amin
    fa = 10.0
    i = 1
    while (np.abs(fa) > 1e-05 and i <= 20):

        a = an
        alp = 1.0 / a
        alpha = 2.0 * np.arcsin(np.sqrt(0.5 * s * alp))
        if (de == 'L'):
            beta = 2.0 * np.arcsin(np.sqrt(0.5 * (s - chord) * alp))
        else:
            beta = - 2.0 * np.arcsin(np.sqrt(0.5 * (s - chord) * alp))
        xi = alpha - beta
        eta = np.sin(alpha) - np.sin(beta)
        fa = ((6.0 * nrev * np.pi + 3.0 * xi - eta) * (np.sin(xi) + eta)
              - 8.0 * (1.0 - np.cos(xi)))
        fadot = (((6.0 * nrev * np.pi + 3.0 * xi - eta)
                  * (np.cos(xi) + np.cos(alpha)) + (3.0 - np.cos(alpha))
                  * (np.sin(xi) + eta) - 8.0 * np.sin(xi))
                 * (- alp * np.tan(0.5 * alpha))
                 + ((6.0 * nrev * np.pi + 3.0 * xi - eta)
                    * (- np.cos(xi) - np.cos(alpha)) + (- 3.0 - np.cos(beta))
                    * (np.sin(xi) + eta) + 8.0 * np.sin(xi))
                 * (- alp * np.tan(0.5 * beta)))
        del_ = fa / fadot
        an = a - del_
        #        fprintf(1, '#2i #8.4f #11.5f #11.5f  #11.5f  #11.5f  #11.5f  #11.5f  #11.5f \n', i, dnu * rad2deg, alpha * rad2deg, beta * rad2deg, xi, eta, fa, fadot, an)
        i = i + 1


    print('iter %2i ' % (i))
    # could update beta one last time with alpha too????
    if (dm == 'S'):
        tmin = (an ** 1.5) * (2.0 * np.pi * nrev + xi - eta) / np.sqrt(mu)
    else:
        tmin = (an ** 1.5) * (2.0 * np.pi * nrev + xi + eta) / np.sqrt(mu)

    #      dm = 'S'
#      de = 'L'
#      nrev = 0
#      blair
#      r1 = [6778.136300000, 0.000000, 0.000000 ]
#      r2 = [-6694.857334274, -1180.483980026, 0.000000 ]
#      v1 = [0.000000, 7.668558568, 0.000000 ]
#      moving
#      r1 = [ -6175.1034, 2757.0706, 1626.6556 ]
#      v1 = [ 2.376641, 1.139677, 7.078097]
#      r2 = [ -1078.007289, 8796.641859, 1890.7135 ]

    #      [tmin, tminp, tminenergy] = lambertminT(r1, r2, dm, de, nrev)

    print('%c  %c %i  %f  %f  %f  \n' % (dm, de, nrev, tmin, tminp, tminenergy))
    return tmin, tminp, tminenergy

# -------------------------------------------------------------------------
# find the minimum psi values for the universal variable lambert problem
# for multi-rev cases
#  inputs          description                    range / units
#    r1          - ijk position vector 1             km
#    r2          - ijk position vector 2             km
#    dm          - direction of motion (long/short)  'L', 'S'
#    you could maybe pass in the max number of revs you're willing to
#    consider...
#    nrev        - multiple revoluions                0, 1, ...
#
#  outputs       :
#    tbi         - 2D array containing the dt (sec), psi (unitless) values for the minimums of
#    nrev = 1, 2, 3, ...
#

def lambertumins(r1=None, r2=None, nrev=None, dm=None):
    small = 1e-08
    mu = 398600.4418

    oomu = 1.0 / np.sqrt(mu)

    sqrtmu = np.sqrt(mu)
    numiter = 20

    # ---- find parameters that are constant for the intiial geometry
    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)
    cosdeltanu = np.dot(r1, r2) / (magr1 * magr2)
    if (dm == 'L'):
        vara = - np.sqrt(magr1 * magr2 * (1.0 + cosdeltanu))
    else:
        vara = np.sqrt(magr1 * magr2 * (1.0 + cosdeltanu))

    # ------------ find the minimum time for a nrev transfer --------------
#   nrev = 0
#   for zz = 0: 4
#       nrev = nrev + 1
# ---- get outer bounds for each nrev case
    lower = 4.0 * nrev ** 2 * np.pi * np.pi
    upper = 4.0 * (nrev + 1.0) ** 2 * np.pi * np.pi
    # ---- streamline since we know it's near the center
    upper = lower + (upper - lower) * 0.6
    lower = lower + (upper - lower) * 0.3
    # ---- initial estimate, just put in center
    psiold = (upper + lower) * 0.5
    c2, c3 = smu.findc2c3(psiold)
    loops = 0
    dtdpsi = 200.0
    while ((np.abs(dtdpsi) >= 0.1) and (loops < numiter)):

        if (np.abs(c2) > small):
            y = magr1 + magr2 - (vara * (1.0 - psiold * c3) / np.sqrt(c2))
        else:
            y = magr1 + magr2
        if (np.abs(c2) > small):
            x = np.sqrt(y / c2)
        else:
            x = 0.0
        sqrty = np.sqrt(y)
        if np.abs(psiold) > 1e-05:
            c2dot = 0.5 / psiold * (1.0 - psiold * c3 - 2.0 * c2)
            c3dot = 0.5 / psiold * (c2 - 3.0 * c3)
            c2ddot = (1.0 / (4.0 * psiold ** 2)
                      * ((8.0 - psiold) * c2 + 5.0 * psiold * c3 - 4.0))
            c3ddot = (1.0 / (4.0 * psiold ** 2)
                      * ((15.0 - psiold) * c3 - 7.0 * c2 + 1.0))
        else:
            c2dot = (- 2.0 / np.math.factorial(4) + 2.0 * psiold / np.math.factorial(6)
                     - 3.0 * psiold ** 2 / np.math.factorial(8)
                     + 4.0 * psiold ** 3 / np.math.factorial(10)
                     - 5.0 * psiold ** 4 / np.math.factorial(12))
            c3dot = (- 1.0 / np.math.factorial(5) + 2.0 * psiold / np.math.factorial(7)
                     - 3.0 * psiold ** 2 / np.math.factorial(9)
                     + 4.0 * psiold ** 3 / np.math.factorial(11)
                     - 5.0 * psiold ** 4 / np.math.factorial(13))
            c2ddot = 0.0
            c3ddot = 0.0
        # now solve this for dt = 0.0
#            dtdpsi = x^3*(c3dot - 3.0*c3*c2dot/(2.0*c2))* oomu + 0.125*vara/sqrt(mu) * (3.0*c3*sqrty/c2 + vara/x)
        dtdpsi = (x ** 3 * (c3dot - 3.0 * c3 * c2dot / (2.0 * c2))
                  * oomu + 0.125 * vara
                  * (3.0 * c3 * sqrty / c2 + vara / x) * oomu)
        q = 0.25 * vara * np.sqrt(c2) - x ** 2 * c2dot
        s1 = - 24.0 * q * x ** 3 * c2 * sqrty * c3dot
        s2 = (36.0 * q * x ** 3 * sqrty * c3 * c2dot
              - 16.0 * x ** 5 * sqrty * c3ddot * c2 ** 2)
        s3 = (24.0 * x ** 5 * sqrty
              * (c3dot * c2dot * c2 + c3 * c2ddot * c2 - c3 * c2dot ** 2)
              - 6.0 * vara * c3dot * y * c2 * x ** 2)
        s4 = (- 0.75 * vara ** 2 * c3 * c2 ** 1.5 * x ** 2
              + 6.0 * vara * c3 * y * c2dot * x ** 2
              + (vara ** 2 * c2 * (0.25 * vara * np.sqrt(c2) - x ** 2 * c2))
              * sqrty / x)
        dtdpsi2 = (- (s1 + s2 + s3 + s4)
                   / (16.0 * sqrtmu * (c2 ** 2 * sqrty * x ** 2)))
        # NR update
        psinew = psiold - dtdpsi / dtdpsi2
        #           fprintf(1, ' #3i #12.4f #12.4f #12.4f #12.4f #11.4f #12.4f #12.4f #11.4f #11.4f \n', loops, y, dtdpsi, psiold, psinew, psinew - psiold, dtdpsi, dtdpsi2, lower, upper)
        psiold = psinew
        c2, c3 = smu.findc2c3(psiold)
        # don't need for the loop iterations
# dtnew = (x^3*c3 + vara*sqrty) * oomu
        loops = loops + 1


    # calculate once at the end
    dtnew = (x ** 3 * c3 + vara * sqrty) * oomu
    tof = dtnew
    psib = psinew
    print(' %3i %12.4f %12.4f %12.4f %12.4f %11.4f %12.4f %12.4f %11.4f %11.4f \n'
          % (loops, y, dtdpsi, psiold, psinew, psinew - psiold, dtdpsi, dtdpsi2, lower, upper))
    #       fprintf(1, ' nrev #3i  dtnew #12.5f psi #12.5f  lower #10.3f upper #10.3f #10.6f #10.6f \n', nrev, dtnew, psiold, lower, upper, c2, c3)
#   end # for checking multi rev cases

    return psib, tof


#
# form the minimum time and universal variable matrix
#
# [tbiu, tbilu] = lambgettbiu(r1, r2, 3)

def lambgettbiu(r1=None, r2=None, order=None):
    tbi = np.zeros((order, 2))
    #tbi = [0 0 0 0 0 0 0 0 0 0]
    for i in range(order):
        psib, tof = lambertumins(r1, r2, i + 1, 'S')
        tbi[i, 0] = psib
        tbi[i, 1] = tof

    tbil = np.zeros((order, 2))
    for i in range(order):
        psib, tof = lambertumins(r1, r2, i + 1, 'L')
        tbil[i, 0] = psib
        tbil[i, 1] = tof

    return tbi, tbil

# ------------------------------------------------------------------------------
#
#                           function lambhodograph
#
# this function accomplishes 180 deg transfer (and 360 deg) for lambert problem.
#
#  author        : david vallado                  719-573-2600   22 may 2017
#
#  inputs          description                    range / units
#    r1          - ijk position vector 1          km
#    r2          - ijk position vector 2          km
#    dtsec       - time between r1 and r2         s
#    dnu         - true anomaly change            rad
#
#  outputs       :
#    v1t         - ijk transfer velocity vector   km / s
#    v2t         - ijk transfer velocity vector   km / s
#
#  references    :
#    Thompson JGCD 2013 v34 n6 1925
#    Thompson AAS GNC 2018
# [v1t, v2t] = lambhodograph(r1, v1, r2, p, a, ecc, dnu, dtsec)
# ------------------------------------------------------------------------------

def lambhodograph(r1=None, v1=None, r2=None, p=None,
                  ecc=None, dnu=None, dtsec=None):
    mu = 398600.4418

    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)
    eps = 0.001 / magr2

    a = mu * (1.0 / magr1 - 1.0 / p)

    b = (mu * ecc / p) ** 2 - a ** 2
    if b <= 0.0:
        x1 = 0.0
    else:
        x1 = - np.sqrt(b)

    # 180 deg, and multiple 180 deg transfers
    if np.abs(np.sin(dnu)) < eps:
        nvec = np.cross(r1, v1) / smu.mag(np.cross(r1, v1))
        if ecc < 1.0:
            ptx = twopi * np.sqrt(p ** 3 / (mu * (1.0 - ecc ** 2) ** 3))
            if (np.mod(dtsec, ptx) > ptx * 0.5):
                x1 = - x1
        print('less than\n' % ())
    else:
        # more common path?
        y2a = mu / p - x1 * np.sin(dnu) + a * np.cos(dnu)
        y2b = mu / p + x1 * np.sin(dnu) + a * np.cos(dnu)
        if np.abs(mu / magr2 - y2b) < np.abs(mu / magr2 - y2a):
            x1 = - x1
        # depending on the cross product, this will be normal or in plane,
# or could even be a fan
        nvec = np.cross(r1, r2) / smu.mag(np.cross(r1, r2))
        if (np.mod(dnu, twopi) > np.pi):
            nvec = - nvec
        #       fprintf(1, 'gtr than\n')

    v1t = (np.sqrt(mu * p) / magr1) * ((x1 / mu) * r1 + np.cross(nvec, r1) / magr1)
    x2 = x1 * np.cos(dnu) + a * np.sin(dnu)
    v2t = (np.sqrt(mu * p) / magr2) * ((x2 / mu) * r2 + np.cross(nvec, r2) / magr2)
    return v1t, v2t




# ------------------------------------------------------------------------------
#
#                           function kepler
#
#  this function solves keplers problem for orbit determination and returns a
#    future geocentric equatorial (ijk) position and velocity vector.  the
#    solution uses universal variables.
#
#  author        : david vallado                  719-573-2600   22 jun 2002
#
#  revisions
#    vallado     - fix some mistakes                             13 apr 2004
#
#  inputs          description                    range / units
#    ro          - ijk position vector - initial  km
#    vo          - ijk velocity vector - initial  km / s
#    dtsec       - length of time to propagate    s
#
#  outputs       :
#    r           - ijk position vector            km
#    v           - ijk velocity vector            km / s
#    error       - error flag                     'ok', ...
#
#  locals        :
#    f           - f expression
#    g           - g expression
#    fdot        - f dot expression
#    gdot        - g dot expression
#    xold        - old universal variable x
#    xoldsqrd    - xold squared
#    xnew        - new universal variable x
#    xnewsqrd    - xnew squared
#    znew        - new value of z
#    c2new       - c2(psi) function
#    c3new       - c3(psi) function
#    dtsec       - change in time                 s
#    timenew     - new time                       s
#    rdotv       - result of ro dot vo
#    a           - semi or axis                   km
#    alpha       - reciprocol  1/a
#    sme         - specific mech energy           km2 / s2
#    period      - time period for satellite      s
#    s           - variable for parabolic case
#    w           - variable for parabolic case
#    h           - angular momentum vector
#    temp        - temporary real*8 value
#    i           - index
#
#  coupling      :
#    mag         - magnitude of a vector
#    findc2c3    - find c2 and c3 functions
#
#  references    :
#    vallado       2004, 95-103, alg 8, ex 2-4
#
# [r, v] = kepler  (ro, vo, dtsec)
# ------------------------------------------------------------------------------


def kepler(ro, vo, dtseco):

  # -------------------------  implementation   -----------------
  # set constants and intermediate printouts
  velkmps = math.sqrt(mu / re)
  numiter = 50

  if sh.show =='y':
      print(' ro %16.8f %16.8f %16.8f km \n'% (ro[0], ro[1], ro[2]))
      print(' vo %16.8f %16.8f %16.8f km/s \n'% (vo[0], vo[1], vo[2]))
      print(' ro %16.8f %16.8f %16.8f ER \n'
            % (ro[0]/re, ro[1]/re, ro[2]/re))
      print(' vo %16.8f %16.8f %16.8f ER/TU \n'
            % (vo[0]/velkmps, vo[1]/velkmps, vo[2]/velkmps))

  # --------------------  initialize values   -------------------
  ktr = 0
  xold = 0.0
  znew = 0.0
  errork = '      ok'
  dtsec = dtseco
  nrev = 0
  sqmu = math.sqrt(mu)

  if (abs(dtseco) > small):
      magro = smu.mag(ro)
      magvo = smu.mag(vo)
      rdotv = np.dot(ro, vo)

      # -------------  find sme, alpha, and a  ------------------
      sme = ((magvo**2)*0.5) - (mu /magro)
      alpha = -sme*2.0/mu

      if (abs(sme) > small):
          a = -mu / (2.0 *sme)
      else:
          a = infinite
      if (abs(alpha) < small):   # parabola
          alpha = 0.0

      if sh.show =='y':
          print(' sme %16.8f  a %16.8f alp  %16.8f ER \n'
                % (sme/(mu/re), a/re, alpha*re))
          print(' sme %16.8f  a %16.8f alp  %16.8f km \n' % (sme, a, alpha))
          print(' ktr      xn        psi           r          xn+1        dtn \n')

      # ------------   setup initial guess for x  ---------------
      # -----------------  circle and ellipse -------------------
      if (alpha >= small):
          period = twopi * math.sqrt(abs(a)**3.0/mu)
          # ------- next if needed for 2body multi-rev ----------
          if (abs(dtseco) > abs(period)):
              # including the truncation will produce vertical lines that are parallel
              # (plotting chi vs time)
              dtsec = math.fmod(dtseco, period)
              nrev = math.floor(dtseco/period)
          xold = sqmu *dtsec * alpha
      else:
          # --------------------  parabola  ---------------------
          if (abs(alpha) < small):
              h = np.cross(ro, vo)
              magh = smu.mag(h)
              p = magh*magh/mu
              s = 0.5 * (halfpi - math.atan(3.0 * math.sqrt(mu / (p*p*p)) * dtsec))
              w = math.atan(math.tan(s) ** (1.0 / 3.0))
              xold = math.sqrt(p) * (2.0 * math.cot(2.0 * w))
              alpha = 0.0
          else:
              # ------------------  hyperbola  ------------------
              temp = -2.0 * mu * dtsec / \
                  (a*(rdotv + np.sign(dtsec)*math.sqrt(-mu*a)* \
                  (1.0 -magro*alpha)))
              xold = np.sign(dtsec) * math.sqrt(-a) *math.log(temp)

      ktr = 1
      dtnew = -10.0
      # conv for dtsec to x units
      tmp = 1.0 / sqmu

      while ((abs(dtnew*tmp - dtsec) >= small) and (ktr < numiter)):
          xoldsqrd = xold*xold
          znew = xoldsqrd * alpha

          # ------------- find c2 and c3 functions --------------
          c2new, c3new = smu.findc2c3(znew)

          # ------- use a newton iteration for new values -------
          rval = xoldsqrd*c2new + rdotv*tmp *xold*(1.0 -znew*c3new) + \
              magro*(1.0  - znew*c2new)
          dtnew = xoldsqrd*xold*c3new + rdotv*tmp*xoldsqrd*c2new + \
              magro*xold*(1.0  - znew*c3new)

          # ------------- calculate new value for x -------------
          temp1 = (dtsec*sqmu - dtnew) / rval
          xnew = xold + temp1

          # ----- check if the univ param goes negative. if so, use bissection
          if (xnew < 0.0 and dtsec > 0.0):
              xnew = xold*0.5

          if sh.show =='y':
              print('kep %3i %11.7f %11.7f %11.7f %11.7f %11.7f \n' \
                  % (ktr, xold, znew, rval, xnew, dtnew*tmp))
              print('%3i %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f %11.7f \n' \
                  % (ktr, xold/math.sqrt(re), znew, c2new, c3new, rval/re, xnew/math.sqrt(re), dtnew/math.sqrt(mu)))

          ktr = ktr + 1
          xold = xnew

      if (ktr >= numiter):
          errork = 'knotconv'
          print('kep not conv in %2i iter %11.3f \n' % (numiter, dtseco))
          r = np.array([0.0, 0.0, 0.0])
          v = np.array([0.0, 0.0, 0.0])
          for i in range(3):
              v[i] = 0.0
              r[i] = v[i]
      else:
          # --- find position and velocity vectors at new time --
          xnewsqrd = xnew*xnew
          f = 1.0  - (xnewsqrd*c2new / magro)
          g = dtsec - xnewsqrd*xnew*c3new/math.sqrt(mu)

          r = np.array([0.0, 0.0, 0.0])
          v = np.array([0.0, 0.0, 0.0])
          for i in range(3):
              r[i] = f*ro[i] + g*vo[i]
          magr = smu.mag(r)
          gdot = 1.0  - (xnewsqrd*c2new / magr)
          fdot = (math.sqrt(mu)*xnew / (magro*magr)) * (znew*c3new-1.0)
          for i in range(3):
              v[i] = fdot*ro[i] + gdot*vo[i]

          temp = f*gdot - fdot*g
          if (abs(temp-1.0) > 0.00001):
              errork = 'fandg'

          if sh.show =='y':
              print('f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n' % (f, g, fdot, gdot))
              tusec = math.sqrt(re**3/mu)
              print('f %16.8f g %16.8f fdot %16.8f gdot %16.8f \n' % (f, g/tusec, fdot*tusec, gdot))
              print('r1 %16.8f %16.8f %16.8f km \n' % (r[0], r[1], r[2]))
              print('v1 %16.8f %16.8f %16.8f km/s \n' % (v[0], v[1], v[2]))
              print('r1 %16.8f %16.8f %16.8f ER \n'
                    % (r[0]/re, r[1]/re, r[2]/re))
              print('v1 %16.8f %16.8f %16.8f ER/TU \n'
                    % (vo[0]/velkmps, vo[1]/velkmps, vo[2]/velkmps))

  else:
      # ----------- set vectors to incoming since 0 time --------
      r = np.array([0.0, 0.0, 0.0])
      v = np.array([0.0, 0.0, 0.0])
      for i in range(3):
          r[i] = ro[i]
          v[i] = vo[i]


  return r, v, errork




# ------------------------------------------------------------------------------
#
#                           function pkepler
#
#  this function propagates a satellite's position and velocity vector over
#    a given time period accounting for perturbations caused by j2.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    ro          - original position vector       km
#    vo          - original velocity vector       km/sec
#    ndot        - time rate of change of n       rad/sec
#    nddot       - time accel of change of n      rad/sec2
#    dtsec       - change in time                 sec
#
#  outputs       :
#    r           - updated position vector        km
#    v           - updated velocity vector        km/sec
#
#  locals        :
#    p           - semi-paramter                  km
#    a           - semior axis                    km
#    ecc         - eccentricity
#    incl        - inclination                    rad
#    argp        - argument of periapsis          rad
#    argpdot     - change in argument of perigee  rad/sec
#    raan       - longitude of the asc node      rad
#    raandot    - change in raan                rad
#    e0          - eccentric anomaly              rad
#    e1          - eccentric anomaly              rad
#    m           - mean anomaly                   rad/sec
#    mdot        - change in mean anomaly         rad/sec
#    arglat      - argument of latitude           rad
#    arglatdot   - change in argument of latitude rad/sec
#    truelon     - true longitude of vehicle      rad
#    truelondot  - change in the true longitude   rad/sec
#    lonper     - longitude of periapsis         rad
#    lonperodot  - longitude of periapsis change  rad/sec
#    n           - mean angular motion            rad/sec
#    nuo         - true anomaly                   rad
#    j2op2       - j2 over p sqyared
#    sinv, cosv   - sine and cosine of nu
#
#  coupling:
#    rv2coe      - orbit elements from position and velocity vectors
#    coe2rv      - position and velocity vectors from orbit elements
#    newtonm     - newton rhapson to find nu and eccentric anomaly
#
#  references    :
#    vallado       2007, 687, alg 64
#
# [r, v] = pkepler(ro, vo, dtsec, ndot, nddot)
# ----------------------------------------------------------------------------- }


def pkepler(ro, vo, dtsec, ndot, nddot):

    p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper = sc.rv2coe(ro, vo)
    print(p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper)
    print('          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg\n')
    print('coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n' % \
            (p, a, ecc, incl * rad2deg, raan * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg))
    print("arglat, truelon, lonper: ", arglat, truelon, lonper)

    n = math.sqrt(mu/(a*a*a))

     # ------------- find the value of j2 perturbations -------------
    j2op2 = (n*1.5*re**2*j2) / (p*p)
#     nbar = n*(1.0 + j2op2*sqrt(1.0-ecc*ecc)* (1.0 - 1.5*sin(incl)*sin(incl)))
    raandot = -j2op2 * math.cos(incl)
    argpdot = j2op2 * (2.0-2.5*math.sin(incl)*math.sin(incl))
    mdot = n

    a = a - 2.0*ndot*dtsec * a / (3.0*n)
    ecc = ecc - 2.0*(1.0 - ecc)*ndot*dtsec / (3.0*n)
    p = a*(1.0 - ecc*ecc)

     # ----- update the orbital elements for each orbit type --------
    if ecc < small:
        # -------------  circular equatorial  ----------------
        if (incl < small) or (abs(incl-math.pi) < small):
            truelondot = raandot + argpdot + mdot
            truelon = truelon  + truelondot * dtsec
            truelon = math.fmod(truelon, twopi)
        else:
          # -------------  circular inclined    --------------
            raan = raan + raandot * dtsec
            raan = math.fmod(raan, twopi)
            arglatdot = argpdot + mdot
            arglat = arglat + arglatdot * dtsec
            arglat = math.fmod(arglat, twopi)
          # -- elliptical, parabolic, hyperbolic equatorial ---
    elif (incl < small) or (abs(incl-math.pi) < small):
        lonperdot = raandot + argpdot
        lonper = lonper + lonperdot * dtsec
        lonper = math.fmod(lonper, twopi)
        m = m + mdot*dtsec + ndot*dtsec*dtsec + nddot*dtsec*dtsec*dtsec
        m = math.fmod(m, twopi)
        e0, nu = smu.newtonm(ecc, m)
    else:
        # --- elliptical, parabolic, hyperbolic inclined --
        raan = raan + raandot * dtsec
        raan = math.fmod(raan, twopi)
        argp = argp  + argpdot  * dtsec
        argp = math.fmod(argp, twopi)
        m = m + mdot*dtsec + ndot*dtsec*dtsec + nddot*dtsec*dtsec*dtsec
        m = np.remainder(m, twopi)
        e0, nu = smu.newtonm(ecc, m)

     # ------------- use coe2rv to find new vectors ---------------
    r, v = sc.coe2rv(p, ecc, incl, raan, argp, nu, arglat, truelon, lonper)
    r = r.T
    v = v.T
#        print('r    %15.9f%15.9f%15.9f' % (r[0], r[1], r[2]))
#        print(' v %15.10f%15.10f%15.10f\n' % ([0], v[1], v[2]))

    return r, v


# ------------------------------------------------------------------------------
#
#                           function pkeplerj4
#
#  this function propagates a satellite's position and velocity vector over
#    a given time period accounting for perturbations caused by j2, j2^2, and j4.
#    to use without drag, set ndot and nddot to 0.0 as inputs. the "xxxe1"
#    variables are a mostly common denominator form of the escobal equations
#
#  author        : david vallado                  719-573-2600    1 aug 2018
#
#  inputs          description                    range / units
#    ro          - original position vector       km
#    vo          - original velocity vector       km/sec
#    ndot        - time rate of change of n       rad/sec
#    nddot       - time accel of change of n      rad/sec2
#    dtsec       - change in time                 sec
#
#  outputs       :
#    r           - updated position vector        km
#    v           - updated velocity vector        km/sec
#
#  locals        :
#    p           - semi-paramter                  km
#    a           - semior axis                    km
#    ecc         - eccentricity
#    incl        - inclination                    rad
#    argp        - argument of periapsis          rad
#    argpdot     - change in argument of perigee  rad/sec
#    raan       - longitude of the asc node      rad
#    raandot    - change in raan                rad
#    e0          - eccentric anomaly              rad
#    e1          - eccentric anomaly              rad
#    m           - mean anomaly                   rad/sec
#    mdot        - change in mean anomaly         rad/sec
#    arglat      - argument of latitude           rad
#    arglatdot   - change in argument of latitude rad/sec
#    truelon     - true longitude of vehicle      rad
#    truelondot  - change in the true longitude   rad/sec
#    lonper     - longitude of periapsis         rad
#    lonperodot  - longitude of periapsis change  rad/sec
#    n           - mean angular motion            rad/sec
#    nuo         - true anomaly                   rad
#    j2op2       - j2 over p sqyared
#    sinv, cosv   - sine and cosine of nu
#
#  coupling:
#    rv2coe      - orbit elements from position and velocity vectors
#    coe2rv      - position and velocity vectors from orbit elements
#    newtonm     - newton rhapson to find nu and eccentric anomaly
#
#  references    :
#    vallado       2007, 687, alg 64
#
# [r, v] = pkeplerj4(ro, vo, dtsec, ndot, nddot)
# ----------------------------------------------------------------------------- }

def pkeplerj4(ro=None, vo=None, dtsec=None, ndot=None, nddot=None):
    # egm-08

    j2 = 0.00108262617
    j3 = - 2.53241052e-06
    j4 = - 1.6198976e-06
    j6 = - 5.40666576e-07
    p, a, ecc, incl, raan, argp, nu, m, arglat, truelon, lonper = sc.rv2coe(ro, vo)
    #     fprintf(1, '          p km       a km      ecc      incl deg     raan deg     argp deg      nu deg      m deg      arglat   truelon    lonper\n')
#     fprintf(1, 'coes #11.4f#11.4f#13.9f#13.7f#11.5f#11.5f#11.5f#11.5f#11.5f#11.5f#11.5f\n', ...
#             p, a, ecc, incl * rad2deg, raan * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg, ...
#             arglat * rad2deg, truelon * rad2deg, lonper * rad2deg)
    n = np.sqrt(mu / (a * a * a))
    cosargp = np.cos(argp)
    sinargp = np.sin(argp)
    cosi = np.cos(incl)
    sini = np.sin(incl)
    cosraan = np.cos(raan)
    sinraan = np.sin(raan)
    beta2 = (1.0 - ecc ** 2)
    sqrtbeta = np.sqrt(beta2)
    # ------------- find the value of j2 perturbations -------------
# merson approach
#     nndot = 0.75*j2*n*(re/p)^2*sqrtbeta*(2.0 - 3.0*sini^2)...
#         +3.0/512.0*j2^2*n*(re/p)^4 * sqrtbeta *(320.0*ecc*ecc - 280.0*ecc^4 + (1600.0 - 1568.0*ecc*ecc + 328*ecc^4)*sini*sini ...
#             + (-2096.0 + 1072.0*ecc*ecc + 79.0*ecc^4)*sini^4)...
#         -45.0/128.0*j4*n*ecc*ecc*(re/p)^4*sqrtbeta*(-8.0 + 40.0*sini - 35.0*sini*sini)... # should be 8.0 - 40.0*sini^2 + 35.0*sini^4
#         +35.0/2048.0*j6*n*(re/p)^6*sqrtbeta*(-128.0 + 320.0*ecc^2 + 240.0*ecc^4 + ...
#         (1344.0 - 3360.0*ecc^2 - 2520.0*ecc^4)*sini + (-1512.0 + 3780.0*ecc^2 + 2835.0*ecc^4)*sini^2 ...
#         - (-1848.0 + 4620.0*ecc^2 + 3465.0*ecc^4)*sini^4)
#     mdot = n + nndot

    # from escobal pg 371 (the e terms) seems best (closest to STK)
    nbar = (n * (1.5 * j2 * (re / p) ** 2 * sqrtbeta * ((1.0 - 1.5 * sini ** 2))
                 + 3.0 / 128.0 * j2 ** 2 * (re / p) ** 4 * sqrtbeta
                 * (16.0 * sqrtbeta + 25.0 * beta2
                    - 15.0 + (30.0 - 96.0 * sqrtbeta - 90.0 * beta2) * cosi
                    * cosi + (105.0 + 144.0 * sqrtbeta + 25.0 * beta2)
                    * cosi ** 4)
                 - 45.0 / 128.0 * j4 * ecc * ecc * (re / p) ** 4
                 * sqrtbeta * (3.0 - 30.0 * cosi * cosi + 35.0 * cosi ** 4)))
    mdot = n + nbar
    # merson approach
#     raandot = -1.5*n*j2*(re/p)^2*cosi...
#         +3.0/32.0*n*j2^2*(re/p)^4*cosi*(12.0 - 4.0*ecc^2 - sini^2*(80.0 + 5.0*ecc*ecc))...
#         +15.0/32.0*n*j4*(re/p)^4*cosi*(8.0 + 12.0 * ecc*ecc - sini^2*(14.0 + 21.0*ecc*ecc)) ...
#         +105.0/1024.0*j6*n*(re/p)^6*cosi*(64.0 + 160.0*ecc^2 + 120.0*ecc^4 + ...
#         -(288.0 + 720.0*ecc^2 + 540.0*ecc^4)*sini^2 + (264.0 + 660.0*ecc^2 + 495.0*ecc^4)*sini^4)
    raandot = (-1.5 * j2 * (re / p) ** 2 * mdot * cosi
               - 1.5 * j2 * (re / p) ** 2 * mdot * cosi * 1.5 * j2
               * (re / p) ** 2
               * (1.5 + ecc ** 2 / 6.0 - 2.0 * sqrtbeta
                  - (5.0 / 3.0 - 5.0 / 24.0 * ecc * ecc - 3.0 * sqrtbeta)
                  * sini * sini)
               - 35.0 / 8.0 * j4 * (re / p) ** 4 * n * cosi
               * (1.0 + 1.5 * ecc * ecc) * ((12.0 - 21.0 * sini ** 2) / 14.0))
    # common denominators same!
    raandote1 = (-1.5 * j2 * (re / p) ** 2 * mdot * cosi
                 - 9.0 / 96.0 * j2 ** 2 * (re / p) ** 4 * mdot * cosi
                 * (36.0 + 4.0 * ecc ** 2 - 48.0 * sqrtbeta
                    - (40.0 - 5.0 * ecc * ecc - 72.0 * sqrtbeta) * sini * sini)
                 - 35.0 / 112.0 * j4 * (re / p) ** 4 * n * cosi
                 * (1.0 + 1.5 * ecc * ecc) * (12.0 - 21.0 * sini ** 2))
    # merson approach
#     argpdot = 0.75*n*j2*(re/p)^2*(4.0 - 5.0*sini^2) ...
#         + 9.0/384.0*n*j2^2*(re/p)^4*(56.0*ecc*ecc + (760.0 - 36.0*ecc*ecc)*sini^2 - (890.0 + 45.0*ecc*ecc)*sini^4) ...
#         -15.0/128.0*n*j4*(re/p)^4*(64.0 + 72.0*ecc*ecc - (248.0 + 252.0*ecc*ecc)*sini*sini + (196.0 + 189.0*ecc*ecc)*sini^4) ...
#         +105.0/2048.0*j6*n*(re/p)^6*(256.0 + 960.0*ecc^2 + 320.0*ecc^4 + ...
#         -(2048.0 + 6880.0*ecc^2 + 2160.0*ecc^4)*sini^2 + (4128.0 + 13080.0*ecc^2 + 3960.0*ecc^4)*sini^4 ...
#         -(2376.0 + 14520.0*ecc^2 + 2145.0*ecc^4)*sini^6)
    argpdot = (1.5 * j2 * (re / p) ** 2 * mdot * (2.0 - 2.5 * sini ** 2)
               + 9.0 / 4.0 * j2 ** 2 * (re / p) ** 4 * mdot
               * (2.0 - 2.5 * sini ** 2)
               * (2.0 + ecc * ecc / 2.0 - 2.0 * sqrtbeta
                  - (43.0 / 24.0 - ecc * ecc / 48.0 - 3.0 * sqrtbeta)
                  * sini * sini)
               - 45.0 / 36.0 * j2 ** 2 * (re / p) ** 4 * ecc * ecc * n
               * cosi ** 4 - 35.0 / 8.0 * j4 * (re / p) ** 4 * n
               * (12.0 / 7.0 - 93.0 / 14.0 * sini ** 2 + 21.0 / 4.0 * sini ** 4
                  + ecc * ecc * (27.0 / 14.0 - 189.0 / 28.0
                                 * sini * sini + 81.0 / 16.0 * sini ** 4)))
    # same!
    argpdote1 = (0.75 * j2 * (re / p) ** 2 * mdot * (4.0 - 5.0 * sini ** 2)
                 + 9.0 / 192.0 * j2 ** 2 * (re / p) ** 4 * mdot
                 * (2.0 - 2.5 * sini ** 2)
                 * (96.0 + 24.0 * ecc * ecc - 96.0 * sqrtbeta
                    - (86.0 - ecc * ecc - 144.0 * sqrtbeta) * sini * sini)
                 - 45.0 / 36.0 * j2 ** 2 * (re / p) ** 4 * ecc * ecc * n
                 * cosi ** 4 - 35.0 / 896.0 * j4 * (re / p) ** 4 * n
                 * (192.0 - 744.0 * sini ** 2 + 588 * sini ** 4 + ecc * ecc
                    * (216.0 - 756.0 * sini * sini + 567.0 * sini ** 4)))
    raandot = raandote1
    argpdot = argpdote1
    print('diffs   %15.11g  %15.11g  %15.11f  %15.11f \n' %
          (raandote1 - raandot, argpdote1 - argpdot, nddot - nbar, mdot - mdot))
    # can put these back in if estimates of ndot etc are known
#    a = a - 2.0*ndot*dtsec * a / (3.0*n)
#    ecc = ecc - 2.0*(1.0 - ecc)*ndot*dtsec / (3.0*n)
#    p = a*(1.0 - ecc*ecc)

    # ----- update the orbital elements for each orbit type --------
    if ecc < small:
        # -------------  circular equatorial  ----------------
        if (incl < small) or (np.abs(incl - np.pi) < small):
            truelondot = raandot + argpdot + mdot
            truelon = truelon + truelondot * dtsec
            truelon = np.fmod(truelon, twopi)
            print('circ equat\n' % ())
        else:
            # -------------  circular inclined    --------------
            raan = raan + raandot * dtsec
            raan = np.fmod(raan, twopi)
            arglatdot = argpdot + mdot
            arglat = arglat + arglatdot * dtsec
            arglat = np.fmod(arglat, twopi)
            print('circ incl\n' % ())
    else:
        # -- elliptical, parabolic, hyperbolic equatorial ---
        if (incl < small) or (np.abs(incl - np.pi) < small):
            lonperdot = raandot + argpdot
            lonper = lonper + lonperdot * dtsec
            lonper = np.fmod(lonper, twopi)
            m = m + mdot * dtsec + ndot * dtsec * dtsec \
                + nddot * dtsec * dtsec * dtsec
            m = np.fmod(m, twopi)
            e0, nu = smu.newtonm(ecc, m)
            print('other equat\n' % ())
        else:
            # --- elliptical, parabolic, hyperbolic inclined --
            raan = raan + raandot * dtsec
            raan = np.fmod(raan, twopi)
            argp = argp + argpdot * dtsec
            argp = np.fmod(argp, twopi)
            m = m + mdot * dtsec + ndot * dtsec * dtsec \
                + nddot * dtsec * dtsec * dtsec
            m = np.fmod(m, twopi)
            e0, nu = smu.newtonm(ecc, m)

    # ------------- use coe2rv to find new vectors ---------------
    r, v = sc.coe2rv(p, ecc, incl, raan, argp, nu, arglat, truelon, lonper)
    r = r.T
    v = v.T
    #        fprintf(1, 'r    #15.9f#15.9f#15.9f', r)
#        fprintf(1, ' v #15.10f#15.10f#15.10f\n', v)
    return r, v



# ------------------------------------------------------------------------------
#
#                           function anglesdr
#
#  this function solves the problem of orbit determination using three
#    optical sightings.  the solution function uses the double-r technique.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#   8 oct 2007
#
#  inputs          description                               range / units
#    rtasc1       - right ascension #1                          rad
#    rtasc2       - right ascension #2                          rad
#    rtasc3       - right ascension #3                          rad
#    decl1        - declination #1                              rad
#    decl2        - declination #2                              rad
#    decl3        - declination #3                              rad
#    jd1, jdf1    - julian date of 1st sighting                 days from 4713 bc
#    jd2, jdf2    - julian date of 2nd sighting                 days from 4713 bc
#    jd3, jdf3    - julian date of 3rd sighting                 days from 4713 bc
#    rsite1       - eci site position vector                    km
#    rsite2       - eci site position vector                    km
#    rsite3       - eci site position vector                    km
#
#  outputs        :
#    r2           - eci position vector at t2                   km
#    v2           - eci velocity vector at t2                   km / s
#
#  locals         :
#    l1           - line of sight vector for 1st
#    l2           - line of sight vector for 2nd
#    l3           - line of sight vector for 3rd
#    tau          - taylor expansion series about
#                   tau (t - to)
#    tausqr       - tau squared
#    i            - index
#    d            -
#    rho          - range from site to sat at t2                km
#    rhodot       -
#    dmat         -
#    earthrate    - velocity of earth rotation
#    p            -
#    q            -
#    oldr         -
#    oldv         -
#    f1           - f coefficient
#    g1           -
#    f3           -
#    g3           -
#    l2dotrs      -
#
#  coupling       :
#    mag          - magnitude of a vector
#    matmult      - multiply two matrices together
#    angl         - angl between two vectors
#
#  references     :
#    vallado       2007, 439-443
#
# [r2, v2] = anglesdr (decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, rsite1, rsite2, rsite3, re, mu)
# ------------------------------------------------------------------------------

def anglesdr(decl1=None, decl2=None, decl3=None, rtasc1=None,
             rtasc2=None, rtasc3=None, jd1=None, jdf1=None,
             jd2=None, jdf2=None, jd3=None, jdf3=None,
             rsite1=None, rsite2=None, rsite3=None, re=None, mu=None):
    # -------------------------  implementation   -------------------------
#   constastro

    # for sun
#re = 149597870.0
#mu = 1.32712428e11
    magr1in = 2.0 * re
    magr2in = 2.04 * re
    direct = 'y'
    tol = 1e-08 * re

    pctchg = 0.005
    # subtract dates and convert fraction of day to seconds
    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0

    tau13 = (jd1 - jd3) * 86400.0 + (jdf1 - jdf3) * 86400.0
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0
    # form line of sight vectors
    los1 = np.transpose(np.array([np.cos(decl1) * np.cos(rtasc1), np.cos(decl1)
                                  * np.sin(rtasc1), np.sin(decl1)]))
    los2 = np.transpose(np.array([np.cos(decl2) * np.cos(rtasc2), np.cos(decl2)
                                  * np.sin(rtasc2), np.sin(decl2)]))
    los3 = np.transpose(np.array([np.cos(decl3) * np.cos(rtasc3), np.cos(decl3)
                                  * np.sin(rtasc3), np.sin(decl3)]))
    # --------- now we're ready to start the actual double r algorithm ---------
    magr1old = 99999.9
    magr2old = 99999.9
    magrsite1 = smu.mag(rsite1)
    magrsite2 = smu.mag(rsite2)
    magrsite3 = smu.mag(rsite3)
    # take away negatives because escobal defines rs opposite
    cc1 = 2.0 * np.dot(los1, rsite1)
    cc2 = 2.0 * np.dot(los2, rsite2)
    ktr = 0
    # main loop to get three values of the double-r for processing
    while (np.abs(magr1in - magr1old) > (tol or np.abs(magr2in - magr2old)) > tol):

        ktr = ktr + 1
        print('%2i ' % (ktr))
        r2, r3, f1, f2, q1, magr1, magr2, a, deltae32 = \
            smu.doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2,
                        los3, rsite1, rsite2, rsite3, tau12, tau32, direct, re, mu)
        # check intermediate status
        f = 1.0 - a / magr2 * (1.0 - np.cos(deltae32))
        g = tau32 - np.sqrt(a ** 3 / mu) * (deltae32 - np.sin(deltae32))
        v2 = (r3 - f * r2) / g
        p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = \
            sc.rv2coeh(r2, v2, re, mu)
        ###print('coes %11.4f%11.4f%13.9f%13.7f%11.5f%11.5f%11.5f%11.5f\n' % (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg))
        # -------------- re-calculate f1 and f2 with r1 = r1 + delta r1
        magr1o = magr1in
        deltar1 = pctchg * magr1in
        magr1in = magr1in + deltar1
        r2, r3, f1delr1, f2delr1, q2, magr1, magr2, a, deltae32 = \
            smu.doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2,
                        los3, rsite1, rsite2, rsite3, tau12, tau32, direct, re, mu)
        pf1pr1 = (f1delr1 - f1) / deltar1
        pf2pr1 = (f2delr1 - f2) / deltar1
        # ----------------  re-calculate f1 and f2 with r2 = r2 + delta r2
        magr1in = magr1o
        magr2o = magr2in
        deltar2 = pctchg * magr2in
        magr2in = magr2in + deltar2
        r2, r3, f1delr2, f2delr2, q3, magr1, magr2, a, deltae32 = \
            smu.doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1,
                        los2, los3, rsite1, rsite2, rsite3, tau12, tau32, direct, re, mu)
        pf1pr2 = (f1delr2 - f1) / deltar2
        pf2pr2 = (f2delr2 - f2) / deltar2
        #       f = 1.0 - a/magr2*(1.0-cos(deltae32))
#       g = t3 - sqrt(a^3/mu)*(deltae32-sin(deltae32))
#       v2 = (r3 - f*r2)/g
#       [p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r2, v2)
#       fprintf(1, 'coes #11.4f#11.4f#13.9f#13.7f#11.5f#11.5f#11.5f#11.5f\n', ...
#              p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg)
        # ------------ now calculate an update
        magr2in = magr2o
        delta = pf1pr1 * pf2pr2 - pf2pr1 * pf1pr2
        delta1 = pf2pr2 * f1 - pf1pr2 * f2
        delta2 = pf1pr1 * f2 - pf2pr1 * f1
        deltar1 = - delta1 / delta
        deltar2 = - delta2 / delta
        magr1old = magr1in
        magr2old = magr2in
        #  may need to limit the amount of the correction
        if np.abs(deltar1) > magr1in * pctchg:
            print('%11.7f \n' % (deltar1))
            #         deltar1 = sign(deltar1)*magr1in*pctchg
        if np.abs(deltar2) > magr2in * pctchg:
            print('%11.7f \n' % (deltar2))
            #         deltar2 = sign(deltar2)*magr2in*pctchg
        magr1in = magr1in + deltar1
        magr2in = magr2in + deltar2
        print('qs %11.7f  %11.7f  %11.7f \n' % (q1, q2, q3))
        print('magr1o %11.7f delr1 %11.7f magr1 %11.7f %11.7f  \n'
              % (magr1o, deltar1, magr1in, magr1old))
        print('magr2o %11.7f delr2 %11.7f magr2 %11.7f %11.7f  \n'
              % (magr2o, deltar2, magr2in, magr2old))
        #       f = 1.0 - a/magr2*(1.0-cos(deltae32))
#       g = t3 - sqrt(a^3/mu)*(deltae32-sin(deltae32))
#       v2 = (r3 - f*r2)/g
#       [p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper ] = rv2coe (r2, v2)
#       fprintf(1, 'coes #11.4f#11.4f#13.9f#13.7f#11.5f#11.5f#11.5f#11.5f\n', ...
#              p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg)
        print('=============================================== \n' % ())

    # needed to get the r2 set properly since the last one was moving r2
    r2, r3, f1, f2, q1, magr1, magr2, a, deltae32 = \
        smu.doubler(cc1, cc2, magrsite1, magrsite2, magr1in, magr2in, los1, los2, los3,
                    rsite1, rsite2, rsite3, tau12, tau32, direct, re, mu)
    if a and mu and magr2:
      f = 1.0 - a / magr2 * (1.0 - np.cos(deltae32))
      g = tau32 - np.sqrt(a ** 3 / mu) * (deltae32 - np.sin(deltae32))
      v2 = (r3 - f * r2) / g
    else:
      v2=None
    return r2, v2

# ------------------------------------------------------------------------------
#
#                           function anglesg
#
#  this function solves the problem of orbit determination using three
#    optical sightings.  the solution function uses the gaussian technique.
#    there are lots of debug statements in here to test various options.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  23 dec 2003
#   8 oct 2007
#
#  inputs          description                    range / units
#    re           - radius earth, sun, etc        km
#    mu           - grav param earth, sun etc     km3/s2
#    rtasc1       - right ascension #1                          rad
#    rtasc2       - right ascension #2                          rad
#    rtasc3       - right ascension #3                          rad
#    decl1        - declination #1                              rad
#    decl2        - declination #2                              rad
#    decl3        - declination #3                              rad
#    jd1, jdf1    - julian date of 1st sighting                 days from 4713 bc
#    jd2, jdf2    - julian date of 2nd sighting                 days from 4713 bc
#    jd3, jdf3    - julian date of 3rd sighting                 days from 4713 bc
#    rsite1       - eci site position vector                    km
#    rsite2       - eci site position vector                    km
#    rsite3       - eci site position vector                    km
#
#  outputs        :
#    r            - ijk position vector at t2     km
#    v            - ijk velocity vector at t2     km / s
#
#  locals         :
#    l1           - line of sight vector for 1st
#    l2           - line of sight vector for 2nd
#    l3           - line of sight vector for 3rd
#    tau          - taylor expansion series about
#                   tau (t - to)
#    tausqr       - tau squared
#    t21t23       - (t2-t1) * (t2-t3)
#    t31t32       - (t3-t1) * (t3-t2)
#    i            - index
#    d            -
#    rho          - range from site to sat at t2  km
#    rhodot       -
#    dmat         -
#    rs1          - site vectors
#    rs2          -
#    rs3          -
#    earthrate    - velocity of earth rotation
#    p            -
#    q            -
#    oldr         -
#    oldv         -
#    f1           - f coefficient
#    g1           -
#    f3           -
#    g3           -
#    l2dotrs      -
#
#  coupling       :
#    mag          - magnitude of a vector
#    detrminant   - evaluate the determinant of a matrix
#    factor       - find roots of a polynomial
#    matmult      - multiply two matrices together
#    gibbs        - gibbs method of orbit determination
#    hgibbs       - herrick gibbs method of orbit determination
#    angl         - angl between two vectors
#
#  references     :
#    vallado       2007, 429-439
#
# [r2, v2] = anglesg (decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, rs1, rs2, rs3, re, mu)
# ------------------------------------------------------------------------------

def anglesg(decl1=None, decl2=None, decl3=None, rtasc1=None,
            rtasc2=None, rtasc3=None, jd1=None, jdf1=None,
            jd2=None, jdf2=None, jd3=None, jdf3=None, rs1=None,
            rs2=None, rs3=None):
    # -------------------------  implementation   -------------------------
    ddpsi = 0.0

    ddeps = 0.0
    magr1in = 2.0 * re

    magr2in = 2.01 * re
    direct = 'y'

    # ---------- set middle to 0, find decls to others -----------
    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0

    tau13 = (jd1 - jd3) * 86400.0 + (jdf1 - jdf3) * 86400.0
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0
    print('jd123 %14.6f %14.6f %14.6f %14.6f %14.6f %14.6f  \n'
          % (jd1, jdf1, jd2, jdf2, jd3, jdf3))
    print('tau12 %14.6f tau13  %14.6f tau32  %14.6f \n' % (tau12, tau13, tau32))
    # ----------------  find line of sight unit vectors  ---------------
    los1 = np.zeros(3)
    los1[0] = np.cos(decl1) * np.cos(rtasc1)
    los1[1] = np.cos(decl1) * np.sin(rtasc1)
    los1[2] = np.sin(decl1)
    los2 = np.zeros(3)
    los2[0] = np.cos(decl2) * np.cos(rtasc2)
    los2[1] = np.cos(decl2) * np.sin(rtasc2)
    los2[2] = np.sin(decl2)
    los3 = np.zeros(3)
    los3[0] = np.cos(decl3) * np.cos(rtasc3)
    los3[1] = np.cos(decl3) * np.sin(rtasc3)
    los3[2] = np.sin(decl3)
    # topo to body fixed (ecef)
#     latgd = 40.0/rad
#     lon = -110.0/rad
#     l1
#     outv = rot2(l1, -pi + latgd)
#     l1 = rot3(outv, -lon)
#     outv = rot2(l2, -pi + latgd)
#     l2 = rot3(outv, -lon)
#     outv = rot2(l3, -pi + latgd)
#     l3 = rot3(outv, -lon)
#     l1

    #     # take the middle trans from eecef to eci
#     tm = [-0.830668624503591  -0.556765707115059   0.001258429966118...
#    0.556766123167794  -0.830669298539751  -0.000023583565998...
#    0.001058469658016   0.000681061045186   0.999999207898604]


    # l1 = tm*l1'
# l2 =tm*l2'
# l3 =tm*l3'

    # ------------- find l matrix and determinant -----------------
    los1
    vs = np.transpose(np.array([0, 0, 0]))
    aecef = np.transpose(np.array([0, 0, 0]))
    #[l1eci, vs3, aeci] = ecef2eci(l1', vs, aecef, (jd1-2451545.0)/36525.0, jd1, 0.0, 0.0, 0.0, 0, ddpsi, ddeps)
#[l2eci, vs3, aeci] = ecef2eci(l2', vs, aecef, (jd2-2451545.0)/36525.0, jd2, 0.0, 0.0, 0.0, 0, ddpsi, ddeps)
#[l3eci, vs3, aeci] = ecef2eci(l3', vs, aecef, (jd3-2451545.0)/36525.0, jd3, 0.0, 0.0, 0.0, 0, ddpsi, ddeps)

    l1eci = los1
    l2eci = los2
    l3eci = los3
    # leave these as they come since the topoc radec are already eci
    print(l1eci)
    # --------- called lmati since it is only used for determ -----
    lmat = np.zeros((3, 3))
    rsmat = np.zeros((3, 3))
    for i in range(3):
        lmat[i, 0] = l1eci[i]
        lmat[i, 1] = l2eci[i]
        lmat[i, 2] = l3eci[i]
        rsmat[i, 0] = rs1[i]
        rsmat[i, 1] = rs2[i]
        rsmat[i, 2] = rs3[i]

    print(lmat)
    rmt = rsmat.T / re
    print('rsmat eci %11.7f  %11.7f  %11.7f km \n'
          % (rmt[0, 0], rmt[0, 1], rmt[0, 2]))
    # the order is right, but to print out, need '
    print('rsmat eci %11.7f  %11.7f  %11.7f km \n'
          % (rmt[1, 0], rmt[1, 1], rmt[1, 2]))
    print('rsmat eci %11.7f  %11.7f  %11.7f km \n'
          % (rmt[2, 0], rmt[2, 1], rmt[2, 2]))
    print(lmat)
    print('this should be the inverse of what the code finds later\n' % ())
    li = lmat.T
    # alt way of Curtis not seem to work ------------------
    p1 = np.cross(los2, los3)
    p2 = np.cross(los1, los3)
    p3 = np.cross(los1, los2)
    # both are the same
    dx = np.dot(los1, p1)
    #dx = dot(los3, p3)
    lcmat = np.zeros((3, 3))
    lcmat[0, 0] = np.dot(rs1, p1)
    lcmat[1, 0] = np.dot(rs2, p1)
    lcmat[2, 0] = np.dot(rs3, p1)
    lcmat[0, 1] = np.dot(rs1, p2)
    lcmat[1, 1] = np.dot(rs2, p2)
    lcmat[2, 1] = np.dot(rs3, p2)
    lcmat[0, 2] = np.dot(rs1, p3)
    lcmat[1, 2] = np.dot(rs2, p3)
    lcmat[2, 2] = np.dot(rs3, p3)
    tau31 = (jd3 - jd1) * 86400.0
    print(lcmat)
    aa = 1 / dx * (- lcmat[0, 1] * tau32 / tau31 + lcmat[1, 1]
                   + lcmat[2, 1] * tau12 / tau31)
    bb = 1 / (6.0 * dx) * (lcmat[0, 1] * (tau32 * tau32 - tau31 * tau31)
                           * tau32 / tau31 + lcmat[2, 1]
                           * (tau31 * tau31 - tau12 * tau12)
                           * tau12 / tau31)
    # alt way of Curtis not seem to work ------------------

    d = np.linalg.det(lmat)
    print(d)
    lmati = np.zeros((3, 3))
    # ------------------ now assign the inverse -------------------
    lmati[0, 0] = (l2eci[1] * l3eci[2] - l2eci[2] * l3eci[1]) / d
    lmati[1, 0] = (- l1eci[1] * l3eci[2] + l1eci[2] * l3eci[1]) / d
    lmati[2, 0] = (l1eci[1] * l2eci[2] - l1eci[2] * l2eci[1]) / d
    lmati[0, 1] = (- l2eci[0] * l3eci[2] + l2eci[2] * l3eci[0]) / d
    lmati[1, 1] = (l1eci[0] * l3eci[2] - l1eci[2] * l3eci[0]) / d
    lmati[2, 1] = (- l1eci[0] * l2eci[2] + l1eci[2] * l2eci[0]) / d
    lmati[0, 2] = (l2eci[0] * l3eci[1] - l2eci[1] * l3eci[0]) / d
    lmati[1, 2] = (- l1eci[0] * l3eci[1] + l1eci[1] * l3eci[0]) / d
    lmati[2, 2] = (l1eci[0] * l2eci[1] - l1eci[1] * l2eci[0]) / d
    print(lmati)
    lir = lmati * rsmat
    # ------------ find f and g series at 1st and 3rd obs ---------
#   speed by assuming circ sat vel for udot here ??
#   some similartities in 1/6t3t1 ...
# --- keep separated this time ----
    a1 = tau32 / (tau32 - tau12)
    a1u = ((tau32 * ((tau32 - tau12) ** 2 - tau32 * tau32))
           / (6.0 * (tau32 - tau12)))
    a3 = - tau12 / (tau32 - tau12)
    a3u = (- (tau12 * ((tau32 - tau12) ** 2 - tau12 * tau12))
           / (6.0 * (tau32 - tau12)))
    print('a1/a3 %11.7f  %11.7f  %11.7f  %11.7f \n' % (a1, a1u, a3, a3u))
    # --- form initial guess of r2 ----
    dl1 = lir[1, 0] * a1 - lir[1, 1] + lir[1, 2] * a3
    dl2 = lir[1, 0] * a1u + lir[1, 2] * a3u
    print(dl1)
    print(dl2)
    # ------- solve eighth order poly not same as laplace ---------
    magrs2 = smu.mag(rs2)
    l2dotrs = np.dot(los2, rs2)
    print('magrs2 %11.7f  %11.7f  \n' % (magrs2, l2dotrs))
    poly = np.zeros(9)
    poly[0] = 1.0

    poly[1] = 0.0
    poly[2] = - (dl1 * dl1 + 2.0 * dl1 * l2dotrs + magrs2 ** 2)
    poly[3] = 0.0
    poly[4] = 0.0
    poly[5] = - 2.0 * mu * (l2dotrs * dl2 + dl1 * dl2)
    poly[6] = 0.0
    poly[7] = 0.0
    poly[8] = - mu * mu * dl2 * dl2
    print(poly)
    rootarr = np.roots(poly)
    print(rootarr)
    #fprintf(1, 'rootarr #11.7f \n', rootarr)

    # ------------------ select the correct root ------------------
    bigr2 = - 99999990.0
    # change from 1
    for j in range(8):
        if (rootarr[j] > bigr2) and np.isreal(rootarr[j]):
            bigr2 = rootarr[j]

    print(bigr2)
    # ------------ solve matrix with u2 better known --------------
    u = mu / (bigr2 * bigr2 * bigr2)
    c1 = a1 + a1u * u
    c2 = - 1.0
    c3 = a3 + a3u * u
    print('u %17.14f c1 %11.7f c3 %11.7f %11.7f \n' % (u, c1, c2, c3))
    cmat = np.zeros((3, 3))
    cmat[0, 0] = - c1
    cmat[1, 0] = - c2
    cmat[2, 0] = - c3
    rhomat = lir * cmat
    rhoold1 = rhomat[0, 0] / c1
    rhoold2 = rhomat[1, 0] / c2
    rhoold3 = rhomat[2, 0] / c3
    print('rhoold %11.7f %11.7f %11.7f \n' % (rhoold1, rhoold2, rhoold3))
    #   fprintf(1, 'rhoold #11.7f #11.7f #11.7f \n', rhoold1/re, rhoold2/re, rhoold3/re)

    r1 = np.zeros(3)
    r2 = np.zeros(3)
    r3 = np.zeros(3)
    for i in range(3):
        r1[i] = rhomat[0, 0] * l1eci[i] / c1 + rs1[i]
        r2[i] = rhomat[1, 0] * l2eci[i] / c2 + rs2[i]
        r3[i] = rhomat[2, 0] * l3eci[i] / c3 + rs3[i]

    print('r1 %11.7f %11.7f %11.7f \n' % (r1[0], r1[1], r1[2]))
    print('r2 %11.7f %11.7f %11.7f \n' % (r2[0], r2[1], r2[2]))
    print('r3 %11.7f %11.7f %11.7f \n' % (r3[0], r3[1], r3[2]))
    # -------- loop through the refining process ------------  while () for
    print('now refine the answer \n' % ())
    rho2 = infinite
    ll = 0
    while ((np.abs(rhoold2 - rho2) > 1e-12) and (ll <= 0)):

        ll = ll + 1
        print(' iteration %3i \n' % (ll))
        rho2 = rhoold2
        # ---------- now form the three position vectors ----------
        for i in range(3):
            r1[i] = rhomat[0, 0] * l1eci[i] / c1 + rs1[i]
            r2[i] = - rhomat[1, 0] * l2eci[i] + rs2[i]
            r3[i] = rhomat[2, 0] * l3eci[i] / c3 + rs3[i]
        magr1 = smu.mag(r1)
        magr2 = smu.mag(r2)
        magr3 = smu.mag(r3)
        v2, theta, theta1, copa, error = gibbsh(r1, r2, r3, re, mu)
        print('r1 %16.14f %16.14f %16.14f %11.7f %11.7f %16.14f \n'
              % (r1[0], r1[1], r1[2], theta * rad2deg, theta1 * rad2deg, copa * rad2deg))
        print('r2 %11.7f %11.7f %11.7f \n' % (r2[0], r2[1], r2[2]))
        print('r3 %11.7f %11.7f %11.7f \n' % (r3[0], r3[1], r3[2]))
        print('w gibbs km/s       v2 %11.7f %11.7f %11.7f \n'
              % (v2[0], v2[1], v2[2]))
        # check if too close obs
        if ((str(error) == str('          ok')) and
                ((np.abs(theta) < 1.0 * deg2rad) or (np.abs(theta1) < 1.0 * deg2rad))):
            p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = \
                sc.rv2coeh(r2, v2, re, mu)
            ###print('coes init ans %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n' % (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg))
            # --- hgibbs to get middle vector ----
            v2, theta, theta1, copa, error = \
                hgibbs(r1, r2, r3, jd1 + jdf1, jd2 + jdf2, jd3 + jdf3)
            print('using hgibbs: ' % ())
        p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = \
            sc.rv2coeh(r2, v2, re, mu)
        ###print('coes init ans %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f\n' % (p, a, ecc, incl * rad2deg, omega * rad2deg, argp * rad2deg, nu * rad2deg, m * rad2deg))
        #fprintf(1, 'dr #11.7f m #11.7f m/s \n', 1000*smu.mag(r2-r2ans), 1000*smu.mag(v2-v2ans))
        if (ll <= 8):
            # --- now get an improved estimate of the f and g series --
            u = mu / (magr2 * magr2 * magr2)
            rdot = np.dot(r2, v2) / magr2
            udot = (- 3.0 * mu * rdot) / (magr2 ** 4)
            print('u %17.15f rdot  %11.7f udot %11.7f \n' % (u, rdot, udot))
            tausqr = tau12 * tau12
            f1 = 1.0 - 0.5 * u * tausqr - (1.0 / 6.0) * udot * tausqr * tau12
            #                       - (1.0/24.0) * u*u*tausqr*tausqr
#                       - (1.0/30.0)*u*udot*tausqr*tausqr*tau1
            g1 = (tau12 - (1.0 / 6.0) * u * tau12 * tausqr
                  - (1.0 / 12.0) * udot * tausqr * tausqr)
            #                       - (1.0/120.0)*u*u*tausqr*tausqr*tau1
#                       - (1.0/120.0)*u*udot*tausqr*tausqr*tausqr
            tausqr = tau32 * tau32
            f3 = 1.0 - 0.5 * u * tausqr - (1.0 / 6.0) * udot * tausqr * tau32
            #                       - (1.0/24.0) * u*u*tausqr*tausqr
#                       - (1.0/30.0)*u*udot*tausqr*tausqr*tau3
            g3 = (tau32 - (1.0 / 6.0) * u * tau32 * tausqr
                  - (1.0 / 12.0) * udot * tausqr * tausqr)
            #                       - (1.0/120.0)*u*u*tausqr*tausqr*tau3
#                       - (1.0/120.0)*u*udot*tausqr*tausqr*tausqr
            print('f1 %11.7f g1 %11.7f f3 %11.7f g3 %11.7f \n' % (f1, g1, f3, g3))
        else:
            # -------- use exact method to find f and g -----------
            theta = smu.angl(r1, r2)
            theta1 = smu.angl(r2, r3)
            f1 = 1.0 - ((magr1 * (1.0 - np.cos(theta)) / p))
            g1 = (magr1 * magr2 * np.sin(- theta)) / np.sqrt(p)
            f3 = 1.0 - ((magr3 * (1.0 - np.cos(theta1)) / p))
            g3 = (magr3 * magr2 * np.sin(theta1)) / np.sqrt(p)
        c1 = g3 / (f1 * g3 - f3 * g1)
        c3 = - g1 / (f1 * g3 - f3 * g1)
        print(' c1 %11.7f c3 %11.7f %11.7f \n' % (c1, c2, c3))
        # ----- solve for all three ranges via matrix equation ----
        cmat[0, 0] = - c1
        cmat[1, 0] = - c2
        cmat[2, 0] = - c3
        rhomat = lir * cmat
        print(rhomat)
        print('rhomat %11.7f %11.7f %11.7f \n'
              % (rhomat[0, 0], rhomat[0, 1], rhomat[0, 2]))
        #        fprintf(1, 'rhomat #11.7f #11.7f #11.7f \n', rhomat/re)
        rhoold1 = rhomat[0, 0] / c1
        rhoold2 = rhomat[1, 0] / c2
        rhoold3 = rhomat[2, 0] / c3
        print('rhoold %11.7f %11.7f %11.7f \n' % (rhoold1, rhoold2, rhoold3))
        #   fprintf(1, 'rhoold #11.7f #11.7f #11.7f \n', rhoold1/re, rhoold2/re, rhoold3/re)
        for i in range(3):
            r1[i] = rhomat[0, 0] * l1eci[i] / c1 + rs1[i]
            r2[i] = rhomat[1, 0] * l2eci[i] / c2 + rs2[i]
            r3[i] = rhomat[2, 0] * l3eci[i] / c3 + rs3[i]
        print('r1 %11.7f %11.7f %11.7f \n' % (r1[0], r1[1], r1[2]))
        print('r2 %11.7f %11.7f %11.7f \n' % (r2[0], r2[1], r2[2]))
        print('r3 %11.7f %11.7f %11.7f \n' % (r3[0], r3[1], r3[2]))
        print('====================next loop \n' % ())
        # ----------------- check for convergence -----------------
        print('rhoold while  %16.14f %16.14f \n' % (rhoold2, rho2))


    # ---------------- find all three vectors ri ------------------
    for i in range(3):
        r1[i] = rhomat[0, 0] * l1eci[i] / c1 + rs1[i]
        r2[i] = - rhomat[1, 0] * l2eci[i] + rs2[i]
        r3[i] = rhomat[2, 0] * l3eci[i] / c3 + rs3[i]

    return r2, v2

# ------------------------------------------------------------------------------
#
#                           function anglesl
#
#  this function solves the problem of orbit determination using three
#    optical sightings and the method of laplace.
#
#  author        : david vallado                  719-573-2600   24 apr 2003
#
#  inputs          description                    range / units
#    re           - radius earth, sun, etc        km
#    mu           - grav param earth, sun etc     km3/s2
#    rtasc1       - right ascension #1                          rad
#    rtasc2       - right ascension #2                          rad
#    rtasc3       - right ascension #3                          rad
#    decl1        - declination #1                              rad
#    decl2        - declination #2                              rad
#    decl3        - declination #3                              rad
#    jd1, jdf1    - julian date of 1st sighting                 days from 4713 bc
#    jd2, jdf2    - julian date of 2nd sighting                 days from 4713 bc
#    jd3, jdf3    - julian date of 3rd sighting                 days from 4713 bc
#    rsite1       - eci site position vector                    km
#    rsite2       - eci site position vector                    km
#    rsite3       - eci site position vector                    km
#
#  outputs        :
#    r            - ijk position vector           km
#    v            - ijk velocity vector           km / s
#
#  locals         :
#    l1           - line of sight vector for 1st
#    l2           - line of sight vector for 2nd
#    l3           - line of sight vector for 3rd
#    ldot         - 1st derivative of l2
#    lddot        - 2nd derivative of l2
#    rs2dot       - 1st derivative of rs2 - vel
#    rs2ddot      - 2nd derivative of rs2
#    t12t13       - (t1-t2) * (t1-t3)
#    t21t23       - (t2-t1) * (t2-t3)
#    t31t32       - (t3-t1) * (t3-t2)
#    i            - index
#    d            -
#    d1           -
#    d2           -
#    d3           -
#    d4           -
#    oldr         - previous iteration on r
#    rho          - range from site to satellite at t2
#    rhodot       -
#    dmat         -
#    d1mat        -
#    d2mat        -
#    d3mat        -
#    d4mat        -
#    earthrate    - angular rotation of the earth
#    l2dotrs      - vector l2 dotted with rs
#    temp         - temporary vector
#    temp1        - temporary vector
#    small        - tolerance
#    roots        -
#
#  coupling       :
#    mag          - magnitude of a vector
#    determinant  - evaluate the determinant of a matrix
#    cross        - cross product of two vectors
#    unit         - unit vector
#    factor       - find the roots of a polynomial
#
#  references     :
#    vallado       2001, 413-417
#
# [r2, v2] = anglesl (decl1, decl2, decl3, rtasc1, rtasc2, rtasc3, jd1, jdf1, jd2, jdf2, jd3, jdf3, rs1, rs2, rs3)
# ------------------------------------------------------------------------------

def anglesl(decl1=None, decl2=None, decl3=None, rtasc1=None,
            rtasc2=None, rtasc3=None, jd1=None, jdf1=None, jd2=None,
            jdf2=None, jd3=None, jdf3=None, diffsites=None, rs1=None,
            rs2=None, rs3=None):
    # -------------------------  implementation   -------------------------
#     omegaearth = 0.000072921158553  # earth rad/s
#     omegaearth = 0.017202791208627  # sun rad/s
#     omegaearth = 2.0 * pi/365.24221897  # au / day
    earthrate = np.zeros(3)
    earthrate[0] = 0.0
    earthrate[1] = 0.0
    earthrate[2] = earthrot
    #        tuday = 58.132440906
#        mu   = 1.32712428e11
# need to switch these for interplanetary

    #        constant

    small = 1e-08
    # ---------- set middle to 0, find deltas to others -----------

    # test problem///////////////////////////////////////////////////////
#      los1 = [-0.425365 0.777650 0.462953]  # just in km
#      los2 =[-0.825702 0.259424 0.500914]
#      los3 = [ -0.947067 -0.129576 0.293725]

    # los1 =[-0.425364592304 , 0.777650239833, 0.462952554914]
# los2 =[ -0.825702365309, 0.259423566604, 0.500914181287]
# los3 =[-0.947067028031, -0.129575647726, 0.2937246941]

    t1 = - 1200

    t2 = 0
    t3 = 1200
    tau12 = t1 - t2
    tau13 = t1 - t3
    tau32 = t3 - t2
    # test problem///////////////////////////////////////////////////////

    tau12 = (jd1 - jd2) * 86400.0 + (jdf1 - jdf2) * 86400.0

    tau13 = (jd1 - jd3) * 86400.0 + (jdf1 - jdf3) * 86400.0
    tau32 = (jd3 - jd2) * 86400.0 + (jdf3 - jdf2) * 86400.0
    # --------------- find line of sight vectors ------------------
# should be eci...
    los1 = np.zeros(3)
    los1[0] = np.cos(decl1) * np.cos(rtasc1)
    los1[1] = np.cos(decl1) * np.sin(rtasc1)
    los1[2] = np.sin(decl1)
    los2 = np.zeros(3)
    los2[0] = np.cos(decl2) * np.cos(rtasc2)
    los2[1] = np.cos(decl2) * np.sin(rtasc2)
    los2[2] = np.sin(decl2)
    los3 = np.zeros(3)
    los3[0] = np.cos(decl3) * np.cos(rtasc3)
    los3[1] = np.cos(decl3) * np.sin(rtasc3)
    los3[2] = np.sin(decl3)
    # same- they're both unit vectors
#     l1
#     unit(l1)

    los1
    los2
    los3
    # -------------------------------------------------------------
#       using lagrange interpolation formula to derive an expression
#       for l(t), substitute t = t2 and differentiate to obtain the
#       derivatives of l.
# -------------------------------------------------------------
    s1 = - tau32 / (tau12 * tau13)
    s2 = (tau12 + tau32) / (tau12 * tau32)

    s3 = - tau12 / (- tau13 * tau32)
    s4 = 2.0 / (tau12 * tau13)
    s5 = 2.0 / (tau12 * tau32)
    s6 = 2.0 / (- tau13 * tau32)
    ldot = np.zeros(3)
    lddot = np.zeros(3)
    for i in range(3):
        ldot[i] = s1 * los1[i] + s2 * los2[i] + s3 * los3[i]
        lddot[i] = s4 * los1[i] + s5 * los2[i] + s6 * los3[i]

    ldot
    lddot
    #    ldotmag = smu.mag(ldot)
#    lddotmag = smu.mag(lddot)
# should these unit vectors use a diff name????????//
#    ldot = unit(ldot)
#    lddot = unit(lddot)
#    ldot
#    lddot

    # ------------------- find 2nd derivative of rs ---------------
#     temp = cross(rs1, rs2)
#     temp1 = cross(rs2, rs3)

    #     #      needs a different test xxxx##
#     if ((smu.mag(temp) > small) & (smu.mag(temp1) > small))
#         #          fix this test here

    if diffsites == 'n':
        # ------------ all sightings from one site -----------------
        rs2dot = np.cross(earthrate, rs2)
        rs2ddot = np.cross(earthrate, rs2dot)
    else:
        # ---------- each sighting from a different site ----------
        rs2dot = np.zeros(3)
        rs2ddot = np.zeros(3)
        for i in range(3):
            rs2dot[i] = s1 * rs1[i] + s2 * rs2[i] + s3 * rs3[i]
            rs2ddot[i] = s4 * rs1[i] + s5 * rs2[i] + s6 * rs3[i]

    rs2dot
    rs2ddot
    dmat = np.zeros((3, 3))
    dmat1 = np.zeros((3, 3))
    dmat2 = np.zeros((3, 3))
    dmat3 = np.zeros((3, 3))
    dmat4 = np.zeros((3, 3))
    for i in range(3):
        dmat[i, 0] = los2[i]
        dmat[i, 1] = ldot[i]
        dmat[i, 2] = lddot[i]
        dmat
        # ----------------  position determinants -----------------
        dmat1[i, 0] = los2[i]
        dmat1[i, 1] = ldot[i]
        dmat1[i, 2] = rs2ddot[i]
        dmat2[i, 0] = los2[i]
        dmat2[i, 1] = ldot[i]
        dmat2[i, 2] = rs2[i]
        # ------------  velocity determinants ---------------------
        dmat3[i, 0] = los2[i]
        dmat3[i, 1] = rs2ddot[i]
        dmat3[i, 2] = lddot[i]
        dmat4[i, 0] = los2[i]
        dmat4[i, 1] = rs2[i]
        dmat4[i, 2] = lddot[i]

    dmat1
    dmat2
    d = 2.0 * np.linalg.det(dmat)
    d1 = np.linalg.det(dmat1)
    d2 = np.linalg.det(dmat2)
    d3 = np.linalg.det(dmat3)
    d4 = np.linalg.det(dmat4)
    print('d, di  %11.6g %11.6g %11.6g %11.6g %11.6g \n' % (d, d1, d2, d3, d4))

    # ---------------  iterate to find rho magnitude ----------------
#     magr = 1.5   # first guess
#     writeln('input initial guess for magr ')
#     readln(magr)
#     i = 1
#     repeat
#         oldr = magr
#         rho = -2.0*d1/d - 2.0*d2/(magr*magr*magr*d)
#         magr = sqrt(rho*rho + 2.0*rho*l2dotrs + rs2(4)*rs2(4))
#         inc[i]
#         magr = (oldr - magr) / 2.0             # simple bissection
#         writeln(fileout, 'rho guesses ', i:2, 'rho ', rho:14:7, ' magr ', magr:14:7, oldr:14:7)
# seems to converge, but wrong numbers
#         inc[i]
#     until (abs(oldr-magr) < small) | (i .ge. 30)

    if (np.abs(d) > 1.0 - 14):
        poly = np.zeros(9)
        # --------------- solve eighth order poly -----------------
        l2dotrs = np.dot(los2, rs2)
        poly[0] = 1.0
        poly[1] = 0.0
        poly[2] = (l2dotrs * 4.0 * d1 / d - 4.0 * d1 * d1 / (d * d)
                   - smu.mag(rs2) * smu.mag(rs2))
        poly[3] = 0.0
        poly[4] = 0.0
        poly[5] = mu * (l2dotrs * 4.0 * d2 / d - 8.0 * d1 * d2 / (d * d))
        poly[6] = 0.0
        poly[7] = 0.0
        poly[8] = - 4.0 * mu * mu * d2 * d2 / (d * d)
        rootarr = np.roots(poly)
        ###poly(2)
        ###poly(5)
        ###poly(8)
        rootarr
        x = rootarr[0]
        print("rootarr:")
        print(rootarr)

        ###poly(1) * x ** 8 + poly(3) * x ** 6 + poly(6) * x ** 3 + poly(9)
        # ------------------ find correct (xx) root ----------------
        bigr2 = 0.0
        for j in range(8):
            #                if (abs(roots(j, 2)) < small)
#                    writeln('root ', j, roots(j, 1), ' + ', roots(j, 2), 'j')
#        temproot = roots(j, 1)*roots(j, 1)
#        temproot = temproot*temproot*temproot*temproot +
#                  poly(3)*temproot*temproot*temproot + poly(6)*roots(j, 1)*temproot + poly(9)
#                    writeln(fileout, 'root ', j, roots(j, 1), ' + ', roots(j, 2), 'j  value = ', temproot)
            if (rootarr[j] > bigr2):
                bigr2 = rootarr[j]
            #                  end
        print('bigr2 %11.7f  %11.7f er \n' % (bigr2, bigr2 / re))
        #  fprintf(1, 'keep this root ? ')
#        input (bigr2)
        rho = - 2.0 * d1 / d - 2.0 * mu * d2 / (bigr2 * bigr2 * bigr2 * d)
        print('rho %11.7f  %11.7f er \n' % (rho, rho / re))
        r2 = np.zeros(3)
        # --------- find the middle position vector ---------------
        for k in range(3):
            r2[k] = rho * los2[k] + rs2[k]
        magr2 = smu.mag(r2)
        # ---------------- find rhodot magnitude ------------------
        rhodot = - d3 / d - mu * d4 / (magr2 * magr2 * magr2 * d)
        #        writeln(fileout, 'rho ', rho:14:7)
#        writeln(fileout, 'rhodot ', rhodot:14:7)
        # -------------- find middle velocity vector --------------
        v2 = np.zeros(3)
        for i in range(3):
            v2[i] = rhodot * los2[i] + rho * ldot[i] + rs2dot[i]
    else:
        print('determinant value was zero %11.7f ' % (d))

    return r2, v2


#
# ------------------------------------------------------------------------------
#
#                           function gibbs
#
#  this function performs the gibbs method of orbit determination.  this
#    method determines the velocity at the middle point of the 3 given position
#    vectors.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    r1          - ijk position vector #1         km
#    r2          - ijk position vector #2         km
#    r3          - ijk position vector #3         km
#
#  outputs       :
#    v2          - ijk velocity vector for r2     km / s
#    theta       - angl between vectors           rad
#    error       - flag indicating success        'ok', ...
#
#  locals        :
#    tover2      -
#    l           -
#    small       - tolerance for roundoff errors
#    r1mr2       - magnitude of r1 - r2
#    r3mr1       - magnitude of r3 - r1
#    r2mr3       - magnitude of r2 - r3
#    p           - p vector     r2 x r3
#    q           - q vector     r3 x r1
#    w           - w vector     r1 x r2
#    d           - d vector     p + q + w
#    n           - n vector (r1)p + (r2)q + (r3)w
#    s           - s vector
#                    (r2-r3)r1+(r3-r1)r2+(r1-r2)r3
#    b           - b vector     d x r2
#    theta1      - temp angl between the vectors   rad
#    pn          - p unit vector
#    r1n         - r1 unit vector
#    dn          - d unit vector
#    nn          - n unit vector
#    i           - index
#
#  coupling      :
#    mag         - magnitude of a vector
#    cross       - cross product of two vectors
#    dot         - dot product of two vectors
#    unit        - unit vector
#    angl       - angl between two vectors
#
#  references    :
#    vallado       2007, 456, alg 52, ex 7-5
#
# [v2, theta, theta1, copa, error] = function gibbs(r1, r2, r3)
# ------------------------------------------------------------------------------

def gibbs(r1=None, r2=None, r3=None):
    # -------------------------  implementation   -------------------------
    small = 1e-06
    theta = 0.0
    error = '          ok'
    theta1 = 0.0
    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)
    magr3 = smu.mag(r3)
    v2 = np.zeros(3)

    p = np.cross(r2, r3)
    q = np.cross(r3, r1)
    w = np.cross(r1, r2)
    pn = smu.unit(p)
    r1n = smu.unit(r1)
    copa = np.arcsin(np.dot(pn, r1n))
    if (np.abs(copa) > 0.017452406):
        error = 'not coplanar'

    # --------------- | don't continue processing --------------
    d = p + q + w
    magd = smu.mag(d)
    n = magr1 * p + magr2 * q + magr3 * w
    magn = smu.mag(n)
    nn = smu.unit(n)
    dn = smu.unit(d)
    # -------------------------------------------------------------
#       determine if  the orbit is possible.  both d and n must be in
#         the same direction, and non-zero.
# -------------------------------------------------------------
    if (((np.abs(magd) < small) or (np.abs(magn) < small)) or
            (np.dot(nn, dn) < small)):
        error = '  impossible'
    else:
        theta = smu.angl(r1, r2)
        theta1 = smu.angl(r2, r3)
        # ----------- perform gibbs method to find v2 -----------
        r1mr2 = magr1 - magr2
        r3mr1 = magr3 - magr1
        r2mr3 = magr2 - magr3
        s = r1mr2 * r3 + r3mr1 * r2 + r2mr3 * r1
        b = np.cross(d, r2)
        l = np.sqrt(mu / (magd * magn))
        tover2 = l / magr2
        v2 = tover2 * b + l * s

    #    fprintf(1, 'p     #11.7f   #11.7f  #11.7f km2 \n', p)
#    fprintf(1, 'n     #11.7f   #11.7f  #11.7f km3 \n', n)
#    fprintf(1, 'd     #11.7f   #11.7f  #11.7f km3 \n', d)
#    fprintf(1, 's     #11.7f   #11.7f  #11.7f km2 \n', s)
#    fprintf(1, 'theta     #11.7f   #11.7f deg \n', theta*180/pi, theta1*180/pi)
#    fprintf(1, 'b     #11.7f   #11.7f  #11.7f km3 \n', b)
#    fprintf(1, 'l     #11.7f  /kms \n', l)
    return v2, theta, theta1, copa, error


# ------------------------------------------------------------------------------
#
#                           function gibbsh
#
#  this function performs the gibbs method of orbit determination.  this
#    method determines the velocity at the middle point of the 3 given position
#    vectors.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    r1          - ijk position vector #1         km
#    r2          - ijk position vector #2         km
#    r3          - ijk position vector #3         km
#
#  outputs       :
#    v2          - ijk velocity vector for r2     km / s
#    theta       - angl between vectors           rad
#    error       - flag indicating success        'ok', ...
#
#  locals        :
#    tover2      -
#    l           -
#    small       - tolerance for roundoff errors
#    r1mr2       - magnitude of r1 - r2
#    r3mr1       - magnitude of r3 - r1
#    r2mr3       - magnitude of r2 - r3
#    p           - p vector     r2 x r3
#    q           - q vector     r3 x r1
#    w           - w vector     r1 x r2
#    d           - d vector     p + q + w
#    n           - n vector (r1)p + (r2)q + (r3)w
#    s           - s vector
#                    (r2-r3)r1+(r3-r1)r2+(r1-r2)r3
#    b           - b vector     d x r2
#    theta1      - temp angl between the vectors   rad
#    pn          - p unit vector
#    r1n         - r1 unit vector
#    dn          - d unit vector
#    nn          - n unit vector
#    i           - index
#
#  coupling      :
#    mag         - magnitude of a vector
#    cross       - cross product of two vectors
#    dot         - dot product of two vectors
#    unit        - unit vector
#    angl       - angl between two vectors
#
#  references    :
#    vallado       2007, 456, alg 52, ex 7-5
#
# [v2, theta, theta1, copa, error] = function gibbs(r1, r2, r3)
# ------------------------------------------------------------------------------

def gibbsh(r1=None, r2=None, r3=None, re=None, mu=None):
    # -------------------------  implementation   -------------------------
    small = 1e-06
    theta = 0.0
    error = '          ok'
    theta1 = 0.0
    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)
    magr3 = smu.mag(r3)
    v2 = np.zeros(3)
    for i in range(1):
        v2[i] = 0.0

    p = np.cross(r2, r3)
    q = np.cross(r3, r1)
    w = np.cross(r1, r2)
    pn = np.linalg.norm(p)
    r1n = np.linalg.norm(r1)
    copa = np.arcsin(np.dot(pn, r1n))
    if (np.abs(np.dot(r1n, pn)) > 0.017452406):
        error = 'not coplanar'

    # --------------- | don't continue processing --------------
    d = p + q + w
    magd = smu.mag(d)
    n = magr1 * p + magr2 * q + magr3 * w
    magn = smu.mag(n)
    nn = np.linalg.norm(n)
    dn = np.linalg.norm(d)
    # -------------------------------------------------------------
#       determine if  the orbit is possible.  both d and n must be in
#         the same direction, and non-zero.
# -------------------------------------------------------------
    if ((((np.abs(magd) < small) or (np.abs(magn) < small)) or
            (np.dot(nn, dn) < small))):
        error = '  impossible'
    else:
        theta = smu.angl(r1, r2)
        theta1 = smu.angl(r2, r3)
        # ----------- perform gibbs method to find v2 -----------
        r1mr2 = magr1 - magr2
        r3mr1 = magr3 - magr1
        r2mr3 = magr2 - magr3
        s = r1mr2 * r3 + r3mr1 * r2 + r2mr3 * r1
        b = np.cross(d, r2)
        l = np.sqrt(mu / (magd * magn))
        tover2 = l / magr2
        v2 = tover2 * b + l * s

    return v2, theta, theta1, copa, error

#
# ------------------------------------------------------------------------------
#
#                           function hgibbs
#
#  this function implements the herrick-gibbs approximation for orbit
#    determination, and finds the middle velocity vector for the 3 given
#    position vectors.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    r1          - ijk position vector #1         km
#    r2          - ijk position vector #2         km
#    r3          - ijk position vector #3         km
#    use seconds to provide more accuracy
#    t1         - time julian date of 1st sighting    days from 4713 bc sec
#    t2         - time julian date of 2nd sighting    days from 4713 bc sec
#    t3         - time julian date of 3rd sighting    days from 4713 bc sec
#
#  outputs       :
#    v2          - ijk velocity vector for r2     km / s
#    theta       - angl between vectors          rad
#    error       - flag indicating success        'ok', ...
#
#  locals        :
#    dt21        - time delta between r1 and r2   s
#    dt31        - time delta between r3 and r1   s
#    dt32        - time delta between r3 and r2   s
#    p           - p vector    r2 x r3
#    pn          - p unit vector
#    r1n         - r1 unit vector
#    theta1      - temporary angl between vec    rad
#    tolangle    - tolerance angl  (1 deg)       rad
#    term1       - 1st term for hgibbs expansion
#    term2       - 2nd term for hgibbs expansion
#    term3       - 3rd term for hgibbs expansion
#    i           - index
#
#  coupling      :
#    mag         - magnitude of a vector
#    cross       - cross product of two vectors
#    dot         - dot product of two vectors
#    unit        - unit vector
#    lncom3      - combination of three scalars and three vectors
#    angl       - angl between two vectors
#
#  references    :
#    vallado       2007, 462, alg 52, ex 7-4
#
# [v2, theta, theta1, copa, error ] = hgibbs (r1, r2, r3, jd1, jd2, jd3)
# ------------------------------------------------------------------------------

def hgibbs(r1=None, r2=None, r3=None, jd1=None, jd2=None, jd3=None):
    # -------------------------  implementation   -------------------------
    error = '          ok'
    theta = 0.0
    theta1 = 0.0
    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)
    magr3 = smu.mag(r3)
    v2 = np.zeros(3)

    tolangle = 0.01745329251994
    dt21 = (jd2 - jd1) * 86400.0
    dt31 = (jd3 - jd1) * 86400.0

    dt32 = (jd3 - jd2) * 86400.0
    p = np.cross(r2, r3)
    pn = smu.unit(p)
    r1n = smu.unit(r1)
    copa = np.arcsin(np.dot(pn, r1n))
    if (np.abs(copa) > 0.017452406):
        error = 'not coplanar'

    # --------------------------------------------------------------
#       check the size of the angles between the three position vectors.
#       herrick gibbs only gives "reasonable" answers when the
#       position vectors are reasonably close.  10 deg is only an estimate.
# --------------------------------------------------------------
    theta = smu.angl(r1, r2)
    theta1 = smu.angl(r2, r3)
    if (((theta > tolangle) or (theta1 > tolangle))):
        error = '   angl > 1'

    # ----------- perform herrick-gibbs method to find v2 ---------
    term1 = - dt32 * (1.0 / (dt21 * dt31) + mu / (12.0 * magr1 * magr1 * magr1))
    term2 = (dt32 - dt21) * (1.0 / (dt21 * dt32)
                             + mu / (12.0 * magr2 * magr2 * magr2))
    term3 = dt21 * (1.0 / (dt32 * dt31) + mu / (12.0 * magr3 * magr3 * magr3))
    v2 = term1 * r1 + term2 * r2 + term3 * r3
    #     fprintf(1, 'p     #11.7f   #11.7f  #11.7f km2 \n', p)
#     fprintf(1, 'theta     #11.7f   #11.7f deg \n', theta*180/pi, theta1*180/pi)
    return v2, theta, theta1, copa, error



# ------------------------------------------------------------------------------
#
#                           function hillsv
#
#  this function calculates initial velocity for hills equations.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    r           - initial position vector of int m
#    alt         - altitude of tgt satellite      km
#    dts         - desired time                   s
#
#  outputs       :
#    v           - initial velocity vector of int m / s
#
#  locals        :
#    numkm       -
#    denom       -
#    nt          - angular velocity times time    rad
#    omega       -
#    sinnt       - sine of nt
#    cosnt       - cosine of nt
#    radius      - magnitude of range vector      km
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 410, eq 6-60, ex 6-15
#
# [v] = hillsv(r, alt, dts)
# ------------------------------------------------------------------------------


def hillsv(r, alt, dts):

    # --------------------  implementation   ----------------------

    radius = re + alt
    omega = math.sqrt(mu / (radius*radius*radius))
    nt = omega*dts
    cosnt = math.cos(nt)
    sinnt = math.sin(nt)

    # --------------- determine initial velocity ------------------
    numkm = ((6.0*r[0]*(nt-sinnt)-r[1])*omega*sinnt \
            - 2.0*omega*r[0]*(4.0-3.0*cosnt)*(1.0-cosnt))
    denom = (4.0*sinnt-3.0*nt)*sinnt + 4.0*(1.0-cosnt) \
            *(1.0-cosnt)

    v = np.zeros(3)
    if (abs(denom) > 0.000001):
        v[1] = numkm / denom
    else:
        v[1] = 0.0
    if (abs(sinnt) > 0.000001):
        v[0] = -(omega*r[0]*(4.0-3.0*cosnt) \
              +2.0*(1.0-cosnt)*v[1]) / (sinnt)
    else:
        v[0] = 0.0
    v[2] = -r[2]*omega*smu.cot(nt)

    return v



# ------------------------------------------------------------------------------
#
#                           function hillsr
#
#  this function calculates various position information for hills equations.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    r           - init rel position of int       m or km
#    v           - init rel velocity of int       m or km/s
#    alt         - altitude of tgt satellite      km
#    dts         - desired time                   s
#
#  outputs       :
#    rinit       - final rel position of int      m or km
#    vinit       - final rel velocity of int      m or km/s
#
#  locals        :
#    nt          - angular velocity times time    rad
#    omega       -
#    sinnt       - sine of nt
#    cosnt       - cosine of nt
#    radius      - magnitude of range vector      km
#
#  coupling      :
#
#
#  references    :
#    vallado       2007, 397, alg 47, ex 6-14
#
#  [rint, vint] = hillsr(r, v, alt, dts)
# ------------------------------------------------------------------------------


def hillsr(r, v, alt, dts):

    # --------------------  implementation   ----------------------

    radius = re + alt # in km
    omega = math.sqrt(mu / (radius*radius*radius)) # rad/s
    nt = omega * dts
    cosnt = math.cos(nt)
    sinnt = math.sin(nt)

    # --------------- determine new positions  --------------------
    rint = np.array([0.0, 0.0, 0.0])
    rint[0] = (v[0]/omega) * sinnt - \
             ((2.0*v[1]/omega) + 3.0*r[0]) * cosnt + \
             ((2.0*v[1]/omega) + 4.0*r[0])
    rint[1] = (2.0*v[0]/omega) * cosnt + \
             ((4.0*v[1]/omega) + 6.0*r[0]) * sinnt + \
             (r[1] - (2.0*v[0]/omega)) - \
             (3.0*v[1] + 6.0*omega*r[0])*dts
    rint[2] = r[2]*cosnt + (v[2]/omega)*sinnt

    # --------------- determine new velocities  -------------------
    vint = np.array([0.0, 0.0, 0.0])
    vint[0] = v[0]*cosnt + (2.0*v[1]+3.0*omega*r[0])*sinnt
    vint[1] = -2.0*v[0]*sinnt + (4.0*v[1] \
              +6.0*omega*r[0])*cosnt - (3.0*v[1]+6.0*omega*r[0])
    vint[2] = -r[2]*omega*sinnt + v[2]*cosnt

    return rint, vint


# -----
#   vallado       2007, 370, alg 46, ex 6-10
# function [ttrans, tphase, dvphase, dvtrans1, dvtrans2, aphase ] = noncoplr(phasei, aint, atgt, ktgt, kint, arglatint, nodeint, truelon, deltai)
#------ }

def noncoplr(phasei=None, aint=None, atgt=None, ktgt=None, kint=None,
             arglatint=None, nodeint=None, truelon=None, deltai=None):
    mu = 1.0

    angvelint = np.sqrt(mu / (aint * aint * aint))
    angveltgt = np.sqrt(mu / (atgt * atgt * atgt))
    atrans = (aint + atgt) * 0.5
    ttrans = np.pi * np.sqrt((atrans * atrans * atrans) / mu)
    deltatnode = phasei / angvelint
    lead = angveltgt * ttrans
    omeganode = angveltgt * deltatnode
    phasenew = nodeint + np.pi - (truelon + omeganode)
    leadnew = np.pi + phasenew
    tphase = (leadnew - lead + twopi * ktgt) / angveltgt
    aphase = (mu * (tphase / (twopi * kint)) ** 2) ** (1.0 / 3.0)
    # -----------------  find deltav's  ----------------- }
    vint = np.sqrt(mu / aint)
    vphase = np.sqrt(2.0 * mu / aint - mu / aphase)
    dvphase = vphase - vint
    vtrans1 = np.sqrt(2.0 * mu / aint - mu / atrans)
    dvtrans1 = vtrans1 - vphase
    vtrans2 = np.sqrt(2.0 * mu / atgt - mu / atrans)
    vtgt = np.sqrt(mu / atgt)
    dvtrans2 = np.sqrt(vtgt * vtgt + vtrans2 * vtrans2
                       - 2.0 * vtgt * vtrans2 * np.cos(deltai))
    dvtotal = dvphase + dvtrans1 + dvtrans2
    ttotal = deltatnode + ttrans + tphase
    print(' angvelint %11.7f %11.7f rad/s  \n' % (angvelint, angvelint / tusec))
    print(' angveltgt %11.7f %11.7f rad/s  \n' % (angveltgt, angveltgt / tusec))
    print(' atrans %11.7f %11.7f km \n' % (atrans, atrans * 6378.137))
    print(' ttrans %11.7f %11.7f min \n' % (ttrans, ttrans * tumin))
    print(' deltatnode  %11.7f %11.7f min \n' % (deltatnode, deltatnode * tumin))
    print(' lead  %11.7f \n' % (lead * rad2deg))
    print(' omeganode  %11.7f \n' % (omeganode * rad2deg))
    print(' phasenew  %11.7f \n' % (phasenew * rad2deg))
    print(' leadnew  %11.7f \n' % (leadnew * rad2deg))
    print(' tphase %11.7f %11.7f min \n' % (tphase, tphase * tumin))
    print(' aphase  %11.7f %11.7f km \n' % (aphase, aphase * re))
    print(' vint   %11.7f %11.7f km/s \n' % (vint, vint * velkmps))
    print(' vphase  %11.7f %11.7f km/s \n' % (vphase, vphase * velkmps))
    print(' dvphase  %11.7f  %11.7f \n' % (dvphase, dvphase * velkmps))
    print(' vtrans1   %11.7f %11.7f km/s \n' % (vtrans1, vtrans1 * velkmps))
    print(' dvtrans1   %11.7f %11.7f km/s \n' % (dvtrans1, dvtrans1 * velkmps))
    print(' vtrans2   %11.7f  %11.7f km/s \n' % (vtrans2, vtrans2 * velkmps))
    print(' vtgt  %11.7f %11.7f km/s \n' % (vtgt, vtgt * velkmps))
    print(' dvtrans2   %11.7f %11.7f km/s \n' % (dvtrans2, dvtrans2 * velkmps))
    print(' dvtotal   %11.7f %11.7f km/s \n' % (dvtotal, dvtotal * velkmps))
    print(' ttotal   %11.7f %11.7f min \n' % (ttotal, ttotal * tumin))
    return ttrans, tphase, dvphase, dvtrans1, dvtrans2, aphase

# ------------------------------------------------------------------------------
#
#                           procedure rendz
#
#  this procedure calculates parameters for a hohmann transfer rendezvous.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    rcs1        - radius of circular orbit int   er
#    rcs2        - radius of circular orbit tgt   er
#    einit       - ecc of first orbit
#    efinal      - ecc of final orbit
#    nuinit      - true anomaly of first orbit    0 or pi rad
#    nufinal     - true anomaly of final orbit    0 or pi rad
#    phasei      - initial phase angle (tgt-int)  +(ahead) or -(behind) rad
#    numrevs     - number of revs to wait
#    ktgt        -
#    kint        -
#
#  outputs       :
#    phasef      - final phase angle              rad
#    waittime    - wait before next intercept opp tu
#    deltav      - change in velocity             er/tu
#
#  locals        :
#    dttutrans   - time of flight of trans orbit  tu
#    atrans      - semimajor axis of trans orbit  er
#    angveltgt   - angular velocity of target     rad / tu
#    angvelint   - angular velocity of int        rad / tu
#    leadang     - lead angle                     rad
#
#  coupling      :
#    power       - raise a base to a power
#
#  references    :
#    vallado       2007, 364, alg 44, alg 45, ex 6-8, ex 6-9
#function [ phasef, waittime, deltav] = rendz(rcs1, rcs3, phasei, einit, efinal, nuinit, nufinal, ktgt, kint)
# ----------------------------------------------------------------------------- }

def rendz(rcs1=None, rcs3=None, phasei=None, einit=None, efinal=None,
          nuinit=None, nufinal=None, ktgt=None, kint=None):

    mu = 1.0
    angvelint = np.sqrt(mu / (rcs1 * rcs1 * rcs1))
    angveltgt = np.sqrt(mu / (rcs3 * rcs3 * rcs3))
    vint = np.sqrt(mu / rcs1)
    print(' angvelint %11.7f %11.7f rad/s  \n' % (angvelint, angvelint / tusec))
    print(' angveltgt %11.7f %11.7f rad/s  \n' % (angveltgt, angveltgt / tusec))
    # ---------- check for satellites in the same orbits ----------- }
    if np.abs(angvelint - angveltgt) < 1e-06:
        periodtrans = (ktgt * twopi + phasei) / angveltgt
        atrans = (periodtrans / (twopi * kint)) ** (2.0 / 3.0)
        rp = 2.0 * atrans - rcs1
        if rp < 1.0:
            print(' error - the transfer orbit intersects the earth ' % ())
        vtrans = np.sqrt(((2.0 * mu) / rcs1) - (mu / atrans)) - np.sqrt(1/rcs1)
        deltav = 2.0 * vtrans
        phasef = phasei
        waittime = periodtrans
        leadang = 0.0
        print('periodtrans %11.7f %11.7f s \n' % (periodtrans, periodtrans * tusec / 60))
        print('vint %11.7f %11.7f km/s \n' % (vint, vint * velkmps))
        print('vtrans %11.7f %11.7f km/s \n' % (vtrans, vtrans * velkmps))
        print('deltavinit %11.7f %11.7f km/s \n' % ((vtrans - vint), (vtrans - vint) * velkmps))
        print('atrans %11.7f %11.7f km \n' % (atrans, atrans * 6378.137))
    else:
        # ---- different orbits
        atrans = (rcs1 + rcs3) / 2.0
        dttutrans = np.pi * np.sqrt(atrans * atrans * atrans / mu)
        leadang = angveltgt * dttutrans
        phasef = np.pi - leadang
        phasef * rad2deg
        if phasef < 0.0:
            phasef = phasef + np.pi
        waittime = (phasef - phasei + 2.0 * np.pi * ktgt) / (angvelint - angveltgt)
        a1 = (rcs1 * (1.0 + einit * np.cos(nuinit))) / (1.0 - einit * einit)
        a2 = (rcs1 + rcs3) / 2.0
        a3 = (rcs3 * (1.0 + efinal * np.cos(nufinal))) / (1.0 - efinal * efinal)
        sme1 = -mu / (2.0 * a1)
        sme2 = -mu / (2.0 * a2)
        sme3 = -mu / (2.0 * a3)
        # -----------------  find delta v at point a  ------------------ }
        vinit = np.sqrt(2.0 * ((mu / rcs1) + sme1))
        vtransa = np.sqrt(2.0 * ((mu / rcs1) + sme2))
        deltava = np.abs(vtransa - vinit)
        # -----------------  find delta v at point b  ------------------ }
        vfinal = np.sqrt(2.0 * ((mu / rcs3) + sme3))
        vtransb = np.sqrt(2.0 * ((mu / rcs3) + sme2))
        deltavb = np.abs(vfinal - vtransb)
        deltav = deltava + deltavb
        print('atrans %11.7f %11.7f km \n' % (atrans, atrans * 6378.137))
        ttrans = np.pi * np.sqrt(atrans ** 3 / 1.0)
        print('ttrans %11.7f %11.7f km \n' % (ttrans, ttrans * tusec / 60.0))
        print('leadang %11.7f %11.7f  \n' % (leadang, leadang * rad2deg))
        print('vinit %11.7f %11.7f km/s \n' % (vinit, vinit * velkmps))
        print('vtransa %11.7f %11.7f km/s \n' % (vtransa, vtransa * velkmps))
        print('vfinal %11.7f %11.7f km/s \n' % (vfinal, vfinal * velkmps))
        print('vtransb %11.7f %11.7f km/s \n' % (vtransb, vtransb * velkmps))
        print('deltava %11.7f %11.7f km/s \n' % (deltava, deltava * velkmps))
        print('deltavb %11.7f %11.7f km/s \n' % (deltavb, deltavb * velkmps))

    return phasef, waittime, deltav



# ------------------------------------------------------------------------------
#
#                           procedure combined
#
#  this procedure calculates the delta v's for a hohmann transfer for either
#    circle to circle, or ellipse to circle.
#
#  author        : david vallado                  719-573-2600   5 may  2012
#
#  inputs          description                    range / units
#    rinit       - initial position magnitude     er
#    rfinal      - final position magnitude       er
#    einit       - eccentricity of first orbit
#    efinal      - eccentricity of final orbit
#    nuinit      - true anomaly of first orbit    0 or pi rad
#    nufinal     - true anomaly of final orbit    0 or pi rad, opp of nuinit
#
#  outputs       :
#    deltava     - change in velocity at point a  er / tu
#    deltavb     - change in velocity at point b  er / tu
#    dttu        - time of flight for the trans   tu
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
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 352-359, ex 6-7
# [deltai1, deltava, deltavb, gam1, gam2] = combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai)
#
# ------------------------------------------------------------------------------

def combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai):
    # --------------------  initialize values   ------------------- }
    mu = 1.0 # canonical

    if sh.show == 'y':
        print('rinit %11.7f %11.7f  rfinal %11.7f  %11.7f \n'
              % (rinit, rinit*6378.137 , rfinal, rfinal*6378.137))

    ainit = (rinit * (1.0 + einit * math.cos(nuinit))) / (1.0 - einit * einit)
    atran = (rinit + rfinal) * 0.5
    afinal = (rfinal * (1.0 + efinal * math.cos(nufinal))) / (1.0 - efinal * efinal)
    sme1 = -mu / (2.0*ainit)
    sme2 = -mu / (2.0*atran)

    if sh.show == 'y':
        print('ainit %11.7f %11.7f \n' % (ainit, ainit*6378.1363))
        print('atran %11.7f %11.7f \n' % (atran, atran*6378.1363))
        print('afinal %11.7f %11.7f \n'% (afinal, afinal*6378.1363))

    # -----------------  find delta v at point a  ----------------- }
    vinit = math.sqrt(2.0*(mu/rinit + sme1))
    vtransa = math.sqrt(2.0*(mu/rinit + sme2))
    #     fpa2a = math.atan((e2*math.sin(nu2a)) / (1.0 + e2*math.cos(nu2a)))
    #     fpa1 = math.atan((einit*math.sin(nuinit)) / (1.0 + einit*math.cos(nuinit)))
    #     deltava = math.sqrt(vtransa*vtransa + vinit*vinit - 2.0*vtransa*vinit* ...
    #                     (math.sin(fpa2a)*math.sin(fpa1)+math.cos(fpa2a)*math.cos(fpa1)*math.cos(deltai)))

    # -----------------  find delta v at point b  ----------------- }
    vfinal = math.sqrt(mu/rfinal)  # assumes circular
    vtransb = math.sqrt(2.0*(mu/rfinal + sme2))
    #     fpa2b = math.atan((e2*math.sin(nu2b)) / (1.0 + e2*math.cos(nu2b)))
    #     fpa3 = math.atan((efinal*math.sin(nufinal)) / (1.0 + efinal*math.cos(nufinal)))

    if sh.show == 'y':
        vkmps = 7.905366149846074
        print('vinit %11.7f %11.7f  vfinal %11.7f  %11.7f \n'
              % (vinit, vinit*vkmps, vfinal, vfinal*vkmps))
        print('vtransa %11.7f %11.7f  vtransb %11.7f  %11.7f \n'
              % (vtransa, vtransa*vkmps, vtransb, vtransb*vkmps))

    #     deltavb = math.sqrt(vtransb*vtransb + vfinal*vfinal - 2.0*vtransb*vfinal* ...
    #(math.sin(fpa2b)*math.sin(fpa3)+math.cos(fpa2b)*math.cos(fpa3)*math.cos(deltai)))

    # -------------- find proportions of inclination change ---------------
    # ----------------- this is the approximate approach ------------------
    ratio = rfinal/rinit
    s = 1.0/deltai * math.atan(math.sin(deltai)/(ratio**1.5 + math.cos(deltai)))
    if sh.show == 'y':
        print(' s %11.7f \n' % s)
    deltai1 = s*deltai
    deltai2 = (1.0-s)*deltai

    deltava = math.sqrt(vinit**2  + vtransa**2 - 2.0*vinit*vtransa*math.cos(deltai1))
    deltavb = math.sqrt(vfinal**2 + vtransb**2 - 2.0*vfinal*vtransb*math.cos(deltai2))

    dttu = math.pi * math.sqrt((atran * atran * atran)/mu)

    # ----------------- figure orientation of the firings -----------------
    gam1 = math.acos(-(vinit**2+deltava**2-vtransa**2) / (2.0*vinit*deltava))
    gam2 = math.acos(-(vtransb**2+deltavb**2-vfinal**2) / (2.0*vtransb*deltavb))

    return deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2



# ------------------------------------------------------------------------------
#
#                           procedure iandnode
#
#  this procedure calculates the delta v's for a change in inclination and
#    right ascension of the ascending node.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    vinit       - initial velocity vector        er/tu
#    iinit       - initial inclination            rad
#    fpa         - flight path angle              rad
#    deltaraan   - change in node                 rad
#    deltai      - change in inclination          rad
#    rfinal      - final position magnitude       er
#
#  outputs       :
#    ifinal      - final inclination              rad
#    deltav      - change in velocity             er/tu
#
#  locals        :
#    arglat      - argument of latitude           rad
#    arglat1     - final argument of latitude     rad
#    theta       -
#
#  coupling      :
#    acos      - arc cosine function
#
#  references    :
#    vallado       2007, 350, alg 41, ex 6-6
#function [deltav] = iandnode(iinit, deltaraan, ifinal, vinit, fpa)
# ----------------------------------------------------------------------------- }


def iandnode(iinit=None, deltaraan=None, ifinal=None, vinit=None, fpa=None):
    deltai = iinit - ifinal
    # variables for speed
    cosdraan = np.cos(deltaraan)
    cosii = np.cos(iinit)
    sinii = np.sin(iinit)
    cosif = np.cos(ifinal)
    sinif = np.sin(ifinal)
    cost = cosii * cosif + sinii * sinif * cosdraan
    theta = np.arccos(cost)
    sint = np.sin(theta)
    deltav = 2.0 * vinit * np.cos(fpa) * np.sin(0.5 * theta)
    arglat = np.arccos((sinif * np.cos(deltaraan) - cost * sinii) / (sint * cosii))
    arglat1 = np.arccos((cosii * sinif - sinii * cosif * cosdraan) / sint)
    print(' theta   %11.7f  \n' % (theta * rad2deg))
    print(' arglat   %11.7f  %11.7f  \n' % (arglat * rad2deg, arglat1 * rad2deg))
    return deltav

# ------------------------------------------------------------------------------
#
#                           procedure nodeonly
#
#  this procedure calculates the delta v's for a change in longitude of
#    ascending node only.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    deltaomega  - change in node                 rad
#    ecc         - ecc of first orbit
#    vinit       - initial velocity vector        er/tu
#    fpa         - flight path angle              rad
#    incl        - inclination                    rad
#
#
#  outputs       :
#    ifinal      - final inclination              rad
#    deltav      - change in velocity             er/tu
#
#  locals        :
#    vfinal      - final velocity vector          er/tu
#    arglat      - argument of latitude           rad
#    arglat1     - final argument of latitude     rad
#    nuinit      - initial true anomaly           rad
#    theta       -
#
#  coupling      :
#    asin      - arc sine function
#    acos      - arc cosine function
#
#  references    :
#    vallado       2007, 349, alg 40, ex 6-5
# function [ifinal, deltav ] = nodeonly(iinit, ecc, deltaomega, vinit, fpa, incl)
# ----------------------------------------------------------------------------- }

def nodeonly(iinit=None, ecc=None, deltaomega=None, vinit=None, fpa=None, incl=None):
    if ecc > 1e-07:
        # ------------------------- elliptical --------------------- }
        theta = np.arctan(np.sin(iinit) * np.tan(deltaomega))
        ifinal = np.arcsin(np.sin(theta) / np.sin(deltaomega))
        deltav = 2.0 * vinit * np.cos(fpa) * np.sin(0.5 * theta)
        arglat = np.pi * 0.5
        arglat1 = np.arccos(np.cos(incl) * np.sin(incl)
                            * (1.0 - np.cos(deltaomega)) / np.sin(theta))
    else:
        # -------------------------- circular ---------------------- }
        ifinal = incl
        theta = np.arccos(np.cos(iinit) * np.cos(iinit)
                          + np.sin(iinit) * np.sin(iinit) * np.cos(deltaomega))
        deltav = 2.0 * vinit * np.sin(0.5 * theta)
        arglat = np.arccos(np.tan(iinit)
                           * (np.cos(deltaomega) - np.cos(theta))
                           / np.sin(theta))
        arglat1 = np.arccos(np.cos(incl) * np.sin(incl)
                            * (1.0 - np.cos(deltaomega)) / np.sin(theta))

    print(' theta   %11.7f deg \n' % (theta * rad2deg))
    print(' arglat   %11.7f  %11.7f  \n' % (arglat * rad2deg, arglat1 * rad2deg))
    return ifinal, deltav

# ------------------------------------------------------------------------------
#
#                           procedure ionlychg
#
#  this procedure calculates the delta v's for a change in inclination only.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    deltai      - change in inclination          rad
#    vinit       - initial velocity vector        er/tu
#    fpa         - flight path angle              rad
#
#  outputs       :
#    deltavionly - answer
#
#  locals        :
#    none.
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 346, alg 39, ex 6-4
#function [ deltavionly] = ionlychg(deltai, vinit, fpa)
# ----------------------------------------------------------------------------- }

def ionlychg(deltai=None, vinit=None, fpa=None):
    deltavionly = 2.0 * vinit * np.cos(fpa) * np.sin(0.5 * deltai)
    return deltavionly


# ------------------------------------------------------------------------------
#
#                           procedure onetang
#
#  this procedure calculates the delta v's for a one tangent transfer for either
#    circle to circle, or ellipse to ellipse.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    rinit       - initial position magnitude     er
#    rfinal      - final position magnitude       er
#    einit       - eccentricity of first orbit
#    efinal      - eccentricity of final orbit
#    nuinit      - true anomaly of first orbit    0 or pi rad
#    nu2         - true anomaly of second orbit   same quad as nuinit, rad
#    nufinal     - true anomaly of final orbit    0 or pi rad
#
#  outputs       :
#    deltava     - change in velocity at point a  er / tu
#    deltavb     - change in velocity at point b  er / tu
#    dttu        - time of flight for the transf  tu
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
#    e           - ecc anomaly of trans at b      rad
#    ratio       - ratio of initial to final
#                    orbit radii
#
#  coupling      :
#    atan2       - arc tangent rountine that solves quadrant ambiguities
#
#  references    :
#    vallado       2007, 335, alg 38, ex 6-3
#function [deltava, deltavb, dttu, etran, atran, vtrana, vtranb ] = onetang(rinit, rfinal, einit, efinal, nuinit, nutran)
# ----------------------------------------------------------------------------- }

def onetang(rinit=None, rfinal=None, einit=None, efinal=None,
            nuinit=None, nutran=None):
    # --------------------  initialize values   ------------------- }
    mu = 1.0
    e = 0.0
    deltava = 0.0
    deltavb = 0.0
    dttu = 0.0
    ratio = rinit / rfinal
    if np.abs(nuinit) < 0.01:
        etran = ((ratio - 1.0) / (np.cos(nutran) - ratio))
        etran = ((- rfinal + rinit) / (rfinal * np.cos(nutran) - rinit))
        eainit = 0.0
    else:
        etran = ((ratio - 1.0) / (np.cos(nutran) + ratio))
        etran = ((- rfinal + rinit) / (rfinal * np.cos(nutran) + rinit))
        eainit = np.pi

    if etran >= 0.0:
        ainit = (rinit * (1.0 + einit * np.cos(nuinit))) / (1.0 - einit * einit)
        afinal = (rfinal * (1.0 + efinal * np.cos(nutran))) / (1.0 - efinal * efinal)
        # nutran is used since it = nufinal!! }
#    fprintf(1, ' ainti and final   #11.7f  #11.7f km \n', ainit*re, afinal*re)
#ainit = rinit
#afinal = rfinal
        if np.abs(etran - 1.0) > 1e-06:
            if np.abs(nuinit) < 0.01:
                atran = ((rinit * (1.0 + etran * np.cos(nuinit)))
                         / (1.0 - etran * etran))
            else:
                atran = ((rinit * (1.0 + etran * np.cos(nuinit)))
                         / (1.0 + etran * etran))
                atran = rinit / (1.0 + etran)
        else:
            atran = infinite
        ptran = rinit * (1.0 + etran)
        atran = ptran / (1.0 - etran ** 2)
        # -----------------  find delta v at point a  ----------------- }
        vinit = np.sqrt(mu / rinit)
        vtrana = np.sqrt((2.0 * mu) / rinit - (mu / atran))
        deltava = np.abs(vtrana - vinit)
        # -----------------  find delta v at point b  ----------------- }
        vfinal = np.sqrt((2.0 * mu) / rfinal - (mu / afinal))
        vtranb = np.sqrt((2.0 * mu) / rfinal - (mu / atran))
        fpatranb = np.arctan((etran * np.sin(nutran))
                             / (1.0 + etran * np.cos(nutran)))
        fpafinal = np.arctan((efinal * np.sin(nutran))
                             / (1.0 + efinal * np.cos(nutran)))
        deltavb = np.sqrt(vtranb * vtranb + vfinal * vfinal
                          - 2.0 * vtranb * vfinal * np.cos(fpatranb - fpafinal))
        # ----------------  find transfer time of flight  ------------- }
        if etran < 0.99999:
            sinv = ((np.sqrt(1.0 - etran * etran) * np.sin(nutran))
                    / (1.0 + etran * np.cos(nutran)))
            cosv = (etran + np.cos(nutran)) / (1.0 + etran * np.cos(nutran))
            e = math.atan2(sinv, cosv)
            dttu = (np.sqrt((atran * atran * atran) / mu)
                    * (e - etran * np.sin(e) - (eainit - etran * np.sin(eainit))))
        else:
            if np.abs(etran - 1.0) < 1e-06:
                # parabolic dttu }
                pass
            else:
                # hyperbolic dttu }
                pass
    else:
        print('the one tangent burn is not possible for this case ' % ())

    #      fprintf(1, ' atran   #11.7f  #11.7f km #11.7f \n', atran, atran*re, ptran*re)
#      fprintf(1, ' etran   #11.7f  \n', etran)
#      fprintf(1, ' vinit   #11.7f  #11.7f km/s \n', vinit, vinit*velkmps)
#      fprintf(1, ' phi  #11.7f  #11.7f km/s \n', fpatranb * rad2deg, fpafinal * rad2deg)
#      fprintf(1, ' E  #11.7f  \n', e * rad2deg)
#      fprintf(1, ' vtrana  #11.7f  #11.7f km/s \n', vtrana, vtrana*velkmps)
#      fprintf(1, ' vtranb  #11.7f  #11.7f km/s \n', vtranb, vtranb*velkmps)
#      fprintf(1, ' vfinal  #11.7f  #11.7f km/s \n', vfinal, vfinal*velkmps)
    return deltava, deltavb, dttu, etran, atran, vtrana, vtranb



# ------------------------------------------------------------------------------
#
#                           procedure hohmann
#
#  this procedure calculates the delta v's for a hohmann transfer for either
#    circle to circle, or ellipse to ellipse.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    rinit       - initial position magnitude     er
#    rfinal      - final position magnitude       er
#    einit       - eccentricity of first orbit
#    efinal      - eccentricity of final orbit
#    nuinit      - true anomaly of first orbit    0 or pi rad
#    nufinal     - true anomaly of final orbit    0 or pi rad
#
#  outputs       :
#    deltava     - change in velocity at point a  er / tu
#    deltavb     - change in velocity at point b  er / tu
#    dttu        - time of flight for the trans   tu
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
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 327, alg 36, ex 6-1
#
#function [deltava, deltavb, dttu ] = hohmann (rinit, rfinal, einit, efinal, nuinit, nufinal)
# ----------------------------------------------------------------------------- }

def hohmann(rinit=None, rfinal=None, einit=None, efinal=None,
            nuinit=None, nufinal=None):
    # --------------------  initialize values   ------------------- }
    mu = 1.0

    ainit = (rinit * (1.0 + einit * np.cos(nuinit))) / (1.0 - einit * einit)
    atran = (rinit + rfinal) / 2.0
    afinal = ((rfinal * (1.0 + efinal * np.cos(nufinal)))
              / (1.0 - efinal * efinal))
    deltava = 0.0
    deltavb = 0.0
    dttu = 0.0
    if (einit < 1.0) or (efinal < 1.0):
        # -----------------  find delta v at point a  -------------- }
        vinit = np.sqrt((2.0 * mu) / rinit - (mu / ainit))
        vtrana = np.sqrt((2.0 * mu) / rinit - (mu / atran))
        deltava = np.abs(vtrana - vinit)
        # -----------------  find delta v at point b  -------------- }
        vfinal = np.sqrt((2.0 * mu) / rfinal - (mu / afinal))
        vtranb = np.sqrt((2.0 * mu) / rfinal - (mu / atran))
        deltavb = np.abs(vfinal - vtranb)
        # ----------------  find transfer time of flight  ---------- }
        dttu = np.pi * np.sqrt((atran * atran * atran) / mu)
        print(' atran   %11.7f  %11.7f km \n' % (atran, atran * re))
        print(' vinit   %11.7f  %11.7f km/s \n' % (vinit, vinit * velkmps))
        print(' vtrana  %11.7f  %11.7f km/s \n' % (vtrana, vtrana * velkmps))
        print(' vtranb  %11.7f  %11.7f km/s \n' % (vtranb, vtranb * velkmps))
        print(' vfinal  %11.7f  %11.7f km/s \n' % (vfinal, vfinal * velkmps))

    return deltava, deltavb, dttu




# ------------------------------------------------------------------------------
#
#                           procedure biellip
#
#  this procedure calculates the delta v's for a bi-elliptic transfer for either
#    circle to circle, or ellipse to ellipse.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    rinit       - initial position magnitude     er
#    r2          - interim orbit magnitude        er
#    rfinal      - final position magnitude       er
#    einit       - eccentricity of first orbit
#    efinal      - eccentricity of final orbit
#    nuinit      - true anomaly of first orbit    0 or pi rad
#    nufinal     - true anomaly of final orbit    0 or pi rad, opp of nuinit
#
#  outputs       :
#    deltava     - change in velocity at point a  er / tu
#    deltavb     - change in velocity at point b  er / tu
#    dttu        - time of flight for the trans   tu
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
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 327, alg 37, ex 6-2
#function [deltava, deltavb, deltavc, dttu ] = biellip(rinit, rb, rfinal, einit, efinal, nuinit, nufinal)
# ----------------------------------------------------------------------------- }


def biellip(rinit, rb, rfinal, einit, efinal, nuinit, nufinal):
     # --------------------  initialize values   ------------------- }
     mu = 1.0 # cannonical units

     ainit = (rinit * (1.0 + einit * math.cos(nuinit))) / (1.0 - einit * einit)
     atran1 = (rinit + rb) * 0.5
     atran2 = (rb + rfinal) * 0.5
     afinal = ((rfinal * (1.0 + efinal * math.cos(nufinal)))
               / (1.0 - efinal * efinal))

     deltava = 0.0
     deltavb = 0.0
     deltavc = 0.0
     dttu = 0.0

     if (einit < 1.0) and (efinal < 1.0):

     # -----------------  find delta v at point a  ----------------- }
         vinit = math.sqrt((2.0 * mu)/rinit - (mu/ainit))
         vtran1a = math.sqrt((2.0 * mu)/rinit - (mu/atran1))
         deltava = abs(vtran1a - vinit)

     # -----------------  find delta v at point b  ----------------- }
         vtran1b = math.sqrt((2.0 * mu)/rb - (mu/atran1))
         vtran2b = math.sqrt((2.0 * mu)/rb - (mu/atran2))
         deltavb = abs(vtran1b - vtran2b)

     # -----------------  find delta v at point c  ----------------- }
         vtran2c = math.sqrt((2.0 * mu)/rfinal - (mu/atran2))
         vfinal = math.sqrt((2.0 * mu)/rfinal - (mu/afinal))
         deltavc = abs(vfinal - vtran2c)

     # ----------------  find transfer time of flight  ------------- }
         dttu = (math.pi * math.sqrt((atran1 * atran1 * atran1)/mu)
                + math.pi * math.sqrt((atran2 * atran2 * atran2)/mu))


     #it is odd that they redefine mu here...
     mu = 398600.4415      # km3/s2
     velkmps = math.sqrt(mu / re)


     t1 = math.pi * math.sqrt((atran1 * atran1 * atran1)/1)*13.446852064
     t2 = math.pi * math.sqrt((atran2 * atran2 * atran2)/1)*13.446852064
     print(' atran1   %11.7f  %11.7f km \n' % (atran1, atran1*re))
     print(' atran2   %11.7f  %11.7f km \n' % (atran2, atran2*re))
     print(' vinit   %11.7f  %11.7f km/s \n' % (vinit, vinit*velkmps))
     print(' vtran1a  %11.7f  %11.7f km/s \n' % (vtran1a, vtran1a*velkmps))
     print(' vtran1b  %11.7f  %11.7f km/s \n' % (vtran1b, vtran1b*velkmps))
     print(' vtran2b  %11.7f  %11.7f km/s \n' % (vtran2b, vtran2b*velkmps))
     print(' vtran2c  %11.7f  %11.7f km/s \n' % (vtran2c, vtran2c*velkmps))
     print(' vfinal  %11.7f  %11.7f km/s \n' % (vfinal, vfinal*velkmps))
     print(' t1  %11.7f t2 %11.7f min \n' % (t1, t2))


     return deltava, deltavb, deltavc, dttu





# ----------------------------------------------------------------------------
#
#                           function iau80in
#
#  this function initializes the nutation matricies needed for reduction
#    calculations. the routine needs the filename of the files as input.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    iar80       - integers for fk5 1980
#    rar80       - reals for fk5 1980             rad
#
#  outputs       :
#    none
#
#  locals        :
#    convrt      - conversion factor to degrees
#    i, j         - index
#
#  coupling      :
#    none        -
#
#  references    :
#
# [iar80, rar80] = iau80in()
# ----------------------------------------------------------------------------- }


def iau80in():
    # ------------------------  implementation   -------------------
    # 0.0001" to rad
    convrt = 0.0001 * math.pi / (180*3600.0)

    fn = os.path.join(os.path.dirname(__file__), "data", "nut80.dat")

    tdat = open(fn, "r").readlines()

    iarr = np.zeros((106, 5), dtype = int)
    rarr = np.zeros((106, 5))

    ii = 0
    for rec in tdat:
      lnd = rec.strip().replace(" ", ", ").replace(", , , , ", ", ").replace(", , , ", ", ").replace(", , ", ", ")
      lid = [int(xx) for xx in lnd.split(", ")[0:5]]
      iarr[ii, 0:5] = np.array([lid])

      rid = [float(xx)*convrt for xx in lnd.split(", ")[5:10]]
      rarr[ii, 0:5] = np.array([rid])

      ii += 1

    #print("iau80 float constant is: ", convrt)


    return iarr, rarr

#
# ----------------------------------------------------------------------------
#
#                           function iau06era
#
#  this function calulates the transformation matrix that accounts for the
#    effects of sidereal time via the earth rotation angle.
#
#  author        : david vallado                  719-573-2600   16 jul 2004
#
#  inputs          description                    range / units
#    jdut1       - julian date of ut1             days
#
#  outputs       :
#    st          - transformation matrix for pef-ire
#
#  locals        :
#    tdut1       - julian centuries of ut1        days
#    era         - earth rotation angle           rad
#
#  coupling      :
#
#  references    :
#    vallado       2004, 212
#
# [st] = iau06era (jdut1)
# ----------------------------------------------------------------------------

def iau06era(jdut1=None):
    # julian centuries of ut1
    tut1d = jdut1 - 2451545.0
    era = twopi * (0.779057273264 + 1.0027378119113546 * tut1d)
    era = np.fmod(era, twopi)
    if sh.iauhelp == 'y':
      print('era%11.7f  \n' % (era * 180 / np.pi))

    # transformation matrix
    st = np.zeros((3, 3))
    st[0, 0] = np.cos(era)
    st[0, 1] = - np.sin(era)
    st[0, 2] = 0.0
    st[1, 0] = np.sin(era)
    st[1, 1] = np.cos(era)
    st[1, 2] = 0.0
    st[2, 0] = 0.0
    st[2, 1] = 0.0
    st[2, 2] = 1.0
    return st



# -----------------------------------------------------------------------------
#
#                           function iau06gst
#
#  this function finds the iau2006 greenwich sidereal time.
#
#  author        : david vallado                  719-573-2600   16 jul 2004
#
#  revisions
#
#  inputs          description                    range / units
#    jdut1       - julian date of ut1             days from 4713 bc
#    ttt         - julian centuries of tt
#    deltapsi    - change in longitude            rad
#    l           - delaunay element               rad
#    ll          - delaunay element               rad
#    f           - delaunay element               rad
#    d           - delaunay element               rad
#    omega       - delaunay element               rad
#    many others for planetary values             rad
#
#  outputs       :
#    gst         - greenwich sidereal time        0 to twopi rad
#    st          - transformation matrix
#
#  locals        :
#    temp        - temporary variable for reals   rad
#    tut1d       - days from the jan 1, 2000 12 h epoch (ut1)
#
#  coupling      :
#    iau00in     - initialize the data arrays
#
#  references    :
#    vallado       2004, 216
#
# [gst, st] = iau06gst(jdut1, ttt, deltapsi, l, l1, f, d, omega, ...
#            lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate)
# -----------------------------------------------------------------------------

def iau06gst(jdut1=None, ttt=None, deltapsi=None, opt=None):

    #added by jmb to resolve this list of vars
    l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, \
        lonsat, lonurn, lonnep, precrate = smu.fundarg(ttt, opt)

    # " to rad
    convrt = np.pi / (180.0 * 3600.0)
    axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti = iau06in()
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    ttt4 = ttt2 * ttt2
    ttt5 = ttt3 * ttt2
    # mean obliquity of the ecliptic
    epsa = (84381.406 - 46.836769 * ttt - 0.0001831 * ttt2
            + 0.0020034 * ttt3 - 5.76e-07 * ttt4 - 4.34e-08 * ttt5)

    epsa = np.fmod(epsa / 3600.0, 360.0)

    epsa = epsa * deg2rad

    #  evaluate the ee complementary terms
    gstsum0 = 0.0
    for i in range(32, -1, -1):
        tempval = agsti[i, 0] * l + agsti[i, 1] * l1 + agsti[i, 2] * f + \
                  agsti[i, 3] * d + agsti[i, 4] * omega + agsti[i, 5] * lonmer + \
                  agsti[i, 6] * lonven + agsti[i, 7] * lonear + agsti[i, 8] * lonmar + \
                  agsti[i, 9] * lonjup + agsti[i, 10] * lonsat + agsti[i, 11] * lonurn + \
                  agsti[i, 12] * lonnep + agsti[i, 13] * precrate
        gstsum0 = gstsum0 + agst[i, 0] * np.sin(tempval) + agst[i, 1] * np.cos(tempval)

    gstsum1 = 0.0
    for j in range(1): #i have no idea why this was done this way -jmb
        i = 33 + j
        tempval = agsti[i, 0] * l + agsti[i, 1] * l1 + agsti[i, 2] * f + \
                  agsti[i, 3] * d + agsti[i, 4] * omega + agsti[i, 5] * lonmer + \
                  agsti[i, 6] * lonven + agsti[i, 7] * lonear + agsti[i, 8] * lonmar + \
                  agsti[i, 9] * lonjup + agsti[i, 10] * lonsat + agsti[i, 11] * lonurn + \
                  agsti[i, 12] * lonnep + agsti[i, 13] * precrate
        gstsum1 = gstsum1 + agst[i, 0] * ttt * np.sin(tempval) + agst[i, 1] * ttt * np.cos(tempval)

    eect2000 = gstsum0 + gstsum1 * ttt

    # equation of the equinoxes
    ee2000 = deltapsi * np.cos(epsa) + eect2000

    #  earth rotation angle
    tut1d = jdut1 - 2451545.0
    era = twopi * (0.779057273264 + 1.0027378119113546 * tut1d)
    era = np.fmod(era, twopi)

    #  greenwich mean sidereal time, iau 2000.
    gmst2000 = era + (0.014506 + 4612.156534 * ttt + 1.3915817 * ttt2
                      - 4.4e-07 * ttt3 + 2.9956e-05 * ttt4
                      + 3.68e-08 * ttt5) * convrt

    gst = gmst2000 + ee2000

    if sh.iauhelp == 'y':
        print('meanobl %11.7f getsum %11.7f %11.7f eect %11.7f  \n'
              % (epsa * 180 / np.pi, gstsum0 * 180 / np.pi,
                 gstsum1 * 180 / np.pi, eect2000 * 180 / np.pi))
        print('ee2000 %11.7f gmst2000 %11.7f gst %11.7f  \n'
              % (ee2000 * 180 / np.pi, gmst2000 * 180 / np.pi,
                  gst * 180 / np.pi))

    # transformation matrix
    st = np.zeros((3, 3))
    st[0, 0] = np.cos(gst)
    st[0, 1] = - np.sin(gst)
    st[0, 2] = 0.0
    st[1, 0] = np.sin(gst)
    st[1, 1] = np.cos(gst)
    st[1, 2] = 0.0
    st[2, 0] = 0.0
    st[2, 1] = 0.0
    st[2, 2] = 1.0
    return gst, st

# ----------------------------------------------------------------------------
#
#                           function iau06pna
#
#  this function calulates the transformation matrix that accounts for the
#    effects of precession-nutation in the iau2000a theory.
#
#  author        : david vallado                  719-573-2600   16 jul 2004
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    ttt         - julian centuries of tt
#
#  outputs       :
#    nut         - transformation matrix for ire-gcrf
#    deltapsi    - change in longitude            rad
#    l           - delaunay element               rad
#    ll          - delaunay element               rad
#    f           - delaunay element               rad
#    d           - delaunay element               rad
#    omega       - delaunay element               rad
#    many others for planetary values             rad
#
#  locals        :
#    x           - coordinate                     rad
#    y           - coordinate                     rad
#    s           - coordinate                     rad
#    axs0        - real coefficients for x        rad
#    a0xi        - integer coefficients for x
#    ays0        - real coefficients for y        rad
#    a0yi        - integer coefficients for y
#    ass0        - real coefficients for s        rad
#    a0si        - integer coefficients for s
#    apn         - real coefficients for nutation rad
#    apni        - integer coefficients for nutation
#    appl        - real coefficients for planetary nutation rad
#    appli       - integer coefficients for planetary nutation
#    ttt2, ttt3,  - powers of ttt
#    deltaeps    - change in obliquity            rad
#
#  coupling      :
#    iau00in     - initialize the arrays
#    fundarg     - find the fundamental arguments
#    precess     - find the precession quantities
#
#  references    :
#    vallado       2004, 212-214
#
# [ deltapsi, pnb, nut, l, l1, f, d, omega, ...
#   lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
# ] = iau06pna (ttt)
# ----------------------------------------------------------------------------

def iau06pna(ttt=None):
    showit = 'y'
    # " to rad
    convrt = np.pi / (180.0 * 3600.0)
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    ttt4 = ttt2 * ttt2
    ttt5 = ttt3 * ttt2
    # obtain data for calculations from the 2000a theory
    opt = '06'

    l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, \
        lonsat, lonurn, lonnep, precrate = smu.fundarg(ttt, opt)
    # ---- obtain data coefficients
    axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti = iau06in()
    #        [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, ape, apei, agst, agsti] = iau06in
    pnsum = 0.0
    ensum = 0.0
    for i in range(677, -1, -1):
        tempval = (apni[i, 0] * l + apni[i, 1] * l1 + apni[i, 2]
                   * f + apni[i, 3] * d + apni[i, 4] * omega)
        tempval = np.mod(tempval, 2 * np.pi)
        #            pnsum = pnsum + (apn[i, 0] + apn[i, 1]*ttt) * sin(tempval) ...
        #                          + (apn[i, 4] + apn[i, 5]*ttt) * cos(tempval)
        #            ensum = ensum + (apn[i, 2] + apn[i, 3]*ttt) * cos(tempval) ...
        #                          + (apn[i, 6] + apn[i, 7]*ttt) * sin(tempval)
        # iers doesn't include the last few terms
        pnsum = (pnsum + (apn[i, 0] + apn[i, 1] * ttt) * np.sin(tempval)
                 + (apn[i, 4]) * np.cos(tempval))
        ensum = (ensum + (apn[i, 2] + apn[i, 3] * ttt) * np.cos(tempval)
                 + (apn[i, 6]) * np.sin(tempval))

    pplnsum = 0.0
    eplnsum = 0.0
    # data file is already reveresed
    for i in range(687):
        tempval = appli[i, 0] * l + appli[i, 1] * l1 + appli[i, 2] * f + \
                  appli[i, 3] * d + appli[i, 4] * omega + appli[i, 5] * lonmer + \
                  appli[i, 6] * lonven + appli[i, 7] * lonear + appli[i, 8] * lonmar + \
                  appli[i, 9] * lonjup + appli[i, 10] * lonsat + appli[i, 11] * lonurn + \
                  appli[i, 12] * lonnep + appli[i, 13] * precrate
        pplnsum = pplnsum + appl[i, 0] * np.sin(tempval) + appl[i, 1] * np.cos(tempval)
        eplnsum = eplnsum + appl[i, 2] * np.sin(tempval) + appl[i, 3] * np.cos(tempval)

    #  add planetary and luni-solar components.
    deltapsi = pnsum + pplnsum

    deltaeps = ensum + eplnsum
    if showit == 'y':
        print('dpsi %11.7f deltaeps %11.7f \n' % (deltapsi / convrt, deltaeps / convrt))

    # iau2006 approach - does not seem to be correct, close though
# looks like they still use the iau2000a method and adjust
#        pnsum = 0.0
#        # data file is not not reveresed
#        for i = 1358 : -1 : 1
#            tempval = apni[i, 0]*l + apni[i, 1]*l1 + apni[i, 2]*f + apni[i, 3]*d + apni[i, 4]*omega + ...
#                      apni[i, 5]*lonmer  + apni[i, 6]*lonven  + apni[i, 7]*lonear  + apni[i, 8]*lonmar + ...
#                      apni[i, 9]*lonjup + apni[i, 10]*lonsat + apni[i, 11]*lonurn + apni[i, 12]*lonnep + apni[i, 13]*precrate
#            if i > 1320
#                pnsum = pnsum + (apn[i, 0] * sin(tempval) + apn[i, 1] * cos(tempval)) * ttt  #note that sin and cos are reveresed ebtween n and e
#            else
#                pnsum = pnsum + apn[i, 0] * sin(tempval) + apn[i, 1] * cos(tempval)
#            end
#          end

    #        ensum = 0.0
#        # data file is not reveresed
#        for i = 1056 : -1 : 1
#            tempval = apei[i, 0]*l + apei[i, 1]*l1 + apei[i, 2]*f + apei[i, 3]*d + apei[i, 4]*omega + ...
#                      apei[i, 5]*lonmer  + apei[i, 6]*lonven  + apei[i, 7]*lonear  + apei[i, 8]*lonmar + ...
#                      apei[i, 9]*lonjup + apei[i, 10]*lonsat + apei[i, 11]*lonurn + apei[i, 12]*lonnep + apei[i, 13]*precrate
#            if i > 1037
#                ensum = ensum + (ape[i, 0] * cos(tempval) + ape[i, 1] * sin(tempval)) * ttt
#            else
#                ensum = ensum + ape[i, 0] * cos(tempval) + ape[i, 1] * sin(tempval)
#            end
#          end
#          #  add planetary and luni-solar components.
#        deltapsi = pnsum  # rad
#        deltaeps = ensum

    # iau2006 corrections to the iau2000a
    j2d = - 2.7774e-06 * ttt * convrt

    deltapsi = deltapsi + deltapsi * (4.697e-07 + j2d)

    deltaeps = deltaeps + deltaeps * j2d
    if showit == 'y':
        print('dpsi %11.7f deltaeps %11.7f \n'
              % (deltapsi / convrt, deltaeps / convrt))

    prec, psia, wa, ea, xa = precess(ttt, '06')
    if showit == 'y':
        print('prec iau 06 \n' % ())
        print((prec))

    oblo = 84381.406 * convrt

    # ----------------- find nutation matrix ----------------------
# mean to true
    a1 = smu.rot1mat(ea + deltaeps)
    a2 = smu.rot3mat(deltapsi)
    a3 = smu.rot1mat(- ea)
    # j2000 to date (precession)
    a4 = smu.rot3mat(- xa)
    a5 = smu.rot1mat(wa)
    a6 = smu.rot3mat(psia)
    a7 = smu.rot1mat(- oblo)
    # icrs to j2000
    a8 = smu.rot1mat(- 0.0068192 * convrt)
    a9 = smu.rot2mat(0.041775 * np.sin(oblo) * convrt)
    #      a9 = smu.rot2mat(0.0166170*convrt)
    a10 = smu.rot3mat(0.0146 * convrt)
    if sh.iauhelp == 'y':
        print('p e %11.7f  %11.7f  \n'
              % (pnsum * 180 / np.pi, ensum * 180 / np.pi))
        print('p e %11.7f  %11.7f  \n'
              % (pnsum * 3600 * 180 / np.pi, ensum * 3600 * 180 / np.pi))
        print('p e %11.7f  %11.7f  \n'
              % (pplnsum * 180 / np.pi, eplnsum * 180 / np.pi))
        print('p e %11.7f  %11.7f  \n'
              % (pplnsum * 3600 * 180 / np.pi, eplnsum * 3600 * 180 / np.pi))
        print('dpsi %11.7f deps %11.7f  \n'
              % (deltapsi * 180 / np.pi, deltaeps * 180 / np.pi))
        print('dpsi %11.7f deps %11.7f  \n'
              % (deltapsi * 3600 * 180 / np.pi, deltaeps * 3600 * 180 / np.pi))
        print('psia %11.7f wa %11.7f ea %11.7f xa %11.7f  \n'
              % (psia * 180 / np.pi, wa * 180 / np.pi,
                 ea * 180 / np.pi, xa * 180 / np.pi))
        print('psia %11.7f wa %11.7f ea %11.7f xa %11.7f  \n'
              % (psia * 3600 * 180 / np.pi,
                 wa * 3600 * 180 / np.pi, ea * 3600 * 180 / np.pi,
                 xa * 3600 * 180 / np.pi))
        # temp1 = a7*a6*a5*a4
# temp2 = a3*a2*a1
# temp3 = a10*a9*a8

    pnb = a10 * a9 * a8 * a7 * a6 * a5 * a4 * a3 * a2 * a1
    prec = a7 * a6 * a5 * a4
    if showit == 'y':
        print('prec iau 06a alt \n' % ())
        print((prec))

    nut = a3 * a2 * a1
    if showit == 'y':
        print('nut iau 06a \n' % ())
        print((nut))

    frb = a10 * a9 * a8
    if showit == 'y':
        print('frb iau 06a \n' % ())
        print((frb))

    return deltapsi, pnb, prec, nut, l, l1, f, d, omega, \
           lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate

# ----------------------------------------------------------------------------
#
#                           function iau06pnb
#
#  this function calulates the transformation matrix that accounts for the
#    effects of precession-nutation in the iau2000b theory.
#
#  author        : david vallado                  719-573-2600   16 jul 2004
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    ttt         - julian centuries of tt
#
#  outputs       :
#    nut         - transformation matrix for ire-gcrf
#    deltapsi    - change in longitude            rad
#    l           - delaunay element               rad
#    ll          - delaunay element               rad
#    f           - delaunay element               rad
#    d           - delaunay element               rad
#    omega       - delaunay element               rad
#    many others
#
#  locals        :
#    x           - coordinate                     rad
#    y           - coordinate                     rad
#    s           - coordinate                     rad
#    axs0        - real coefficients for x        rad
#    a0xi        - integer coefficients for x
#    ays0        - real coefficients for y        rad
#    a0yi        - integer coefficients for y
#    ass0        - real coefficients for s        rad
#    a0si        - integer coefficients for s
#    apn         - real coefficients for nutation rad
#    apni        - integer coefficients for nutation
#    appl        - real coefficients for planetary nutation rad
#    appli       - integer coefficients for planetary nutation
#    ttt2, ttt3,  - powers of ttt
#    deltaeps    - change in obliquity            rad
#
#  coupling      :
#    iau00in     - initialize the arrays
#    fundarg     - find the fundamental arguments
#    precess     - find the precession coefficients
#
#  references    :
#    vallado       2004, 212-214
#
# [ deltapsi, pnb, nut, l, l1, f, d, omega, ...
#   lonmer, lonven, lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate ...
# ] = iau06pnb (ttt)
# ----------------------------------------------------------------------------

def iau06pnb(ttt=None):
    # " to rad
    convrt = np.pi / (180.0 * 3600.0)
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    ttt4 = ttt2 * ttt2
    ttt5 = ttt3 * ttt2
    # obtain data for calculations form the 2000b theory
    opt = '06' #changed by jmb

    l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, \
        lonurn, lonnep, precrate = smu.fundarg(ttt, opt)
    # ---- obtain data coefficients
    axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti = iau06in()
    #        [axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, ape, apei, agst, agsti] = iau06in

    pnsum = 0.0
    ensum = 0.0
    for i in range(76, -1, -1):
        tempval = (apni[i, 0] * l + apni[i, 1] * l1 + apni[i, 2] * f
                   + apni[i, 3] * d + apni[i, 4] * omega)
        pnsum = (pnsum + (apn[i, 0] + apn[i, 1] * ttt) * np.sin(tempval)
                 + (apn[i, 4] + apn[i, 5] * ttt) * np.cos(tempval))
        ensum = (ensum + (apn[i, 2] + apn[i, 3] * ttt) * np.cos(tempval)
                 + (apn[i, 6] + apn[i, 7] * ttt) * np.sin(tempval))
        #             pnsum = pnsum + (apn[i, 0] + apn[i, 1]*ttt) * sin(tempval) ...
#                           + (apn[i, 4]) * cos(tempval)
#             tempval = apei[i, 0]*l + apei[i, 1]*l1 + apei[i, 2]*f + apei[i, 3]*d + apei[i, 4]*omega
#             ensum = ensum + (ape[i, 2] + ape[i, 3]*ttt) * cos(tempval) ...
#                           + (ape[i, 6]) * sin(tempval)

    # iau2006 approach - does not seem to be correct
# looks like they still use the iau2000a method and adjust
#        pnsum = 0.0
#        # data file is not already reveresed
#        for i = 77 : -1 : 1
#            tempval = apni[i, 0]*l + apni[i, 1]*l1 + apni[i, 2]*f + apni[i, 3]*d + apni[i, 4]*omega
#            if i > 1320
#                pnsum = pnsum + (apn[i, 0] * sin(tempval) + apn[i, 1] * cos(tempval)) * ttt
#            else
#                pnsum = pnsum + apn[i, 0] * sin(tempval) + apn[i, 1] * cos(tempval)
#            end
#          end

    #        ensum = 0.0
#        # data file is already reveresed
#        for i = 77 : -1 : 1
#            tempval = apei[i, 0]*l + apei[i, 1]*l1 + apei[i, 2]*f + apei[i, 3]*d + apei[i, 4]*omega
#            if i > 1037
#                ensum = ensum + (ape[i, 0] * cos(tempval) + ape[i, 1] * sin(tempval)) * ttt
#            else
#                ensum = ensum + ape[i, 0] * cos(tempval) + ape[i, 1] * sin(tempval)
#            end
#          end
#          #  add planetary and luni-solar components.
#        deltapsi = pnsum  # rad
#        deltaeps = ensum

    # ------ form the planetary arguments
    pplnsum = - 0.000135 * convrt

    eplnsum = 0.000388 * convrt
    #  add planetary and luni-solar components.
    deltapsi = pnsum + pplnsum
    deltaeps = ensum + eplnsum
    prec, psia, wa, ea, xa = precess(ttt, '06')
    oblo = 84381.406 * convrt

    # or 84381.406????

    # ----------------- find nutation matrix ----------------------
# mean to true
    a1 = smu.rot1mat(ea + deltaeps)
    a2 = smu.rot3mat(deltapsi)
    a3 = smu.rot1mat(- ea)
    # j2000 to date (precession)
    a4 = smu.rot3mat(- xa)
    a5 = smu.rot1mat(wa)
    a6 = smu.rot3mat(psia)
    a7 = smu.rot1mat(- oblo)
    # icrs to j2000
    a8 = smu.rot1mat(- 0.0068192 * convrt)
    a9 = smu.rot2mat(0.041775 * np.sin(oblo) * convrt)
    #      a9 = smu.rot2mat(0.0166170*convrt)
    a10 = smu.rot3mat(0.0146 * convrt)
    pnb = a10 * a9 * a8 * a7 * a6 * a5 * a4 * a3 * a2 * a1
    prec = a10 * a9 * a8 * a7 * a6 * a5 * a4
    nut = a3 * a2 * a1
    if sh.iauhelp == 'y':
        print('p e %11.7f  %11.7f  \n' % (pnsum * 180 / np.pi, ensum * 180 / np.pi))
        print('dpsi %11.7f deps %11.7f  \n'
              % (deltapsi * 180 / np.pi, deltaeps * 180 / np.pi))
        print('psia %11.7f wa %11.7f ea %11.7f xa %11.7f  \n'
              % (psia * 180 / np.pi, wa * 180 / np.pi, ea * 180 / np.pi,
                 xa * 180 / np.pi))

    # -------------- these are extra not needed for pnb
    if (sh.iaupnhelp == 'y'):
        p = psia + (deltapsi * np.sin(ea) * np.cos(xa)
                    - deltaeps * np.sin(xa)) / np.sin(wa)
        w = wa + deltapsi * np.sin(ea) * np.sin(xa) + deltaeps * np.cos(xa)
        xbar = np.sin(w) * np.sin(p)
        ybar = - np.sin(oblo) * np.cos(w) + np.cos(oblo) * np.sin(w) * np.cos(p)
        x = xbar + (- 0.016617 + 0.0146 * ybar) * convrt
        #            x = xbar + (-0.0417750 + 0.01460*ybar)*convrt  # rad
        y = ybar + (- 0.0068192 - 0.0146 * xbar) * convrt
        # -------- now find a
        a = 0.5 + 0.125 * (x * x + y * y)
        # -------- now find s
        ssum0 = 0.0
        for i in range(32, -1, -1):
            tempval = (a0si[i, 0] * l + a0si[i, 1] * l1
                       + a0si[i, 2] * f + a0si[i, 3] * d + a0si[i, 4] * omega)
            ssum0 = (ssum0 + ass0[i, 0] * np.sin(tempval)
                     + ass0[i, 1] * np.cos(tempval))
        ssum1 = 0.0
        for j in range(2, -1, -1):
            i = 33 + j
            tempval = (a0si[i, 0] * l + a0si[i, 1] * l1
                       + a0si[i, 2] * f + a0si[i, 3] * d + a0si[i, 4] * omega)
            ssum1 = ssum1 + ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)
        ssum2 = 0.0
        for j in range(24, -1, -1):
            i = 33 + 3 + j
            tempval = (a0si[i, 0] * l + a0si[i, 1] * l1
                       + a0si[i, 2] * f + a0si[i, 3] * d + a0si[i, 4] * omega)
            ssum2 = (ssum2 + ass0[i, 0] * np.sin(tempval)
                     + ass0[i, 1] * np.cos(tempval))
        ssum3 = 0.0
        for j in range(3, -1, -1):
            i = 33 + 3 + 25 + j
            tempval = (a0si[i, 0] * l + a0si[i, 1] * l1
                       + a0si[i, 2] * f + a0si[i, 3] * d + a0si[i, 4] * omega)
            ssum3 = (ssum3 + ass0[i, 0] * np.sin(tempval)
                     + ass0[i, 1] * np.cos(tempval))
        ssum4 = 0.0
        for j in range(1):
            i = 33 + 3 + 25 + 4 + j
            tempval = (a0si[i, 0] * l + a0si[i, 1] * l1 + a0si[i, 2]
                       * f + a0si[i, 3] * d + a0si[i, 4] * omega)
            ssum4 = (ssum4 + ass0[i, 0] * np.sin(tempval)
                     + ass0[i, 1] * np.cos(tempval))
        s = (9.4e-05 + 0.00380835 * ttt - 0.00011994 * ttt2
             - 0.07257409 * ttt3 + 2.77e-05 * ttt4 + 1.561e-05 * ttt5)
        s = (- x * y * 0.5 + s * convrt + ssum0 + ssum1 * ttt + ssum2 * ttt2
             + ssum3 * ttt3 + ssum4 * ttt4)
        if sh.iauhelp == 'y':
            print('00pnb  x  %14.12f" y  %14.12f" s %14.12f" a %14.12fdeg \n'
                  % (x / deg2rad * 3600, y / deg2rad * 3600,
                     s / deg2rad * 3600, a / deg2rad))
            #                fprintf(1, 'p #11.7f w #11.7f xbar #11.7f ybar #11.7f deg \n', p*180/pi, w*180/pi, xbar*180/pi, ybar*180/pi)


    return deltapsi, pnb, prec, nut, l, l1, f, d, omega, lonmer, lonven, \
        lonear, lonmar, lonjup, lonsat, lonurn, lonnep, precrate

# ----------------------------------------------------------------------------
#
#                           function iau06xys
#
#  this function calulates the transformation matrix that accounts for the
#    effects of precession-nutation in the iau2006 theory.
#
#  author        : david vallado                  719-573-2600   16 jul 2004
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    ttt         - julian centuries of tt
#    ddx         - eop correction for x           rad
#    ddy         - eop correction for y           rad
#
#  outputs       :
#    nut         - transformation matrix for tirs-gcrf
#    x           - coordinate of cip              rad
#    y           - coordinate of cip              rad
#    s           - coordinate                     rad
#
#  locals        :
#    axs0        - real coefficients for x        rad
#    a0xi        - integer coefficients for x
#    ays0        - real coefficients for y        rad
#    a0yi        - integer coefficients for y
#    ass0        - real coefficients for s        rad
#    a0si        - integer coefficients for s
#    apn         - real coefficients for nutation rad
#    apni        - integer coefficients for nutation
#    appl        - real coefficients for planetary nutation rad
#    appli       - integer coefficients for planetary nutation
#    ttt2, ttt3,  - powers of ttt
#    l           - delaunay element               rad
#    ll          - delaunay element               rad
#    f           - delaunay element               rad
#    d           - delaunay element               rad
#    omega       - delaunay element               rad
#    deltaeps    - change in obliquity            rad
#    many others
#
#  coupling      :
#    iau00in     - initialize the arrays
#    fundarg     - find the fundamental arguments
#
#  references    :
#    vallado       2004, 212-214
#
# [x, y, s, nut] = iau06xys (ttt, ddx, ddy)
# ----------------------------------------------------------------------------

def iau06xys(ttt=None, ddx=None, ddy=None):
    # " to rad
    convrt = np.pi / (180.0 * 3600.0)
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt
    ttt4 = ttt2 * ttt2
    ttt5 = ttt3 * ttt2
    axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti = iau06in()
    opt = '06'

    l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, \
        lonurn, lonnep, precrate = smu.fundarg(ttt, opt)
    # fprintf(1, '\ndelauany #11.7f #11.7f #11.7f #11.7f #11.7f \n planetary #11.7f #11.7f #11.7f #11.7f #11.7f #11.7f #11.7f #11.7f #11.7f \n', ...
#     l/deg2rad, l1/deg2rad, f/deg2rad, d/deg2rad, omega/deg2rad, lonmer/deg2rad, lonven/deg2rad, ...
#     lonear/deg2rad, lonmar/deg2rad, lonjup/deg2rad, lonsat/deg2rad, lonurn/deg2rad, lonnep/deg2rad, precrate/deg2rad)

    # ---------------- first find x
# the iers code puts the constants in here, however
# don't sum constants in here because they're larger than the last few terms
    xsum0 = 0.0
    for i in range(1305, -1, - 1):
        tempval = (a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f
                   + a0xi[i, 3] * d + a0xi[i, 4] * omega + a0xi[i, 5] * lonmer
                   + a0xi[i, 6] * lonven + a0xi[i, 7] * lonear
                   + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup
                   + a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn
                   + a0xi[i, 12] * lonnep + a0xi[i, 13] * precrate)
        xsum0 = xsum0 + axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    xsum1 = 0.0
    # note that the index changes here to j. this is because the a0xi etc
# indicies go from 1 to 1600, but there are 5 groups. the i index counts through each
# calculation, and j takes care of the individual summations. note that
# this same process is used for y and s.
    for j in range(252, -1, - 1):
        i = 1306 + j
        tempval = (a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f
                   + a0xi[i, 3] * d + a0xi[i, 4] * omega + a0xi[i, 5] * lonmer
                   + a0xi[i, 6] * lonven + a0xi[i, 7] * lonear
                   + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup
                   + a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn
                   + a0xi[i, 12] * lonnep + a0xi[i, 13] * precrate)
        xsum1 = xsum1 + axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    xsum2 = 0.0
    for j in range(35, - 1, - 1):
        i = 1306 + 253 + j
        tempval = (a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f
                   + a0xi[i, 3] * d + a0xi[i, 4] * omega + a0xi[i, 5] * lonmer
                   + a0xi[i, 6] * lonven + a0xi[i, 7] * lonear
                   + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup
                   + a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn
                   + a0xi[i, 12] * lonnep + a0xi[i, 13] * precrate)
        xsum2 = xsum2 + axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    xsum3 = 0.0
    for j in range(3, - 1, - 1):
        i = 1306 + 253 + 36 + j
        tempval = (a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f
                   + a0xi[i, 3] * d + a0xi[i, 4] * omega + a0xi[i, 5] * lonmer
                   + a0xi[i, 6] * lonven + a0xi[i, 7] * lonear
                   + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup
                   + a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn
                   + a0xi[i, 12] * lonnep + a0xi[i, 13] * precrate)
        xsum3 = xsum3 + axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    xsum4 = 0.0
    for j in range(1):
        i = 1306 + 253 + 36 + 4 + j
        tempval = (a0xi[i, 0] * l + a0xi[i, 1] * l1 + a0xi[i, 2] * f
                   + a0xi[i, 3] * d + a0xi[i, 4] * omega + a0xi[i, 5] * lonmer
                   + a0xi[i, 6] * lonven + a0xi[i, 7] * lonear
                   + a0xi[i, 8] * lonmar + a0xi[i, 9] * lonjup
                   + a0xi[i, 10] * lonsat + a0xi[i, 11] * lonurn
                   + a0xi[i, 12] * lonnep + a0xi[i, 13] * precrate)
        xsum4 = xsum4 + axs0[i, 0] * np.sin(tempval) + axs0[i, 1] * np.cos(tempval)

    x = - 0.016617 + 2004.191898 * ttt - 0.4297829 * ttt2 - 0.19861834 * ttt3 - 7.578e-06 * ttt4 + 5.9285e-06 * ttt5

    x = x * convrt + xsum0 + xsum1 * ttt + xsum2 * ttt2 + xsum3 * ttt3 + xsum4 * ttt4

    if sh.iauhelp == 'y':
        print('xys x %14.12f  %14.12f  %14.12f  %14.12f  %14.12f \n' % (xsum0 / deg2rad, xsum1 / deg2rad, xsum2 / deg2rad, xsum3 / deg2rad, xsum4 / deg2rad))

    # ---------------- now find y
    ysum0 = 0.0
    for i in range(961, - 1, - 1):
        tempval = (a0yi[i, 0] * l + a0yi[i, 1] * l1
                   + a0yi[i, 2] * f + a0yi[i, 3] * d + a0yi[i, 4] * omega
                   + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven
                   + a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar
                   + a0yi[i, 9] * lonjup + a0yi[i, 10] * lonsat
                   + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep
                   + a0yi[i, 13] * precrate)
        ysum0 = ysum0 + ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    ysum1 = 0.0
    for j in range(276, - 1, - 1):
        i = 962 + j
        tempval = (a0yi[i, 0] * l + a0yi[i, 1] * l1
                   + a0yi[i, 2] * f + a0yi[i, 3] * d + a0yi[i, 4] * omega
                   + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven
                   + a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar
                   + a0yi[i, 9] * lonjup + a0yi[i, 10] * lonsat
                   + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep
                   + a0yi[i, 13] * precrate)
        ysum1 = ysum1 + ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    ysum2 = 0.0
    for j in range(29, - 1, - 1):
        i = 962 + 277 + j
        tempval = (a0yi[i, 0] * l + a0yi[i, 1] * l1
                   + a0yi[i, 2] * f + a0yi[i, 3] * d + a0yi[i, 4] * omega
                   + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven
                   + a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar
                   + a0yi[i, 9] * lonjup + a0yi[i, 10] * lonsat
                   + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep
                   + a0yi[i, 13] * precrate)
        ysum2 = ysum2 + ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    ysum3 = 0.0
    for j in range(4, - 1, - 1):
        i = 962 + 277 + 30 + j
        tempval = (a0yi[i, 0] * l + a0yi[i, 1] * l1
                   + a0yi[i, 2] * f + a0yi[i, 3] * d + a0yi[i, 4] * omega
                   + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven
                   + a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar
                   + a0yi[i, 9] * lonjup + a0yi[i, 10] * lonsat
                   + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep
                   + a0yi[i, 13] * precrate)
        ysum3 = ysum3 + ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    ysum4 = 0.0
    for j in range(1):
        i = 962 + 277 + 30 + 5 + j
        tempval = (a0yi[i, 0] * l + a0yi[i, 1] * l1
                   + a0yi[i, 2] * f + a0yi[i, 3] * d + a0yi[i, 4] * omega
                   + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven
                   + a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar
                   + a0yi[i, 9] * lonjup + a0yi[i, 10] * lonsat
                   + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep
                   + a0yi[i, 13] * precrate)
        ysum4 = ysum4 + ays0[i, 0] * np.sin(tempval) + ays0[i, 1] * np.cos(tempval)

    y = (- 0.006951 - 0.025896 * ttt - 22.4072747 * ttt2
         + 0.00190059 * ttt3 + 0.001112526 * ttt4 + 1.358e-07 * ttt5)
    y = (y * convrt + ysum0 + ysum1 * ttt + ysum2 * ttt2
         + ysum3 * ttt3 + ysum4 * ttt4)

    if sh.iauhelp == 'y':
        print('xys y %14.12f  %14.12f  %14.12f  %14.12f  %14.12f \n'
              % (ysum0 / deg2rad, ysum1 / deg2rad, ysum2 / deg2rad,
                 ysum3 / deg2rad, ysum4 / deg2rad))

    # ---------------- now find s
    ssum0 = 0.0
    for i in range(32, - 1, - 1):
        tempval = (a0si[i, 0] * l + a0si[i, 1] * l1
                   + a0si[i, 2] * f + a0si[i, 3] * d + a0si[i, 4] * omega
                   + a0si[i, 5] * lonmer + a0si[i, 6] * lonven
                   + a0si[i, 7] * lonear + a0si[i, 8] * lonmar
                   + a0si[i, 9] * lonjup + a0si[i, 10] * lonsat
                   + a0si[i, 11] * lonurn + a0si[i, 12] * lonnep
                   + a0si[i, 13] * precrate)
        ssum0 = ssum0 + ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    ssum1 = 0.0
    for j in range(2, - 1, - 1):
        i = 33 + j
        tempval = (a0yi[i, 0] * l + a0yi[i, 1] * l1
                   + a0yi[i, 2] * f + a0yi[i, 3] * d + a0yi[i, 4] * omega
                   + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven
                   + a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar
                   + a0yi[i, 9] * lonjup + a0yi[i, 10] * lonsat
                   + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep
                   + a0yi[i, 13] * precrate)
        ssum1 = ssum1 + ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    ssum2 = 0.0
    for j in range(24, - 1, - 1):
        i = 33 + 3 + j
        tempval = (a0yi[i, 0] * l + a0yi[i, 1] * l1
                   + a0yi[i, 2] * f + a0yi[i, 3] * d + a0yi[i, 4] * omega
                   + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven
                   + a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar
                   + a0yi[i, 9] * lonjup + a0yi[i, 10] * lonsat
                   + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep
                   + a0yi[i, 13] * precrate)
        ssum2 = ssum2 + ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    ssum3 = 0.0
    for j in range(3, - 1, - 1):
        i = 33 + 3 + 25 + j
        tempval = (a0yi[i, 0] * l + a0yi[i, 1] * l1
                   + a0yi[i, 2] * f + a0yi[i, 3] * d + a0yi[i, 4] * omega
                   + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven
                   + a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar
                   + a0yi[i, 9] * lonjup + a0yi[i, 10] * lonsat
                   + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep
                   + a0yi[i, 13] * precrate)
        ssum3 = ssum3 + ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    ssum4 = 0.0
    for j in range(1):
        i = 33 + 3 + 25 + 4 + j
        tempval = (a0yi[i, 0] * l + a0yi[i, 1] * l1
                   + a0yi[i, 2] * f + a0yi[i, 3] * d + a0yi[i, 4] * omega
                   + a0yi[i, 5] * lonmer + a0yi[i, 6] * lonven
                   + a0yi[i, 7] * lonear + a0yi[i, 8] * lonmar
                   + a0yi[i, 9] * lonjup + a0yi[i, 10] * lonsat
                   + a0yi[i, 11] * lonurn + a0yi[i, 12] * lonnep
                   + a0yi[i, 13] * precrate)
        ssum4 = ssum4 + ass0[i, 0] * np.sin(tempval) + ass0[i, 1] * np.cos(tempval)

    s = 9.4e-05 + 0.00380865 * ttt - 0.00012268 * ttt2 - 0.07257411 * ttt3 + 2.798e-05 * ttt4 + 1.562e-05 * ttt5

    #            + 0.00000171*ttt*sin(omega) + 0.00000357*ttt*cos(2.0*omega) ...
#            + 0.00074353*ttt2*sin(omega) + 0.00005691*ttt2*sin(2.0*(f-d+omega)) ...
#            + 0.00000984*ttt2*sin(2.0*(f+omega)) - 0.00000885*ttt2*sin(2.0*omega)
    s = - x * y * 0.5 + s * convrt + ssum0 + ssum1 * ttt + ssum2 * ttt2 + ssum3 * ttt3 + ssum4 * ttt4

    if sh.iauhelp == 'y':
        print('06xys before x  %14.12f y  %14.12f s %14.12f a %14.12f rad \n'
              % (x, y, s, a))
        print('xys s %14.12f  %14.12f  %14.12f  %14.12f  %14.12f \n'
              % (ssum0 / deg2rad, ssum1 / deg2rad, ssum2 / deg2rad,
                 ssum3 / deg2rad, ssum4 / deg2rad))

    # add corrections if available
    x = x + ddx
    y = y + ddy
    # ---------------- now find a
    a = 0.5 + 0.125 * (x * x + y * y)

    if sh.iauhelp == 'y':
      print('06xys  x  %14.12f y  %14.12f s %14.12f a %14.12f rad \n'
          % (x, y, s, a))
      print('06xys  x  %14.12f y  %14.12f s %14.12f a %14.12f deg \n'
          % (x / deg2rad, y / deg2rad, s / deg2rad, a / deg2rad))
      print('06xys  x  %14.12f" y  %14.12f" s %14.12f" a %14.12fdeg \n'
          % (x / deg2rad * 3600, y / deg2rad * 3600, s / deg2rad * 3600,
             a / deg2rad))


    # ----------------- find nutation matrix ----------------------
    nut1 = np.zeros((3, 3))
    nut1[0, 0] = 1.0 - a * x * x
    nut1[0, 1] = - a * x * y
    nut1[0, 2] = x
    nut1[1, 0] = - a * x * y
    nut1[1, 1] = 1.0 - a * y * y
    nut1[1, 2] = y
    nut1[2, 0] = - x
    nut1[2, 1] = - y
    nut1[2, 2] = 1.0 - a * (x * x + y * y)
    #nut1

    nut2 = np.eye(3)
    nut2[0, 0] = np.cos(s)
    nut2[1, 1] = np.cos(s)
    nut2[0, 1] = np.sin(s)
    nut2[1, 0] = - np.sin(s)
    nut = nut1 * nut2
    #       the matrix appears to be orthogonal now, so the extra processing is not needed.
#        if (x ~= 0.0) && (y ~= 0.0)
#            e = atan2(y, x)
#          else
#            e = 0.0
#          end
#        d = atan(sqrt((x^2 + y^2) / (1.0-x^2-y^2)))
#        nut1 = smu.rot3mat(-e)*smu.rot2mat(-d)*smu.rot3mat(e+s)

    return x, y, s, nut







# ----------------------------------------------------------------------------
#
#                           function nutation
#
#  this function calulates the transformation matrix that accounts for the
#    effects of nutation.
#
#  author        : david vallado                  719-573-2600   27 jun 2002
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    ttt         - julian centuries of tt
#    ddpsi       - delta psi correction to gcrf   rad
#    ddeps       - delta eps correction to gcrf   rad
#
#  outputs       :
#    deltapsi    - nutation angle                 rad
#    trueeps     - true obliquity of the ecliptic rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       -                                rad
#    nut         - transformation matrix for tod - mod
#
#  locals        :
#    iar80       - integers for fk5 1980
#    rar80       - reals for fk5 1980
#    ttt2        - ttt squared
#    ttt3        - ttt cubed
#    l           -                                rad
#    ll          -                                rad
#    f           -                                rad
#    d           -                                rad
#    deltaeps    - change in obliquity            rad
#
#  coupling      :
#    fundarg     - find fundamental arguments
#
#  references    :
#    vallado       2013, 224-226
#
# [deltapsi, trueeps, meaneps, omega, nut] = nutation  (ttt, ddpsi, ddeps)
# ----------------------------------------------------------------------------


def nutation (ttt, ddpsi, ddeps):

        iar80, rar80 = iau80in()  # coeff in deg

        # ---- determine coefficients for iau 1980 nutation theory ----
        ttt2 = ttt*ttt
        ttt3 = ttt2*ttt

        meaneps = -46.8150 *ttt - 0.00059 *ttt2 + 0.001813 *ttt3 + 84381.448
        meaneps = math.fmod(meaneps/3600.0, 360.0)
        meaneps = meaneps * deg2rad

        l, l1, f, d, omega, lonmer, lonven, lonear, lonmar, lonjup, lonsat, \
            lonurn, lonnep, precrate = smu.fundarg(ttt, '80')
#fprintf(1, 'nut del arg %11.7f  %11.7f  %11.7f  %11.7f  %11.7f  \n', l*180/pi, l1*180/pi, f*180/pi, d*180/pi, omega*180/pi)

        deltapsi = 0.0
        deltaeps = 0.0
        #for i = 106:-1: 1

        for i in range (105, 0, -1):
            tempval = iar80[i, 0]*l + iar80[i, 1]*l1 + iar80[i, 2]*f + \
                     iar80[i, 3]*d + iar80[i, 4]*omega
            deltapsi = deltapsi + (rar80[i, 0]+rar80[i, 1]*ttt) * math.sin(tempval)
            deltaeps = deltaeps + (rar80[i, 2]+rar80[i, 3]*ttt) * math.cos(tempval)

        # --------------- find nutation parameters --------------------
        deltapsi = math.fmod(deltapsi + ddpsi, 2.0 * math.pi)
        deltaeps = math.fmod(deltaeps + ddeps, 2.0 * math.pi)
        trueeps = meaneps + deltaeps

#fprintf(1, 'meaneps %11.7f dp  %11.7f de  %11.7f te  %11.7f  ttt  %11.7f \n', meaneps*180/pi, deltapsi*180/pi, deltaeps*180/pi, trueeps*180/pi, ttt)

        cospsi = math.cos(deltapsi)
        sinpsi = math.sin(deltapsi)
        coseps = math.cos(meaneps)
        sineps = math.sin(meaneps)
        costrueeps = math.cos(trueeps)
        sintrueeps = math.sin(trueeps)

        nut = np.zeros((3, 3))
        nut[0, 0] = cospsi
        nut[0, 1] = costrueeps * sinpsi
        nut[0, 2] = sintrueeps * sinpsi
        nut[1, 0] = -coseps * sinpsi
        nut[1, 1] = costrueeps * coseps * cospsi + sintrueeps * sineps
        nut[1, 2] = sintrueeps * coseps * cospsi - sineps * costrueeps
        nut[2, 0] = -sineps * sinpsi
        nut[2, 1] = costrueeps * sineps * cospsi - sintrueeps * coseps
        nut[2, 2] = sintrueeps * sineps * cospsi + costrueeps * coseps

#         fprintf(1, 'nut matrix \n')
#         nut
#         fprintf(1, 'nut rotations \n')
#         n1 = smu.rot1mat(trueeps)
#         n2 = smu.rot3mat(deltapsi)
#         n3 = smu.rot1mat(-meaneps)
#         nut1 = n3*n2*n1

        return deltapsi, trueeps, meaneps, omega, nut


# ----------------------------------------------------------------------------
#
#                           function precess
#
#  this function calulates the transformation matrix that accounts for the effects
#    of precession. both the 1980 and 2006 theories are handled. note that the
#    required parameters differ a little.
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    ttt         - julian centuries of tt
#    opt         - method option                  '01', '02', '96', '80'
#
#  outputs       :
#    prec        - transformation matrix for mod - j2000 (80 only)
#    psia        - cannonical precession angle    rad    (00 only)
#    wa          - cannonical precession angle    rad    (00 only)
#    ea          - cannonical precession angle    rad    (00 only)
#    xa          - cannonical precession angle    rad    (00 only)
#
#  locals        :
#    ttt2        - ttt squared
#    ttt3        - ttt cubed
#    zeta        - precession angle               rad
#    z           - precession angle               rad
#    theta       - precession angle               rad
#    oblo        - obliquity value at j2000 epoch "%
#
#  coupling      :
#    none        -
#
#  references    :
#    vallado       2004, 214-216, 219-221
#
# [prec, psia, wa, ea, xa] = precess (ttt, opt)
# ----------------------------------------------------------------------------




def precess(ttt, opt):

    # " to rad
    convrt = math.pi / (180.0*3600.0)
    ttt2 = ttt * ttt
    ttt3 = ttt2 * ttt

    prec = np.eye(3)

    # ------------------- fk4 b1950 precession angles --------------------
    if (opt == '50'):

        '''
        % ---- Seidelmann pg 107
        % for these calls, ttt will come in with the current jd
        %            t1 = 0.0 %(ttt - 2433282.42345905)/365242.198782  % set start as B1850, 0.0 to simplify
        %            t2 = (ttt - 2396758.203)/365242.198782  % uses B1850
        %            ttt = t2 - t1
        %            ttt2 = ttt * ttt
        %            ttt3 = ttt * ttt2
        %            fprintf(1, '50prec %15.9f  \n', ttt)
        % exp supp 61 pg 38
        %            psia = 50.3708 + 0.0050 * ttt
        %            wa = 0.0
        %            ea = 0.0
        %            xa = 0.1247 - 0.0188 * ttt
        %            zeta = (23035.545 + 139.720*t1 + 0.060 *t1*t1)*ttt + (30.240-0.270*t1)*ttt2 + 17.995*ttt3 % "
        %            theta = (20051.12 - 85.29*t1 - 0.37 *t1*t1)*ttt + (-42.65-0.37*t1)*ttt2 - 41.80*ttt3
        %            z = (23035.545 + 139.720*t1 + 0.060 *t1*t1)*ttt + (109.480+0.390*t1)*ttt2 + 18.325*ttt3
        '''

        # ---- Newcomb Exp Supp 61 approach, but see GTDS pg 3-17
        # Exp Supp 61 says use 1900? but gtds says use 1950.
        # for these calls, ttt will come in with the current jd
        t1 = 0.0 #(ttt - 2415020.31352)/36524.2198782  # set start as B1900, 0.0 to simplify
        #            t2 = (ttt - 2415020.31352)/36525  # uses B1900
        t2 = (ttt - 2433282.42345905)/36525  # uses B1950
        #            ttt = t2 - t1
        #            ttt2 = ttt * ttt
        #            ttt3 = ttt * ttt2
        # exp supp 61 pg 38
        psia = 50.3708 + 0.0050 * ttt
        wa = 0.0 # not sure which one is which...
        ea = 84428.26 - 46.845*ttt - 0.00059*ttt2 + 0.00181*ttt3
        xa = 0.1247 - 0.0188 * ttt
        print('50prec %15.9f  \n' % ttt)
        # seems like Exp supp 61 is ok with 1900 as epoch, and Seidlemann is ok with listed measr,
        #            zeta = (2304.25 + 1.396*t1)*ttt + 0.302*ttt2 + 0.018*ttt3 # "
        #            theta = (2004.682 - 0.853*t1)*ttt -0.426*ttt2 - 0.042*ttt3
        #            z = (2304.25 + 1.396*t1)*ttt + 1.093*ttt2 + 0.018*ttt3
        # GTDS pg 3-17 using days from 1950 - avoids long rpecession
        # constants...
        zeta = 2304.9969*ttt + 0.302*ttt2 + 0.01808*ttt3 # "
        theta = 2004.2980*ttt -0.425936*ttt2 - 0.0416*ttt3
        z = 2304.9969*ttt + 1.092999*ttt2 + 0.0192*ttt3


        # tp-008 36-45
        # ttt is tropical centruies from 1950 36524.22 days
        prec = np.zeros((3, 3))
        prec[0, 0] = 1.0 - 2.9696e-4 * ttt2 - 1.3e-7 * ttt3
        prec[0, 1] = 2.234941e-2 * ttt + 6.76e-6 * ttt2 - 2.21e-6 * ttt3
        prec[0, 2] = 9.7169e-3 * ttt - 2.07e-6 * ttt2 - 9.6e-7 * ttt3
        prec[1, 0] = -prec[0, 1]
        prec[1, 1] = 1.0 - 2.4975e-4 * ttt2 - 1.5e-7 * ttt3
        prec[1, 2] = - 1.0858e-4 * ttt2
        prec[2, 0] = -prec[0, 2]
        prec[2, 1] = prec[1, 2]
        prec[2, 2] = 1.0 - 4.721e-5 * ttt2


        # pass these back out for testing
        psia = zeta
        wa = theta
        ea = z
        # ------------------- iau 76 precession angles --------------------
    else:
        if (opt == '80'):
            #     fprintf(1, '80prec %15.9f  \n', ttt)
            psia = 5038.7784*ttt - 1.07259*ttt2 - 0.001147*ttt3 # "
            wa = 84381.448 + 0.05127*ttt2 - 0.007726*ttt3
            ea = 84381.448 - 46.8150*ttt - 0.00059*ttt2 + 0.001813*ttt3
            xa = 10.5526*ttt - 2.38064*ttt2 - 0.001125*ttt3

            zeta = 2306.2181*ttt + 0.30188*ttt2 + 0.017998*ttt3 # "
            theta = 2004.3109*ttt - 0.42665*ttt2 - 0.041833*ttt3
            z = 2306.2181*ttt + 1.09468*ttt2 + 0.018203*ttt3
            # ------------------ iau 06 precession angles -------------------
        else:
            oblo = 84381.406 # "
            psia = ((((-0.0000000951 * ttt + 0.000132851) * ttt - 0.00114045)
                      * ttt - 1.0790069) * ttt + 5038.481507) * ttt # "
            wa = ((((0.0000003337 * ttt - 0.000000467) * ttt - 0.00772503)
                      * ttt + 0.0512623) * ttt -    0.025754) * ttt + oblo
            ea = ((((-0.0000000434 * ttt - 0.000000576) * ttt + 0.00200340)
                      * ttt - 0.0001831) * ttt -   46.836769) * ttt + oblo
            xa = ((((-0.0000000560 * ttt + 0.000170663) * ttt - 0.00121197)
                      * ttt - 2.3814292) * ttt +   10.556403) * ttt

            zeta = ((((-0.0000003173 * ttt - 0.000005971) * ttt + 0.01801828)
                      * ttt + 0.2988499) * ttt + 2306.083227) * ttt + 2.650545 # "
            theta = ((((-0.0000001274 * ttt - 0.000007089) * ttt - 0.04182264)
                      * ttt - 0.4294934) * ttt + 2004.191903) * ttt
            z = ((((0.0000002904 * ttt - 0.000028596) * ttt + 0.01826837)
                      * ttt + 1.0927348) * ttt + 2306.077181) * ttt - 2.650545

    # convert units to rad
    psia = psia  * convrt # rad
    wa = wa    * convrt
    ea = ea    * convrt
    xa = xa    * convrt

    zeta = zeta  * convrt
    theta = theta * convrt
    z = z     * convrt
    if (sh.iauhelp == 'y'):
        print('pr %11.7f  %11.7f  %11.7f %11.7fdeg \n'
              % (psia*180/math.pi, wa*180/math.pi, ea*180/math.pi,
                 xa*180/math.pi))
        print('pr %11.7f  %11.7f  %11.7fdeg \n'
              % (zeta*180/math.pi, theta*180/math.pi, z*180/math.pi))
    #print('pr %11.7f  %11.7f  %11.7fdeg \n' % (zeta*180/math.pi, theta*180/math.pi, z*180/math.pi))

    if (opt =='80') or (opt == '06'):
        coszeta = math.cos(zeta)
        sinzeta = math.sin(zeta)
        costheta = math.cos(theta)
        sintheta = math.sin(theta)
        cosz = math.cos(z)
        sinz = math.sin(z)

        # ----------------- form matrix  mod to j2000 -----------------
        prec = np.zeros((3, 3))
        prec[0, 0] = coszeta * costheta * cosz - sinzeta * sinz
        prec[0, 1] = coszeta * costheta * sinz + sinzeta * cosz
        prec[0, 2] = coszeta * sintheta
        prec[1, 0] = -sinzeta * costheta * cosz - coszeta * sinz
        prec[1, 1] = -sinzeta * costheta * sinz + coszeta * cosz
        prec[1, 2] = -sinzeta * sintheta
        prec[2, 0] = -sintheta * cosz
        prec[2, 1] = -sintheta * sinz
        prec[2, 2] = costheta

        '''
        % ----------------- do rotations instead ----------------------
        %             fprintf(1, 'prec matrix \n')
        %             prec
        %             fprintf(1, 'prec rotations z \n')
        %             p1 = rot3mat(z)
        %             p2 = rot2mat(-theta)
        %             p3 = rot3mat(zeta)
        %             prec1 = p3*p2*p1
        %
        %             fprintf(1, 'prec rotations w \n')
        %             a4 = rot3mat(-xa)
        %             a5 = rot1mat(wa)
        %             a6 = rot3mat(psia)
        %             a7 = rot1mat(-ea)
        %             prec2 = a7*a6*a5*a4
    else  #if %(strcmp(opt, '50') ~= 1)
        %             oblo = oblo * convrt % " to rad
        %             a4 = rot3mat(-xa)
        %             a5 = rot1mat(wa)
        %             a6 = rot3mat(psia)
        %             a7 = rot1mat(-oblo)
        %             prec3 = a7*a6*a5*a4
        '''
    return prec, psia, wa, ea, xa


# ------------------------------------------------------------------------------
#
#                           function sight
#
#  this function takes the position vectors of two satellites and determines
#    if there is line-of-sight between the two satellites.  an oblate earth
#    with radius of 1 er is assumed.  the process forms the equation of
#    a line between the two vectors.  differentiating and setting to zero finds
#    the minimum value, and when plugged back into the original line equation,
#    gives the minimum distance.  the parameter tmin is allowed to range from
#    0.0  to 1.0 .  scale the k-component to account for oblate earth because it's
#    the only qunatity that changes.
#
#  author        : david vallado                  719-573-2600   31 oct 2003
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    r1          - position vector of the 1st sat km
#    r2          - position vector of the 2nd sat km
#    whichkind   - spherical or ellipsoidal earth 's', 'e'*default
#
#  outputs       :
#    los         - line of sight                  'yes', 'no ' ##True or False now
#
#  locals        :
#    tr1         - scaled r1 vector
#    tr2         - scaled r2 vector
#    adotb       - dot product of a dot b
#    tmin        - minimum value of t from a to b
#    distsqrd    - min distance squared to earth
#    asqrd       - magnitude of a squared
#    bsqrd       - magnitude of b squared
#
#  coupling:
#
#  references    :
#    vallado       2001, 291-295, alg 35, ex 5-3
#
# [los] = sight (r1, r2, whichkind)
# ------------------------------------------------------------------------------


def sight (r1, r2, whichkind):
    # -------------------------  implementation   -----------------
    tr1 = np.array([0.0, 0.0, 0.0])
    tr2 = np.array([0.0, 0.0, 0.0])
    for i in range(3):
        tr1[i] = r1[i]
        tr2[i] = r2[i]
    magr1 = smu.mag(tr1)
    magr2 = smu.mag(tr2)

    # --------------------- scale z component ---------------------
    if (whichkind == 'e'):
        temp = 1.0 /math.sqrt(1.0 - eccearthsqrd)
    else:
        temp = 1.0
    tr1[2] = tr1[2]*temp
    tr2[2] = tr2[2]*temp
    bsqrd = magr2*magr2
    asqrd = magr1*magr1
    adotb = np.dot(tr1, tr2)

    # ---------------------- find tmin ----------------------------
    distsqrd = 0.0
    if (abs(asqrd + bsqrd - 2.0 *adotb) < 0.0001):
        tmin = 0.0
    else:
        tmin = (asqrd - adotb) / (asqrd + bsqrd - 2.0 *adotb)
    # ----------------------- check los ---------------------------
    if ((tmin < 0.0) or (tmin > 1.0)):
        los = True
    else:
        distsqrd = ((1.0 -tmin)*asqrd + adotb*tmin)/re**2
        if (distsqrd > 1.0):
            los = True
        else:
            los = False
    return los


# ------------------------------------------------------------------------------
#
#                           function sun
#
#  this function calculates the geocentric equatorial position vector
#    the sun given the julian date. Sergey K (2022) has noted that improved results
#    are found assuming the output is in a precessing frame (TEME) and converting to ICRF.
#    this is the low precision formula and is valid for years from 1950 to 2050.
#    accuaracy of apparent coordinates is about 0.01 degrees. notice many of
#    the calculations are performed in degrees, and are not changed until later.
#    this is due to the fact that the almanac uses degrees exclusively in their formulations.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#    vallado     - fix mean lon of sun                            7 may 2004
#
#  inputs          description                            range / units
#    jd          - julian date (UTC)                         days from 4713 bc
#
#  outputs       :
#    rsun        - inertial position vector of the sun       au
#    rtasc       - right ascension                           rad
#    decl        - declination                               rad
#
#  locals        :
#    meanlong    - mean longitude
#    meananomaly - mean anomaly
#    eclplong    - ecliptic longitude
#    obliquity   - mean obliquity of the ecliptic
#    tut1        - julian centuries of ut1 from
#                  jan 1, 2000 12h
#    ttdb        - julian centuries of tdb from
#                  jan 1, 2000 12h
#    hr          - hours                                   0 .. 24              10
#    min         - minutes                                 0 .. 59              15
#    sec         - seconds                                 0.0  .. 59.99          30.00
#    temp        - temporary variable
#    deg         - degrees
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 281, alg 29, ex 5-1
#
# [rsun, rtasc, decl] = sun (jd)
# ------------------------------------------------------------------------------


def sun(jd):

    # -------------------------  implementation   -----------------
    # -------------------  initialize values   --------------------
    tut1 = (jd - 2451545.0) / 36525.0

    if sh.show == 'y':
        print('tut1 %14.9f \n' % tut1)

    meanlong = 280.460 + 36000.77 * tut1
    meanlong = math.fmod(meanlong, 360.0)  #deg

    ttdb = tut1
    meananomaly = 357.5277233 + 35999.05034 *ttdb
    meananomaly = math.fmod(meananomaly*deg2rad , twopi)  #rad
    if (meananomaly < 0.0):
        meananomaly = twopi + meananomaly

    eclplong = meanlong + 1.914666471 * math.sin(meananomaly) \
                + 0.019994643 * math.sin(2.0 * meananomaly) #deg
    eclplong = math.fmod(eclplong, 360.0)  #deg

    obliquity = 23.439291 - 0.0130042 * ttdb  #deg

    eclplong = eclplong * deg2rad
    obliquity = obliquity * deg2rad

    # --------- find magnitude of sun vector, )   components ------
    magr = 1.000140612 - 0.016708617 * math.cos(meananomaly) \
        - 0.000139589 * math.cos(2.0 * meananomaly)    # in au's


    rsun = np.array([0.0, 0.0, 0.0])
    rsun[0] = magr*math.cos(eclplong)
    rsun[1] = magr*math.cos(obliquity)*math.sin(eclplong)
    rsun[2] = magr*math.sin(obliquity)*math.sin(eclplong)

    if sh.show == 'y':
        print('meanlon %11.6f meanan %11.6f eclplon %11.6f obli %11.6f \n'
              % (meanlong, meananomaly/deg2rad, eclplong/deg2rad,
                 obliquity/deg2rad))
        print('rs %11.9f %11.9f %11.9f \n' % (rsun[0], rsun[1], rsun[2]))
        print('magr %14.7f \n' % magr)

    rtasc = math.atan(math.cos(obliquity) * math.tan(eclplong))

    # --- check that rtasc is in the same quadrant as eclplong ----
    if (eclplong < 0.0):
        eclplong = eclplong + twopi    # make sure it's in 0 to 2pi range
    if (abs(eclplong-rtasc) > math.pi*0.5):
        rtasc = rtasc + halfpi*round((eclplong-rtasc)/halfpi)
    decl = math.asin(math.sin(obliquity)*math.sin(eclplong))

    return rsun, rtasc, decl






# ---------------------------------------------------------------------------
#
#                           function site
#
#  this function finds the position and velocity vectors for a site.  the
#    answer is returned in the geocentric equatorial (ecef) coordinate system.
#    note that the velocity is zero because the coordinate system is fixed to
#    the earth.
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#    vallado     - fix velocity vector                           23 jul 2002
#
#  inputs          description                    range / units
#    latgd       - geodetic latitude              -pi/2 to pi/2 rad
#    lon         - longitude of site              -2pi to 2pi rad
#    alt         - altitude                       km
#
#  outputs       :
#    rs          - ecef site position vector      km
#    vs          - ecef site velocity vector      km/s
#
#  locals        :
#    sinlat      - variable containing  sin(lat)  rad
#    temp        - temporary real value
#    rdel        - rdel component of site vector  km
#    rk          - rk component of site vector    km
#    cearth      -
#
#  coupling      :
#    none
#
#  references    :
#    vallado       2001, 404-407, alg 47, ex 7-1
#
# [rs, vs] = site (latgd, lon, alt)
# -----------------------------------------------------------------------------


def site (latgd, lon, alt):

        # EGM-08 constants used here
        flat = 1.0/298.257223563
        earthrot = 7.292115e-5     # rad/s  old 7.29211514670698e-05
        mu = 398600.4415      # km3/s2
        mum = 3.986004415e14   # m3/s2
        # derived constants from the base values
        eccearth = math.sqrt(2.0*flat - flat**2)
        eccearthsqrd = eccearth**2
        # -------------------------  implementation   -----------------
        sinlat = np.sin((latgd))

        # ------  find rdel and rk components of site vector  ---------
        cearth = re / math.sqrt(1.0 - (eccearthsqrd*sinlat*sinlat))
        rdel = (cearth + alt)*np.cos((latgd))
        rk = ((1.0-eccearthsqrd)*cearth + alt)*sinlat

        # ---------------  find site position vector  -----------------
        rs = np.zeros((3))
        rs[0] = rdel * np.cos((lon))
        rs[1] = rdel * np.sin((lon))
        rs[2] = rk
        rs = rs.T

        # ---------------  find site velocity vector  -----------------
        #ome = [0.0 0.0 omegaearth]
        #[vs] = cross(ome, rs)
        vs = np.zeros((3))

        return rs, vs



# ------------------------------------------------------------------------------
#
#                           function sunalmanac
#
#  this function calculates the geocentric equatorial position vector
#    the sun given the julian date.  this is the low precision formula and
#    is valid for years from 1950 to 2050.  accuaracy of apparent coordinates
#    is 0.01  degrees.  notice many of the calculations are performed in
#    degrees, and are not changed until later.  this is due to the fact that
#    the almanac uses degrees exclusively in their formulations.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#    vallado     - fix mean lon of sun                            7 mat 2004
#
#  inputs          description                    range / units
#    jd          - julian date                    days from 4713 bc
#
#  outputs       :
#    rsun        - ijk position vector of the sun au
#    rtasc       - right ascension                rad
#    decl        - declination                    rad
#
#  locals        :
#    meanlong    - mean longitude
#    meananomaly - mean anomaly
#    eclplong    - ecliptic longitude
#    obliquity   - mean obliquity of the ecliptic
#    tut1        - julian centuries of ut1 from
#                  jan 1, 2000 12h
#    ttdb        - julian centuries of tdb from
#                  jan 1, 2000 12h
#    hr          - hours                          0 .. 24              10
#    min         - minutes                        0 .. 59              15
#    sec         - seconds                        0.0  .. 59.99          30.00
#    temp        - temporary variable
#    deg         - degrees
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 281, alg 29, ex 5-1
#
# [rsun, rtasc, decl] = sunalmanac (jd)
# ------------------------------------------------------------------------------


def sunalmanac(jd):

    # -------------------------  implementation   -----------------
    # -------------------  initialize values   --------------------
    tut1 = (jd - 2451545.0) / 36525.0
    print('tut1 %14.9f \n' % tut1)

    meanlong = 280.460 + 36000.771285*tut1
    meanlong = np.fmod(meanlong, 360.0)  #deg

    ttdb = tut1
    meananomaly = 357.528 + 35999.050957 * ttdb
    meananomaly = np.fmod(meananomaly*deg2rad, twopi)  #rad
    if (meananomaly < 0.0):
        meananomaly = twopi + meananomaly

    eclplong = meanlong + 1.915 * math.sin(meananomaly) \
                + 0.020 * math.sin(2.0 *meananomaly) #deg
    eclplong = np.fmod(eclplong, 360.0)  #deg

    obliquity = 23.439 - 0.01461 * ttdb  #deg

    eclplong = eclplong * deg2rad
    obliquity = obliquity * deg2rad

    # --------- find magnitude of sun vector, )   components ------
    magr = 1.00014  - 0.01671 * math.cos(meananomaly) \
        - 0.00014 *math.cos(2.0 * meananomaly)    # in au's

    rsun = np.zeros((3))
    rsun[0] = magr*math.cos(eclplong)
    rsun[1] = magr*math.cos(obliquity) * math.sin(eclplong)
    rsun[2] = magr*math.sin(obliquity) * math.sin(eclplong)

    print('meanlon %11.6f meanan %11.6f eclplon %11.6f obli %11.6f \n'
          % (meanlong, meananomaly/deg2rad,
             eclplong/deg2rad, obliquity/deg2rad))
    print('rs %11.9f %11.9f %11.9f \n' % (rsun[0], rsun[1], rsun[2]))
    print('magr %14.7f \n' % magr)

    rtasc = math.atan(math.cos(obliquity)*math.tan(eclplong))

    # --- check that rtasc is in the same quadrant as eclplong ----
    if (eclplong < 0.0):
        eclplong = eclplong + twopi    # make sure it's in 0 to 2pi range
    if (abs(eclplong-rtasc) > math.pi*0.5):
        rtasc = rtasc + 0.5 * math.pi*round((eclplong-rtasc)/(0.5 * math.pi))
    decl = math.asin(math.sin(obliquity) * math.sin(eclplong))

    return rsun, rtasc, decl



# ------------------------------------------------------------------------------
#
#                           function light
#
#  this function determines if a spacecraft is sunlit or in the dark at a
#    particular time.  an oblate earth and cylindrical shadow is assumed.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    r           - position vector of sat         er
#    jd          - julian date at desired time    days from 4713 bc
#    whichkind   - spherical or ellipsoidal earth 's', 'e'*default
#
#  outputs       :
#    vis         - visibility flag                'yes', 'no '
#
#  locals        :
#    rtasc       - suns right ascension           rad
#    decl        - suns declination               rad
#    rsun        - sun vector                     au
#    auer        - conversion from au to er
#
#  coupling      :
#    sun         - position vector of sun
#    lncom1      - multiple a vector by a constant
#    sight       - does line-of-sight exist beteen vectors
#
#  references    :
#    vallado       2001, 291-295, alg 35, ex 5-6
#
# [lit] = light (r, jd, whichkind)
# ------------------------------------------------------------------------------


def light (r, jd, whichkind):

        # -------------------------  implementation   -------------------------
        rsun, rtasc, decl = sun(jd)
        rsun = auer*rsun

        # ------------ is the satellite in the shadow? ----------------
        lit = sight(rsun, r, whichkind)
        return lit


# -----------------------------------------------------------------------------
#
#                           function sunriset
#
#  this function finds the universal time for sunrise and sunset given the
#    day and site location. use (- lon) for local time
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    jd          - julian date                    days from 4713 bc
#    latgd       - site latitude (south -)        -65 to 65 rad
#    lon         - site longitude (west -)        -2pi to 2pi rad
#    whichkind   - character for which rise/set   's' 'c' 'n' 'a'
#
#  outputs       :
#    utsunrise   - universal time of sunrise      hrs
#    utsunset    - universal time of sunset       hrs
#    error       - error parameter
#
#  locals        :
#    sunangle    - angle between the sun vector
#                  and a point on the earth     rad
#    jdtemp      - julian date for sunrise/set    days from 4713 bc
#    uttemp      - temporary ut time              days
#    tut1        - julian centuries from the
#                  jan 1, 2000 12 h epoch (ut1)
#    ra          - right ascension                rad
#    decl        - declination                    rad
#    meanlonsun  -                                rad
#    meananomalysun                               rad
#    lonecliptic - longitude of the ecliptic      rad
#    obliquity   - obliquity of the ecliptic      rad
#    gst         - for 0 h utc of each day        rad
#    lha         - local hour angle               rad
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - hour                           0 .. 23
#    min         - minute                         0 .. 59
#    sec         - second                         0.0d0 .. 59.999d0
#    opt         - idx to for rise and set calc    1, 2
#
#  coupling      :
#    invjday   - finds the year day mon hr min sec from the julian date
#    jday   - finds the julian date given year, mon day, hr, min, sec
#
#  references    :
#    vallado       2007, 283, Alg 30, Ex 5-2
#
# [utsunrise, utsunset, error] = sunriset(jd, latgd, lon, whichkind)
# -----------------------------------------------------------------------------

def sunriset(jd=None, latgd=None, lon=None, whichkind=None):
    # ------------------------  implementation   ------------------
    # -------------- make sure lon is within +- 180 deg -----------
    if (lon > np.pi):
        lon = lon - twopi

    if (lon < - np.pi):
        lon = lon + twopi

    if (whichkind == 's'):
        sunangle = (90.0 + 50.0 / 60.0) * deg2rad

    if (whichkind == 'c'):
        sunangle = 96.0 * deg2rad

    if (whichkind == 'n'):
        sunangle = 102.0 * deg2rad

    if (whichkind == 'a'):
        sunangle = 108.0 * deg2rad

    year, month, day, hr, min, sec = stu.invjday(jd, 0.0)
    for opt in range(3):
        error = 'ok'
        if (opt == 1):
            jdtemp, jdtempf = stu.jday(year, month, day, 6, 0, 0.0)
        else:
            jdtemp, jdtempf = stu.jday(year, month, day, 18, 0, 0.0)
        jdtemp = jdtemp + jdtempf - lon * rad2deg / 15.0 / 24.0
        jdtemp
        lon * rad2deg / 15.0 / 24.0
        tut1 = (jdtemp - 2451545.0) / 36525.0
        meanlonsun = 280.4606184 + 36000.77005361 * tut1
        #            meanlonsun = 280.460 + 36000.770 * tut1
        meananomalysun = 357.5277233 + 35999.05034 * tut1
        meananomalysun = np.fmod(meananomalysun * deg2rad, twopi)
        if (meananomalysun < 0.0):
            meananomalysun = meananomalysun + twopi
        lonecliptic = (meanlonsun + 1.914666471 * np.sin(meananomalysun)
                       + 0.019994643 * np.sin(2.0 * meananomalysun))
        lonecliptic = np.fmod(lonecliptic * deg2rad, twopi)
        if (lonecliptic < 0.0):
            lonecliptic = lonecliptic + twopi
        obliquity = 23.439291 - 0.0130042 * tut1
        obliquity = obliquity * deg2rad
        print('lonecl %11.7f tut1 %11.7f  obl %11.7f  \n'
              % (lonecliptic * rad2deg, tut1, obliquity * rad2deg))
        ra = np.arctan(np.cos(obliquity) * np.tan(lonecliptic))
        decl = np.arcsin(np.sin(obliquity) * np.sin(lonecliptic))
        if (ra < 0.0):
            ra = ra + twopi
        if (((lonecliptic > np.pi) and (ra < np.pi))):
            ra = ra + np.pi
        if (((lonecliptic < np.pi) and (ra > np.pi))):
            ra = ra - np.pi
        print('mlonsun %11.7f meanansun %11.7f  eclon %11.7f  \n'
              % (meanlonsun + 1080.0, meananomalysun * rad2deg,
                 lonecliptic * rad2deg))
        print('ra %11.7f decl %11.7f  \n' % (ra * rad2deg, decl * rad2deg))
        lha = ((np.cos(sunangle) - np.sin(decl) * np.sin(latgd))
               / (np.cos(decl) * np.cos(latgd)))
        #fprintf(1, 'lha 1st #11.7f   \n', 90.0 + lha*rad2deg)
        if (np.abs(lha) <= 1.0):
            lha = np.arccos(lha)
        else:
            error = 'not ok'
            print('error \n' % ())
        print('lha 2nd  %11.7f   \n' % (lha * rad2deg))
        if (error == 'ok'):
            if (opt == 1):
                lha = twopi - lha
            gst = (1.75336855923327 + 628.331970688841 * tut1
                   + 6.77071394490334e-06 * tut1 * tut1
                   - 4.50876723431868e-10 * tut1 * tut1 * tut1)
            gst = np.fmod(gst, twopi)
            if (gst < 0.0):
                gst = gst + twopi
            print('lha %11.7f gst %11.7f  \n' % (lha * rad2deg, (gst - twopi) * rad2deg))
            uttemp = lha + ra - gst
            print('gst %11.7f uttemp %11.7f uttemp %11.7f   \n'
                  % (gst * rad2deg, uttemp * rad2deg,
                     (twopi + uttemp) * rad2deg))
            uttemp = uttemp * rad2deg / 15.0
            uttemp = np.fmod(uttemp, 24.0)
            if (uttemp < 0.0):
                uttemp = uttemp + 24.0
                error = 'day before'
            if (uttemp > 24.0):
                uttemp = uttemp - 24.0
                error = 'day after'
        else:
            uttemp = 99.99
        if (opt == 1):
            utsunrise = uttemp
        else:
            utsunset = uttemp

    return utsunrise, utsunset, error


##############################################################################################################
# ------------------------------------------------------------------------------
#
#                           function moon
#
#  this function calculates the geocentric equatorial (ijk) position vector
#    for the moon given the julian date.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    jd          - julian date                    days from 4713 bc
#
#  outputs       :
#    rmoon       - ijk position vector of moon    er
#    rtasc       - right ascension                rad
#    decl        - declination                    rad
#
#  locals        :
#    eclplong    - ecliptic longitude
#    eclplat     - eclpitic latitude
#    hzparal     - horizontal parallax
#    l           - geocentric direction cosines
#    m           -             "     "
#    n           -             "     "
#    ttdb        - julian centuries of tdb from
#                  jan 1, 2000 12h
#    hr          - hours                          0 .. 24
#    min         - minutes                        0 .. 59
#    sec         - seconds                        0.0  .. 59.99
#    deg         - degrees
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 290, alg 31, ex 5-3
#
# [rmoon, rtasc, decl] = moon (jd, show)
# ------------------------------------------------------------------------------

def moon(jd=None, show=None):

    # -------------------------  implementation   -----------------
    ttdb = (jd - 2451545.0) / 36525.0
    eclplong = (218.32 + 481267.8813 * ttdb
                + 6.29 * np.sin((134.9 + 477198.85 * ttdb) * deg2rad)
                - 1.27 * np.sin((259.2 - 413335.38 * ttdb) * deg2rad)
                + 0.66 * np.sin((235.7 + 890534.23 * ttdb) * deg2rad)
                + 0.21 * np.sin((269.9 + 954397.7 * ttdb) * deg2rad)
                - 0.19 * np.sin((357.5 + 35999.05 * ttdb) * deg2rad)
                - 0.11 * np.sin((186.6 + 966404.05 * ttdb) * deg2rad))

    eclplat = (5.13 * np.sin((93.3 + 483202.03 * ttdb) * deg2rad)
               + 0.28 * np.sin((228.2 + 960400.87 * ttdb) * deg2rad)
               - 0.28 * np.sin((318.3 + 6003.18 * ttdb) * deg2rad)
               - 0.17 * np.sin((217.6 - 407332.2 * ttdb) * deg2rad))

    hzparal = (0.9508 + 0.0518 * np.cos((134.9 + 477198.85 * ttdb) * deg2rad)
               + 0.0095 * np.cos((259.2 - 413335.38 * ttdb) * deg2rad)
               + 0.0078 * np.cos((235.7 + 890534.23 * ttdb) * deg2rad)
               + 0.0028 * np.cos((269.9 + 954397.7 * ttdb) * deg2rad))

    eclplong = np.fmod(eclplong * deg2rad, twopi)
    eclplat = np.fmod(eclplat * deg2rad, twopi)
    hzparal = np.fmod(hzparal * deg2rad, twopi)
    obliquity = 23.439291 - 0.0130042 * ttdb

    obliquity = obliquity * deg2rad
    if show == 'y':
        360 + eclplong / deg2rad
        eclplat / deg2rad
        hzparal / deg2rad
        obliquity / deg2rad

    # ------------ find the geocentric direction cosines ----------
    l = np.cos(eclplat) * np.cos(eclplong)
    m = (np.cos(obliquity) * np.cos(eclplat) * np.sin(eclplong)
         - np.sin(obliquity) * np.sin(eclplat))
    n = (np.sin(obliquity) * np.cos(eclplat) * np.sin(eclplong)
         + np.cos(obliquity) * np.sin(eclplat))
    # ------------- calculate moon position vector ----------------
    magr = 1.0 / np.sin(hzparal)
    if show == 'y':
        magr * re

    rmoon = np.zeros(3)
    rmoon[0] = magr * l
    rmoon[1] = magr * m
    rmoon[2] = magr * n
    # -------------- find rt ascension and declination ------------
    rtasc = math.atan2(m, l)
    decl = np.arcsin(n)
    return rmoon, rtasc, decl


# -----------------------------------------------------------------------------
#
#                           function moonriset
#
#  this function finds the universal time for moonrise and moonset given the
#    day and site location.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  revisions
#                - david vallado                                 25 jan 2011
#
#  inputs          description                    range / units
#    jd          - julian date                    days from 4713 bc
#    latgd       - site latitude (south -)        -65 to 65 rad
#    lon         - site longitude (west -)        -2pi to 2pi rad
#
#  outputs       :
#    utmoonrise  - universal time of moonrise     hrs
#    utmoonset   - universal time of moonset      hrs
#    moonphaseang- phase angle of the moon        deg
#    error       - error parameter
#
#  locals        :
#    moonangle   - angle between the moon vector
#                  and a point on the earth       rad
#    jdtemp      - julian date for moonrise/set   days from 4713 bc
#    uttemp      - temporary ut time              days
#    tut1        - julian centuries from the
#                  jan 1, 2000 12 h epoch (ut1)
#    rtasc       - right ascension                rad
#    decl        - declination                    rad
#    meanlonmoon -                                rad
#    meananomaly -                                rad
#    eclplong    - longitude of the ecliptic      rad
#    obliquity   - obliquity of the ecliptic      rad
#    rmoon
#    rmoonrs
#    rv
#    rhosat
#    tsry1
#    l, m, n     - direction cosines
#    eclplat
#    moongha, moonghan
#    dgha, dghan
#    lhan
#    lst
#    deltaut, deltautn
#    t, tn
#    hzparal
#    loneclsun
#    loneclmoon
#    ttdb
#    gst         - for 0 h utc of each day        rad
#    lha         - local hour angle               rad
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - hour                           0 .. 23
#    min         - minute                         0 .. 59
#    sec         - second                         0.0  .. 59.999
#    opt         - idx to for rise and set calc    1, 2
#
#  coupling      :
#    invjulian- finds the year day mon hr min sec from the julian date
#    julianday   - finds the julian date given year, mon day, hr, min, sec
#
#  references    :
#    vallado       2007, 292, Alg 32, Ex 5-4
#
# [utmoonrise, utmoonset, moonphaseang, error] = moonrise(jd, latgd, lon)
# -----------------------------------------------------------------------------

def moonrise(jd=None, latgd=None, lon=None, show=None):
    # ------------------------  implementation   ------------------
    error = 'ok'
    # -------------- for once for moonrise (1), ) set (2) ---------
# -------------- make sure lon is within +- 180 deg -----------
    if (lon > np.pi):
        lon = lon - 2.0 * np.pi

    if (lon < - np.pi):
        lon = lon + twopi

    try1 = 1

    opt = 1

    while (opt <= 2):

        year, month, day, hr, min, sec = stu.invjday(jd, 0.0)
        jdtemp, jdtempf = stu.jday(year, month, day, 0, 0, 0.0)
        if (opt == 1):
            uttemp = (6.0 - lon / 15.0) / 24.0
        else:
            uttemp = (18.0 + lon / 15.0) / 24.0
        if (try1 == 2):
            uttemp = 0.5
        i = 0
        tn = uttemp
        t = tn + 10.0
        jdtemp = jdtemp + jdtempf + uttemp
        while ((np.abs(tn - t) >= 0.008) and (i <= 5)):

            ttdb = (jdtemp + jdtempf - 2451545.0) / 36525.0
            eclplong = (218.32 + 481267.8813 * ttdb
                        + 6.29 * np.sin((134.9 + 477198.85 * ttdb) * deg2rad)
                        - 1.27 * np.sin((259.2 - 413335.38 * ttdb) * deg2rad)
                        + 0.66 * np.sin((235.7 + 890534.23 * ttdb) * deg2rad)
                        + 0.21 * np.sin((269.9 + 954397.7 * ttdb) * deg2rad)
                        - 0.19 * np.sin((357.5 + 35999.05 * ttdb) * deg2rad)
                        - 0.11 * np.sin((186.6 + 966404.05 * ttdb) * deg2rad))
            eclplat = (5.13 * np.sin((93.3 + 483202.03 * ttdb) * deg2rad)
                       + 0.28 * np.sin((228.2 + 960400.87 * ttdb) * deg2rad)
                       - 0.28 * np.sin((318.3 + 6003.18 * ttdb) * deg2rad)
                       - 0.17 * np.sin((217.6 - 407332.2 * ttdb) * deg2rad))
            eclplong = np.fmod(eclplong * deg2rad, twopi)
            eclplat = np.fmod(eclplat * deg2rad, twopi)
            if show == 'y':
                print('%2d %2d ecpllon %11.7f ecllat %11.7f '
                      % (opt, try1, eclplong / deg2rad, eclplat / deg2rad))
            obliquity = 23.439291 - 0.0130042 * ttdb
            obliquity = obliquity * deg2rad
            # ------- find the geocentric direction cosines -------
            l = np.cos(eclplat) * np.cos(eclplong)
            m = (np.cos(obliquity) * np.cos(eclplat) * np.sin(eclplong)
                 - np.sin(obliquity) * np.sin(eclplat))
            n = (np.sin(obliquity) * np.cos(eclplat) * np.sin(eclplong)
                 + np.cos(obliquity) * np.sin(eclplat))
            if show == 'y':
                print('l %11.7f m %11.7f n %11.7f ' % (l, m, n))
            rtasc = math.atan2(m, l)
            # - check that rtasc is in the same quadrant as eclplong
            if (eclplong < 0.0):
                eclplong = eclplong + twopi
            if (np.abs(eclplong - rtasc) > np.pi * 0.5):
                rtasc = (rtasc + 0.5 * np.pi
                         * np.rint(0.5 + (eclplong - rtasc) / (0.5 * np.pi)))
            decl = np.arcsin(n)
            lst, gst = stu.lstime(lon, jdtemp + jdtempf)
            if show == 'y':
                print('ra %8.5f dcl %8.5f lst %8.5f jdtemp %8.5f \n'
                      % (rtasc / deg2rad, decl / deg2rad, lst / deg2rad,
                         jdtemp + jdtempf))
            moonghan = lst - lon - rtasc
            #cdav
            if (i == 0):
                lha = moonghan + lon
                dgha = 347.81 * deg2rad
            else:
                dgha = (moonghan - moongha) / deltaut
            if (dgha < 0.0):
                dgha = dgha + twopi / np.abs(deltaut)
            if show == 'y':
                print('mn gha %11.7f  dgha  %11.7f '
                      % (moonghan / deg2rad, dgha / deg2rad))
            lhan = ((0.00233 - np.sin(latgd) * np.sin(decl))
                    / (np.cos(latgd) * np.cos(decl)))
            if show == 'y':
                print('lhan  %11.7f rad ' % (lhan))
            #fprintf('lhan  #11.7f deg \n', lhan/deg2rad)
            if (lhan > 1.0):
                lhan = 0.0
            if (lhan < - 1.0):
                lhan = - 1.0
            lhan = np.arccos(lhan)
            if (opt == 1):
                lhan = twopi - lhan
            if show == 'y':
                print('lhan1 %11.7f ' % (lhan / deg2rad))
            if (np.abs(dgha) > 0.0001):
                deltaut = (lhan - lha) / dgha
            else:
                deltaut = 1.0
                error = 'error1 dgha is too small'
            if show == 'y':
                print('deltaut %11.7f tn  %11.7f ' % (deltaut, tn))
            t = tn
            if (np.abs(deltaut) > 0.5):
                if (np.abs(dgha) > 0.001):
                    if (deltaut < 0.0):
                        deltaut = deltaut + twopi / dgha
                        if (np.abs(deltaut) > 0.51):
                            i = 6
                    else:
                        deltaut = deltaut - twopi / dgha
                        if (np.abs(deltaut) > 0.51):
                            i = 6
                else:
                    error = 'error2 dgha is too small'
            tn = uttemp + deltaut
            jdtemp = jdtemp + jdtempf - uttemp + tn
            i = i + 1
            moongha = moonghan
            if show == 'y':
                print('deltaut %11.7f  jd  %14.4f  tn  %11.7f \n'
                      % (deltaut, jdtemp + jdtempf, tn))

        uttemp = tn * 24.0
        if (i > 5):
            uttemp = 9999.99
        if (uttemp > 24.0) and (uttemp < 9999.0):
            #fprintf('rem #11.7f ', uttemp)
            uttemp = np.fmod(uttemp, 24.0)
            #fprintf('#11.7f ', uttemp)
        if (uttemp < 0.0):
            uttemp = uttemp + 24.0
        #    if (uttemp > 900)
#        uttemp = 24.0
#    end
        if (opt == 1):
            utmoonrise = uttemp
        if (opt == 2):
            utmoonset = uttemp
        # update the iteration and check for solution
        try1 = try1 + 1
        if ((i > 5) and (try1 < 3)):
            if show == 'y':
                print('try1 #2 %4d' % (opt))
        else:
            if ((i > 5) and (try1 > 2)):
                if (opt == 1):
                    error = 'no rise'
                    uttemp = 9999.99
                if (opt == 2):
                    error = 'no set'
                    uttemp = 9997.99
            opt = opt + 1
            try1 = 1


    # ------------- determine phase angle of the moon -------------
    meanlong = 280.4606184 + 36000.77005361 * ttdb
    meanlong = np.fmod(meanlong, 360.0)
    meananomaly = 357.5277233 + 35999.05034 * ttdb
    meananomaly = np.fmod(meananomaly * deg2rad, twopi)
    if (meananomaly < 0.0):
        meananomaly = twopi + meananomaly

    loneclsun = (meanlong + 1.914666471 * np.sin(meananomaly)
                 + 0.019994643 * np.sin(2.0 * meananomaly))
    loneclmoon = (218.32 + 481267.8813 * ttdb
                  + 6.29 * np.sin((134.9 + 477198.85 * ttdb) * deg2rad)
                  - 1.27 * np.sin((259.2 - 413335.38 * ttdb) * deg2rad)
                  + 0.66 * np.sin((235.7 + 890534.23 * ttdb) * deg2rad)
                  + 0.21 * np.sin((269.9 + 954397.7 * ttdb) * deg2rad)
                  - 0.19 * np.sin((357.5 + 35999.05 * ttdb) * deg2rad)
                  - 0.11 * np.sin((186.6 + 966404.05 * ttdb) * deg2rad))
    loneclmoon = np.fmod(loneclmoon, 360.0)
    moonphaseang = loneclmoon - loneclsun
    if (moonphaseang < 0.0):
        moonphaseang = 360.0 + moonphaseang

    return utmoonrise, utmoonset, moonphaseang, error

# -----------------------------------------------------------------------------
#
#                           function moonriset
#
#  this function finds the universal time for moonrise and moonset given the
#    day and site location. notice this routine kicks out if the event
#    occurs on the next or previous day.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  revisions
#                - david vallado                                 25 jan 2011
#
#  inputs          description                    range / units
#    jd          - julian date                    days from 4713 bc
#    latgd       - site latitude (south -)        -65 to 65 rad
#    lon         - site longitude (west -)        -2pi to 2pi rad
#
#  outputs       :
#    utmoonrise  - universal time of moonrise     hrs
#    utmoonset   - universal time of moonset      hrs
#    moonphaseang- phase angle of the moon        deg
#    error       - error parameter
#
#  locals        :
#    moonangle   - angle between the moon vector
#                  and a point on the earth       rad
#    jdtemp      - julian date for moonrise/set   days from 4713 bc
#    uttemp      - temporary ut time              days
#    tut1        - julian centuries from the
#                  jan 1, 2000 12 h epoch (ut1)
#    rtasc       - right ascension                rad
#    decl        - declination                    rad
#    meanlonmoon -                                rad
#    meananomaly -                                rad
#    eclplong    - longitude of the ecliptic      rad
#    obliquity   - obliquity of the ecliptic      rad
#    rmoon
#    rmoonrs
#    rv
#    rhosat
#    tsry1
#    l, m, n     - direction cosines
#    eclplat
#    moongha, moonghan
#    dgha, dghan
#    lhan
#    lst
#    deltaut, deltautn
#    t, tn
#    hzparal
#    loneclsun
#    loneclmoon
#    ttdb
#    gst         - for 0 h utc of each day        rad
#    lha         - local hour angle               rad
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - hour                           0 .. 23
#    min         - minute                         0 .. 59
#    sec         - second                         0.0  .. 59.999
#    opt         - idx to for rise and set calc    1, 2
#
#  coupling      :
#    invjulian- finds the year day mon hr min sec from the julian date
#    julianday   - finds the julian date given year, mon day, hr, min, sec
#
#  references    :
#    vallado       2007, 292, Alg 32, Ex 5-4
#
# [utmoonrise, utmoonset, moonphaseang, error] = moonrise(jd, latgd, lon)
# -----------------------------------------------------------------------------

def moonrise2(jd=None, latgd=None, lon=None, show=None):
    # ------------------------  implementation   ------------------
    error = 'ok'
    # -------------- for once for moonrise (1), ) set (2) ---------
# -------------- make sure lon is within +- 180 deg -----------
    if (lon > np.pi):
        lon = lon - 2.0 * np.pi

    if (lon < - np.pi):
        lon = lon + twopi

    try1 = 1

    opt = 1

    while (opt <= 2):

        year, month, day, hr, min, sec = stu.invjday(jd, 0.0)
        jdtemp, jdtempf = stu.jday(year, month, day, 0, 0, 0.0)
        if (opt == 1):
            uttemp = (6.0 - lon / 15.0) / 24.0
        else:
            uttemp = (18.0 + lon / 15.0) / 24.0
        if (try1 == 2):
            uttemp = 0.5
        i = 0
        tn = uttemp
        t = tn + 10.0
        jdtemp = jdtemp + jdtempf + uttemp
        deltaut = 10.0
        ktlim = 15
        if show == 'y':
            print('   ecplon    ecllat       l          m           n       rtasc          decl        gmst       mgha      dgha    lhan     lhan1   deltaut     tn    deltaut    jd    tn \n' % ())
        while ((np.abs(deltaut) >= 8e-05) and (i < ktlim)):

            ttdb = (jdtemp + jdtempf - 2451545.0) / 36525.0
            eclplong = (218.32 + 481267.8813 * ttdb
                        + 6.29 * np.sin((134.9 + 477198.85 * ttdb) * deg2rad)
                        - 1.27 * np.sin((259.2 - 413335.38 * ttdb) * deg2rad)
                        + 0.66 * np.sin((235.7 + 890534.23 * ttdb) * deg2rad)
                        + 0.21 * np.sin((269.9 + 954397.7 * ttdb) * deg2rad)
                        - 0.19 * np.sin((357.5 + 35999.05 * ttdb) * deg2rad)
                        - 0.11 * np.sin((186.6 + 966404.05 * ttdb) * deg2rad))
            eclplat = (5.13 * np.sin((93.3 + 483202.03 * ttdb) * deg2rad)
                       + 0.28 * np.sin((228.2 + 960400.87 * ttdb) * deg2rad)
                       - 0.28 * np.sin((318.3 + 6003.18 * ttdb) * deg2rad)
                       - 0.17 * np.sin((217.6 - 407332.2 * ttdb) * deg2rad))
            eclplong = np.fmod(eclplong * deg2rad, twopi)
            eclplat = np.fmod(eclplat * deg2rad, twopi)
            if show == 'y':
                print('%2i %2i %11.7f %11.7f ' % (opt, i, eclplong / deg2rad, eclplat / deg2rad))
            obliquity = 23.439291 - 0.0130042 * ttdb
            obliquity = obliquity * deg2rad
            # ------- find the geocentric direction cosines -------
            l = np.cos(eclplat) * np.cos(eclplong)
            m = (np.cos(obliquity) * np.cos(eclplat) * np.sin(eclplong)
                 - np.sin(obliquity) * np.sin(eclplat))
            n = (np.sin(obliquity) * np.cos(eclplat) * np.sin(eclplong)
                 + np.cos(obliquity) * np.sin(eclplat))
            if show == 'y':
                print(' %11.7f %11.7f %11.7f ' % (l, m, n))
            rtasc = math.atan2(m, l)
            # - check that rtasc is in the same quadrant as eclplong
#        if (eclplong < 0.0)
#            eclplong = eclplong + twopi
#        end
            if (np.abs(eclplong - rtasc) > np.pi * 0.5):
                rtasc = rtasc + 0.5 * np.pi * np.rint(0.5 + (eclplong - rtasc)
                                                      / (0.5 * np.pi))
            decl = np.arcsin(n)
            lst, gmst = stu.lstime(lon, jdtemp)
            if show == 'y':
                print(' %d %d %d ' % (rtasc, decl, gmst / deg2rad))
            moonghan = gmst - rtasc
            #cdav
            lha = moonghan + lon
            #           # find elevation directly at this time
#           el = asin(sin(latgd) * sin(rtasc) + cos(latgd) * cos(rtasc) * cos(lha))
            if (i == 0):
                dgha = 347.81 * deg2rad
            else:
                dgha = (moonghan - moongha) / deltaut
            #        if (dgha < 0.0)
#            dgha = dgha + twopi / abs(deltaut)
#        end
            if show == 'y':
                print(' %11.7f  %11.7f ' % (moonghan / deg2rad, dgha / deg2rad))
            lhan = (0.00233 - (np.sin(latgd) * np.sin(decl))
                    / (np.cos(latgd) * np.cos(decl)))
            if show == 'y':
                print(' %11.7f rad ' % (lhan))
            #fprintf('lhan  #11.7f deg \n', lhan/deg2rad)
#       if (lhan > 1.0)
#           lhan = 0.0
#       end
#       if (lhan < -1.0)
#           lhan = -1.0
#       end
            if np.abs(lhan) > 1.0:
                deltaut = 1.0
                uttemp = 9999.9
                i = ktlim
            else:
                lhan = np.arccos(lhan)
                if (opt == 1):
                    lhan = twopi - lhan
                if show == 'y':
                    print(' %11.7f ' % (lhan / deg2rad))
                deltaut = (lhan - lha) / dgha
                if show == 'y':
                    print(' % 11.7f  %11.7f ' % (deltaut, tn))
                #                if (abs(deltaut) > 0.5)
#                    if (abs(dgha) > 0.001)
                if (deltaut < - 0.5):
                    deltaut = deltaut + twopi / dgha
                    #                            if (abs(deltaut) > 0.51)
#                                uttemp = -9999.9  # day before
#                                i = ktlim # end this trial
#                            end
                else:
                    if (deltaut > 0.5):
                        deltaut = deltaut - twopi / dgha
                        #                                if (abs(deltaut) > 0.51)
#                                    i = ktlim # end this trial needed
#                                    uttemp = 9999.9
#                                end
                #                    else
#                        error = 'error2 dgha is too small'
#                        fprintf('dgha too small #11.7f ', dgha)
#                    end
#                end
                #                 if (abs(deltaut) > 0.5)
#                     if (abs(dgha) > 0.001)
#                         if (deltaut < 0.0) # event is day before
#                             deltaut = deltaut + twopi / dgha
#                             if (abs(deltaut) > 0.51)
#                                 uttemp = -9999.9  # day before
#                                 i = ktlim # end this trial
#                             end
#                         else
#                             deltaut = deltaut - twopi / dgha
#                             if (abs(deltaut) > 0.51)
#                                 i = ktlim # end this trial needed
#                                 uttemp = 9999.9
#                             end
#                         end
#                     else
#                         error = 'error2 dgha is too small'
#                         fprintf('dgha too small #11.7f ', dgha)
                #                     end
#                 end
                if uttemp + deltaut < 0.0:
                    deltaut = deltaut - 1.0
                    i = ktlim
                    uttemp = - 9999.9
                i = i + 1
                moongha = moonghan
            # change tn-t to uttemp for iteration towards 0.0
            jdtemp = jdtemp + jdtempf + deltaut
            uttemp = uttemp + deltaut
            tn = uttemp
            if show == 'y':
                print(' % 11.7f  %14.4f  %11.7f \n'
                      % (deltaut, jdtemp + jdtempf, tn))

        uttemp = tn * 24.0
        if (tn > 2.0):
            i = ktlim
            uttemp = 9999.99
        #if (uttemp > 24.0) && (uttemp < 9999.0)
# fprintf('rem #11.7f ', uttemp)
        #  uttemp = np.fmod(uttemp, 24.0)
# fprintf('#11.7f ', uttemp)
#end
#        if (uttemp < 0.0)
#            uttemp = uttemp + 24.0
#        end
#    if (uttemp > 900)
#        uttemp = 24.0
#        fprintf('this cannot happen\n')
#    end
        if (opt == 1):
            utmoonrise = uttemp
        if (opt == 2):
            utmoonset = uttemp
        # update the iteration and check for solution
        try1 = try1 + 1
        if ((i > ktlim - 1) and (try1 < 3)):
            if show == 'y':
                print('try1 #2 %4d \n' % (opt))
                print(' uttemp  %11.7f  hrs \n' % (uttemp))
        else:
            #             if ((i > ktlim-1) && (try1 > 2))
#                 if (opt == 1)
#                     error = 'no rise'
#                     uttemp = 9999.99
#                 end
#                 if (opt == 2)
#                     error = 'no set'
#                     uttemp = 9997.99
#                 end
#             end
            if show == 'y':
                print(' uttemp  %11.7f  hrs \n' % (uttemp))
            opt = opt + 1
            try1 = 1


    # ------------- determine phase angle of the moon -------------
    meanlong = 280.4606184 + 36000.77005361 * ttdb
    meanlong = np.fmod(meanlong, 360.0)
    meananomaly = 357.5277233 + 35999.05034 * ttdb
    meananomaly = np.fmod(meananomaly * deg2rad, twopi)
    if (meananomaly < 0.0):
        meananomaly = twopi + meananomaly

    loneclsun = (meanlong + 1.914666471 * np.sin(meananomaly)
                 + 0.019994643 * np.sin(2.0 * meananomaly))
    loneclmoon = (218.32 + 481267.8813 * ttdb
                  + 6.29 * np.sin((134.9 + 477198.85 * ttdb) * deg2rad)
                  - 1.27 * np.sin((259.2 - 413335.38 * ttdb) * deg2rad)
                  + 0.66 * np.sin((235.7 + 890534.23 * ttdb) * deg2rad)
                  + 0.21 * np.sin((269.9 + 954397.7 * ttdb) * deg2rad)
                  - 0.19 * np.sin((357.5 + 35999.05 * ttdb) * deg2rad)
                  - 0.11 * np.sin((186.6 + 966404.05 * ttdb) * deg2rad))
    loneclmoon = np.fmod(loneclmoon, 360.0)
    moonphaseang = loneclmoon - loneclsun
    if (moonphaseang < 0.0):
        moonphaseang = 360.0 + moonphaseang

    return utmoonrise, utmoonset, moonphaseang, error

# -----------------------------------------------------------------------------
#
#                           function moonriset
#
#  this function finds the universal time for moonrise and moonset given the
#    day and site location.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  revisions
#                - david vallado                                 25 jan 2011
#
#  inputs          description                    range / units
#    jd          - julian date                    days from 4713 bc
#    latgd       - site latitude (south -)        -65 to 65 rad
#    lon         - site longitude (west -)        -2pi to 2pi rad
#
#  outputs       :
#    utmoonrise  - universal time of moonrise     hrs
#    utmoonset   - universal time of moonset      hrs
#    moonphaseang- phase angle of the moon        deg
#    error       - error parameter
#
#  locals        :
#    moonangle   - angle between the moon vector
#                  and a point on the earth       rad
#    jdtemp      - julian date for moonrise/set   days from 4713 bc
#    uttemp      - temporary ut time              days
#    tut1        - julian centuries from the
#                  jan 1, 2000 12 h epoch (ut1)
#    rtasc       - right ascension                rad
#    decl        - declination                    rad
#    meanlonmoon -                                rad
#    meananomaly -                                rad
#    eclplong    - longitude of the ecliptic      rad
#    obliquity   - obliquity of the ecliptic      rad
#    rmoon
#    rmoonrs
#    rv
#    rhosat
#    tsry1
#    l, m, n     - direction cosines
#    eclplat
#    moongha, moonghan
#    dgha, dghan
#    lhan
#    lst
#    deltaut, deltautn
#    t, tn
#    hzparal
#    loneclsun
#    loneclmoon
#    ttdb
#    gst         - for 0 h utc of each day        rad
#    lha         - local hour angle               rad
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - hour                           0 .. 23
#    min         - minute                         0 .. 59
#    sec         - second                         0.0  .. 59.999
#    opt         - idx to for rise and set calc    1, 2
#
#  coupling      :
#    invjulian- finds the year day mon hr min sec from the julian date
#    julianday   - finds the julian date given year, mon day, hr, min, sec
#
#  references    :
#    vallado       2007, 292, Alg 32, Ex 5-4
#
# [utmoonrise, utmoonset, moonphaseang, error] = moonrise(jd, latgd, lon)
# -----------------------------------------------------------------------------

def moonrise3(jd=None, latgd=None, lon=None, show=None):
    # ------------------------  implementation   ------------------
    error = 'ok'
    # -------------- for once for moonrise (1), ) set (2) ---------
# -------------- make sure lon is within +- 180 deg -----------
    if (lon > np.pi):
        lon = lon - 2.0 * np.pi

    if (lon < - np.pi):
        lon = lon + twopi

    #   try1 = 1 # try another approach on the current option
# opt = 1-2 for rise/set
    for opt in range(3):
        tolerance = 0.00035
        GHA = 0.0
        LHAn = 0.0
        jdtemp = jd
        deltaUT = 0.001
        deltaJD = 0.0
        loopCount = 1
        while ((np.abs(deltaUT) > tolerance) and (loopCount < 10)):

            ttdb = (jdtemp - 2451545.0) / 36525.0
            eclplong = (218.32 + 481267.8813 * ttdb
                        + 6.29 * np.sin((134.9 + 477198.85 * ttdb) * deg2rad)
                        - 1.27 * np.sin((259.2 - 413335.38 * ttdb) * deg2rad)
                        + 0.66 * np.sin((235.7 + 890534.23 * ttdb) * deg2rad)
                        + 0.21 * np.sin((269.9 + 954397.7 * ttdb) * deg2rad)
                        - 0.19 * np.sin((357.5 + 35999.05 * ttdb) * deg2rad)
                        - 0.11 * np.sin((186.6 + 966404.05 * ttdb) * deg2rad))
            eclplat = (5.13 * np.sin((93.3 + 483202.03 * ttdb) * deg2rad)
                       + 0.28 * np.sin((228.2 + 960400.87 * ttdb) * deg2rad)
                       - 0.28 * np.sin((318.3 + 6003.18 * ttdb) * deg2rad)
                       - 0.17 * np.sin((217.6 - 407332.2 * ttdb) * deg2rad))
            eclplong = np.fmod(eclplong * deg2rad, twopi)
            eclplat = np.fmod(eclplat * deg2rad, twopi)
            if show == 'y':
                print('%2i %2i%8.5f  %11.7f %11.7f ' %
                      (opt, loopCount, ttdb, eclplong / deg2rad,
                       eclplat / deg2rad))
            obliquity = 23.439291 - 0.0130042 * ttdb
            obliquity = obliquity * deg2rad
            # ------- find the geocentric direction cosines -------
            l = np.cos(eclplat) * np.cos(eclplong)
            m = (np.cos(obliquity) * np.cos(eclplat) * np.sin(eclplong)
                 - np.sin(obliquity) * np.sin(eclplat))
            n = (np.sin(obliquity) * np.cos(eclplat) * np.sin(eclplong)
                 + np.cos(obliquity) * np.sin(eclplat))
            rtasc = math.atan2(m, l)
            # - check that rtasc is in the same quadrant as eclplong
            if (eclplong < 0.0):
                eclplong = eclplong + twopi
            if (np.abs(eclplong - rtasc) > np.pi * 0.5):
                rtasc = (rtasc + 0.5 * np.pi
                         * np.rint(0.5 + (eclplong - rtasc) / (0.5 * np.pi)))
            decl = np.arcsin(n)
            lst, gmst = stu.lstime(lon, jdtemp)
            GHAn = gmst - rtasc
            LHA = GHAn + lon
            if (loopCount == 1):
                deltaGHA = 347.81 * deg2rad
            else:
                deltaGHA = ((GHAn - GHA) / deltaUT)
            if (deltaGHA < 0.0):
                deltaGHA = deltaGHA + twopi / np.abs(deltaUT)
            cosLHAn = ((0.00233 - np.sin(latgd) * np.sin(decl))
                       / (np.cos(latgd) * np.cos(decl)))
            if show == 'y':
                print(' coslha %8.5f lmn %11.7f %11.7f %11.7f \n'
                      % (cosLHAn, l, m, n))
            if (np.abs(cosLHAn) > 1.0):
                # No event on this day advance to the next
                deltaUT = 1
                if show == 'y':
                    print('nothing Advancing one day \n' % ())
                print('a' % ())
            else:
                LHAn = np.arccos(cosLHAn)
                if (opt == 1):
                    LHAn = twopi - LHAn
                #if (debugging) out.printf("LHAn - LHA = #f\n", LHAn - LHA)
                deltaUT = (LHAn - LHA) / deltaGHA
                if (deltaUT < - 0.5):
                    deltaUT = deltaUT + twopi / deltaGHA
                else:
                    if (deltaUT > 0.5):
                        deltaUT = deltaUT - twopi / deltaGHA
                if (deltaJD + deltaUT < 0.0):
                    deltaUT = deltaUT + 1
                    print('A' % ())
                    if show == 'y':
                        print('Advancing one day \n' % ())
                GHA = GHAn
            jdtemp = jdtemp + deltaUT
            deltaJD = deltaJD + deltaUT
            loopCount = loopCount + 1
            if show == 'y':
                print('%2i %2i %11.7f hrs %8.5f '
                      % (opt, loopCount, deltaUT * 24, deltaJD * 24))
                print('rtasc %11.7f  decl %8.5f gmst %8.5f GHAn %8.5f dGHA %8.5f  LHAn %8.5f jdtemp %8.5f\n'
                      % (rtasc / deg2rad, decl / deg2rad, gmst / deg2rad,
                         GHAn / deg2rad, deltaGHA / deg2rad, LHAn / deg2rad,
                         jdtemp))

        if opt == 1:
            utmoonrise = deltaJD * 24
        if opt == 2:
            utmoonset = deltaJD * 24

    if show == 'y':
        print('rise %11.7f  set %11.7f  \n' % (utmoonrise, utmoonset))



    # ------------- determine phase angle of the moon -------------
    meanlong = 280.4606184 + 36000.77005361 * ttdb
    meanlong = np.fmod(meanlong, 360.0)
    meananomaly = 357.5277233 + 35999.05034 * ttdb
    meananomaly = np.fmod(meananomaly * deg2rad, twopi)
    if (meananomaly < 0.0):
        meananomaly = twopi + meananomaly

    loneclsun = (meanlong + 1.914666471 * np.sin(meananomaly)
                 + 0.019994643 * np.sin(2.0 * meananomaly))
    loneclmoon = (218.32 + 481267.8813 * ttdb
                  + 6.29 * np.sin((134.9 + 477198.85 * ttdb) * deg2rad)
                  - 1.27 * np.sin((259.2 - 413335.38 * ttdb) * deg2rad)
                  + 0.66 * np.sin((235.7 + 890534.23 * ttdb) * deg2rad)
                  + 0.21 * np.sin((269.9 + 954397.7 * ttdb) * deg2rad)
                  - 0.19 * np.sin((357.5 + 35999.05 * ttdb) * deg2rad)
                  - 0.11 * np.sin((186.6 + 966404.05 * ttdb) * deg2rad))
    loneclmoon = np.fmod(loneclmoon, 360.0)
    moonphaseang = loneclmoon - loneclsun
    if (moonphaseang < 0.0):
        moonphaseang = 360.0 + moonphaseang

    return utmoonrise, utmoonset, moonphaseang, error

# ------------------------------------------------------------------------------
#
#                           function moonill
#
#  this function calculates the illumination due to the moon.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    f           - moon hase angle                rad
#    moonel      - moon elevation                 rad
#
#  outputs       :
#    moonill     - moon illumination
#
#  locals        :
#                -
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2001, 295-297, eq 5-9
#
# [moonillum] = moonill (f, moonel)
# ------------------------------------------------------------------------------

def moonill(f=None, moonel=None):
    # ------------------------  implementation   ------------------
    x = moonel / 90.0
    g = 1.0
    if (moonel >= 20):
        l0 = - 1.95
        l1 = 4.06
        l2 = - 4.24
        l3 = 1.56
    else:
        if (((moonel >= 5.0) and (moonel < 20.0))):
            l0 = - 2.58
            l1 = 12.58
            l2 = - 42.58
            l3 = 59.06
        else:
            if (((moonel > - 0.8) and (moonel < 5.0))):
                l0 = - 2.79
                l1 = 24.27
                l2 = - 252.95
                l3 = 1321.29
            else:
                l0 = 0.0
                l1 = 0.0
                l2 = 0.0
                l3 = 0.0
                f = 0.0
                g = 0.0

    l1 = l0 + l1 * x + l2 * x * x + l3 * x * x * x
    l2 = (- 0.00868 * f - 2.2e-09 * f * f * f * f)
    #       hzparal = 0.9508 + 0.0518*cos((134.9+477198.85*ttdb)*deg2rad)
#                + 0.0095*cos((259.2-413335.38*ttdb)*deg2rad)
#                + 0.0078*cos((235.7+890534.23*ttdb)*deg2rad)
#                + 0.0028*cos((269.9+954397.70*ttdb)*deg2rad)   { deg }
#       hzparal = realmod(hzparal*deg2rad, twopi)
#       l3 = (2.0* power(10.0, (hzparal * rad2deg / 0.951))*g) { use g to eliminate neg el passes }

    moonillum = 10.0 ** (l1 + l2)
    if (((moonillum < - 1e+36) or (moonillum > 0.999))):
        moonillum = 0.0

    return moonillum

# ------------------------------------------------------------------------------
#
#                           function checkhitearth
#
#  this function checks to see if the trajectory hits the earth during the
#    transfer.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  inputs          description                    range / units
#    altPad      - pad for alt above surface       km
#    r1          - initial position vector of int  km
#    v1t         - initial velocity vector of trns km/s
#    r2          - final position vector of int    km
#    v2t         - final velocity vector of trns   km/s
#    nrev        - number of revolutions           0, 1, 2, ...
#
#  outputs       :
#    hitearth    - is earth was impacted           'y' 'n'
#    hitearthstr - is earth was impacted           "y - radii" "no"
#
#  locals        :
#    sme         - specific mechanical energy
#    rp          - radius of perigee               km
#    a           - semimajor axis of transfer      km
#    ecc         - eccentricity of transfer
#    p           - semi-paramater of transfer      km
#    hbar        - angular momentum vector of
#                  transfer orbit
#
#  coupling      :
#    dot         - dot product of vectors
#    mag         - magnitude of a vector
#    cross       - cross product of vectors
#
#  references    :
#    vallado       2013, 503, alg 60
#
# [hitearth, hitearthstr] = checkhitearth (altpad, r1, v1t, r2, v2t, nrev)
# ------------------------------------------------------------------------------

def checkhitearth(altpad=None, r1=None, v1t=None, r2=None, v2t=None, nrev=None):
    # --------------------------  implementation   -----------------
    show = 'n'
    mu = 398600.4418
    hitearth = 'n'
    hitearthstr = 'no'
    magr1 = smu.mag(r1)
    magr2 = smu.mag(r2)
    rpad = 6378.137 + altpad
    #fprintf(1, 'mr1 #11.7f mr2 #11.7f ', magr1, magr2)
# check whether Lambert transfer trajectory hits the Earth
# this check may not be needed depending on input data
    if (magr1 < rpad or magr2 < rpad):
        # hitting earth already at start or stop point
        hitearth = 'y'
        hitearthstr = hitearth + ' initradii'
        if show == 'y':
            print('hitearth? %s \n' % (hitearthstr))
    else:
        rdotv1 = np.dot(r1, v1t)
        rdotv2 = np.dot(r2, v2t)
        # Solve for a
        ainv = 2.0 / magr1 - smu.mag(v1t) ** 2 / mu
        # Find ecos(E)
        ecosea1 = 1.0 - magr1 * ainv
        ecosea2 = 1.0 - magr2 * ainv
        #fprintf(1, 'ec1 #11.7f ec2 #11.7f ainv #11.7f ', ecosea1, ecosea2, ainv)
        # Determine radius of perigee
# 4 distinct cases pass thru perigee
# nrev > 0 you have to check
        if (nrev > 0):
            a = 1.0 / ainv
            # elliptical orbit
            if (a > 0.0):
                esinea1 = rdotv1 / np.sqrt(mu * a)
                ecc = np.sqrt(ecosea1 * ecosea1 + esinea1 * esinea1)
            else:
                # hyperbolic orbit
                esinea1 = rdotv1 / np.sqrt(mu * np.abs(- a))
                ecc = np.sqrt(ecosea1 * ecosea1 - esinea1 * esinea1)
            rp = a * (1.0 - ecc)
            if (rp < rpad):
                hitearth = 'y'
                hitearthstr = hitearth + 'Sub_Earth_nrev'
            # nrev = 0, 3 cases:
# heading to perigee and ending after perigee
# both headed away from perigee, but end is closer to perigee
# both headed toward perigee, but start is closer to perigee
        else:
            if ((rdotv1 < 0.0 and rdotv2 > 0.0) or
                    (rdotv1 > 0.0 and rdotv2 > 0.0 and ecosea1 < ecosea2) or
                    (rdotv1 < 0.0 and rdotv2 < 0.0 and ecosea1 > ecosea2)):
                # parabola
                if (np.abs(ainv) <= 1e-10):
                    hbar = math.cross(r1, v1t)
                    magh = smu.mag(hbar)
                    rp = magh * magh * 0.5 / mu
                    if (rp < rpad):
                        hitearth = 'y'
                        hitearthstr = hitearth + ' Sub_Earth_para'
                else:
                    # for both elliptical & hyperbolic
                    a = 1.0 / ainv
                    esinea1 = rdotv1 / np.sqrt(mu * np.abs(a))
                    if (ainv > 0.0):
                        ecc = np.sqrt(ecosea1 * ecosea1 + esinea1 * esinea1)
                    else:
                        ecc = np.sqrt(ecosea1 * ecosea1 - esinea1 * esinea1)
                    if (ecc < 1.0):
                        rp = a * (1.0 - ecc)
                        if (rp < rpad):
                            hitearth = 'y'
                            hitearthstr = hitearth + ' Sub_Earth_ell'
                    else:
                        # hyperbolic heading towards the earth
                        if (rdotv1 < 0.0 and rdotv2 > 0.0):
                            rp = a * (1.0 - ecc)
                            if (rp < rpad):
                                hitearth = 'y'
                                hitearthstr = hitearth + ' Sub_Earth_hyp'
                    #   fprintf(1, 'hitearth? #s rp #11.7f  #11.7f km \n', hitearthstr, rp*6378.137, rpad*6378.137)
            if show == 'y':
                print('hitearth? %s rp %11.7f km \n' % (hitearth, rp * 6378.0))

    return hitearth, hitearthstr


def ShadowEntryExit(RSun=None, rp=None, a=None, ecc=None, incl=None,
                    raan=None, argp=None, nu=None, mu=None):
    Een = 0.0
    Eex = 0.0
    # Semi-Parameter
    p = a * (1.0 - ecc ** 2)
    # Eccentric Anomaly:
    sinE = (np.sin(np.pi/180*nu) * np.sqrt(1 - ecc ** 2)) / (1 + ecc * np.cos(np.pi/180*nu))
    cosE = (ecc + np.cos(np.pi/180*nu)) / (1 + ecc * np.cos(np.pi/180*nu))
    # (5.3)
    Px = (np.cos(np.pi/180*raan) * np.cos(np.pi/180*argp)
          - np.sin(np.pi/180*raan) * np.sin(np.pi/180*argp)
          * np.cos(np.pi/180*incl))
    Py = (np.cos(np.pi/180*raan) * np.sin(np.pi/180*argp)
          + np.sin(np.pi/180*raan) * np.cos(np.pi/180*argp)
          * np.cos(np.pi/180*incl))
    Pz = np.sin(np.pi/180*raan) * np.sin(np.pi/180*incl)
    P_ = np.array([Px, Py, Pz])
    Qx = (- np.sin(np.pi/180*raan) * np.cos(np.pi/180*argp)
          - np.cos(np.pi/180*raan) * np.sin(np.pi/180*argp)
          * np.cos(np.pi/180*incl))
    Qy = (- np.sin(np.pi/180*raan) * np.sin(np.pi/180*argp)
          + np.cos(np.pi/180*raan) * np.cos(np.pi/180*argp)
          * np.cos(np.pi/180*incl))
    Qz = np.cos(np.pi/180*raan) * np.sin(np.pi/180*incl)
    Q_ = np.array([Qx, Qy, Qz])
    # (5.6)
    beta = np.dot(P_, RSun) / np.linalg.norm(RSun)
    zeta = np.dot(Q_, RSun) / np.linalg.norm(RSun)
    A0 = (((rp / p) ** 4) * (ecc ** 4)
          - 2 * ((rp / p) ** 2) * (zeta ** 2 - beta ** 2) * ecc ** 2
          + (beta ** 2 + zeta ** 2) ** 2)
    A1 = (4 * ((rp / p) ** 4) * ecc ** 3
          - 4 * ((rp / p) ** 2) * (zeta ** 2 - beta ** 2) * ecc)
    A2 = (6 * ((rp / p) ** 4) * ecc ** 2
          - 2 * ((rp / p) ** 2) * (zeta ** 2 - beta ** 2)
          - 2 * ((rp / p) ** 2) * (1 - zeta ** 2) * ecc ** 2
          + 2 * (zeta ** 2 - beta ** 2) * (1 - zeta ** 2)
          - 4 * (beta ** 2) * (zeta ** 2))
    A3 = (4 * ((rp / p) ** 4) * ecc
          - 4 * ((rp / p) ** 2) * (1 - zeta ** 2) * ecc)
    A4 = ((rp / p) ** 4 - 2 * ((rp / p) ** 2) * (1 - zeta ** 2)
          + (1 - zeta ** 2) ** 2)
    r1r, __, r2r, __, r3r, __, r4r, __ = smu.quartic(A0, A1, A2, A3, A4, 'R')
    nu = np.zeros(4)
    check = np.zeros(4)
    nu[0] = np.degrees(np.arccos(r1r))
    nu[1] = np.degrees(np.arccos(r2r))
    nu[2] = np.degrees(np.arccos(r3r))
    nu[3] = np.degrees(np.arccos(r4r))
    k = 1
    for incl in range(np.size(nu)):
        check[incl] = (beta * np.cos(np.pi/180*nu[incl])
                       + zeta * np.sin(np.pi/180*nu[incl]))
        if check[incl] < 0:
            nugood = nu[incl]
            k = k + 1
            print('----------------------------------------------------------------------')
            print('Valid Eccentric Anomaly (beta*cos(nu) + zeta*sin(nu) < 0): ', str(nu[incl]))
            before = (A0 * np.cos(np.pi/180*nugood - 0.01) ** 4
                      + A1 * np.cos(np.pi/180*nugood - 0.01) ** 3
                      + A2 * np.cos(np.pi/180*nugood - 0.01) ** 2
                      + A3 * np.cos(np.pi/180*nugood - 0.01) + A4)
            after = (A0 * np.cos(np.pi/180*nugood + 0.01) ** 4
                     + A1 * np.cos(np.pi/180*nugood + 0.01) ** 3
                     + A2 * np.cos(np.pi/180*nugood + 0.01) ** 2
                     + A3 * np.cos(np.pi/180*nugood + 0.01) + A4)
            if before < 0 and after > 0:
                print('Entering Shadow for This Eccentric Anomaly')
                Een = nu[incl]
            else:
                if before > 0 and after < 0:
                    print('Exiting Shadow for This Eccentric Anomaly')
                    Eex = nu[incl]
                else:
                    print('Error - whether entering or exiting is not clear')

    return Een, Eex

# -----------------------------------------------------------------------------
#
#                            procedure dspace
#
#   this procedure provides deep space contributions to mean elements for
#     perturbing third body.  these effects have been averaged over one
#     revolution of the sun and moon.  for earth resonance effects, the
#     effects have been averaged over no revolutions of the satellite.
#     (mean motion)
#
# Author:
#   Jeff Beck
#   beckja@alumni.lehigh.edu
#   1.0 (aug 6, 2006) - update for paper dav
# original comments from Vallado C++ version:
#   author        : david vallado                  719-573-2600   28 jun 2005
#
#   inputs        :
#     d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433       -
#     dedt        -
#     del1, del2, del3  -
#     didt        -
#     dmdt        -
#     dnodt       -
#     domdt       -
#     irez        - flag for resonance           0-none, 1-one day, 2-half day
#     argpo       - argument of perigee
#     argpdot     - argument of perigee dot (rate)
#     t           - time
#     tc          -
#     gsto        - gst
#     xfact       -
#     xlamo       -
#     no          - mean motion
#     atime       -
#     em          - eccentricity
#     ft          -
#     argpm       - argument of perigee
#     inclm       - inclination
#     xli         -
#     mm          - mean anomaly
#     xni         - mean motion
#     nodem      - right ascension of ascending node
#
#   outputs       :
#     atime       -
#     em          - eccentricity
#     argpm       - argument of perigee
#     inclm       - inclination
#     xli         -
#     mm          - mean anomaly
#     xni         -
#     nodem      - right ascension of ascending node
#     dndt        -
#     nm          - mean motion
#
#   locals        :
#     delt        -
#     ft          -
#     theta       -
#     x2li        -
#     x2omi       -
#     xl          -
#     xldot       -
#     xnddt       -
#     xndt        -
#     xomi        -
#
#   coupling      :
#     none        -
#
#   references    :
#     hoots, roehrich, norad spacetrack report #3 1980
#     hoots, norad spacetrack report #6 1986
#     hoots, schumacher and glover 2004
#     vallado, crawford, hujsak, kelso  2006
#  ----------------------------------------------------------------------------*/

def dspace(d2201=None, d2211=None, d3210=None, d3222=None,
           d4410=None, d4422=None, d5220=None, d5232=None,
           d5421=None, d5433=None, dedt=None, del1=None, del2=None,
           del3=None, didt=None, dmdt=None, dnodt=None, domdt=None,
           irez=None, argpo=None, argpdot=None, t=None, tc=None,
           gsto=None, xfact=None, xlamo=None, no=None, atime=None,
           em=None, argpm=None, inclm=None, xli=None, mm=None,
           xni=None, nodem=None, nm=None):
    fasx2 = 0.13130908
    fasx4 = 2.8843198
    fasx6 = 0.37448087
    g22 = 5.7686396
    g32 = 0.95240898
    g44 = 1.8014998
    g52 = 1.050833
    g54 = 4.4108898
    rptim = 0.0043752690880113
    stepp = 720.0
    stepn = - 720.0
    step2 = 259200.0
    # /* ----------- calculate deep space resonance effects ----------- */
    dndt = 0.0
    theta = np.fmod(gsto + tc * rptim, twopi)
    em = em + dedt * t
    inclm = inclm + didt * t
    argpm = argpm + domdt * t
    nodem = nodem + dnodt * t
    mm = mm + dmdt * t
    # //   sgp4fix for negative inclinations
# //   the following if statement should be commented out
# //  if (inclm < 0.0)
# // {
# //    inclm = -inclm
# //    argpm = argpm - pi
# //    nodem = nodem + pi
# //  }

    # /* - update resonances : numerical (euler-maclaurin) integration - */
# /* ------------------------- epoch restart ----------------------  */

    # //   sgp4fix for propagator problems
# //   the following integration works for negative time steps and periods
# //   the specific changes are unknown because the original code was so convoluted

    # // sgp4fix take out atime = 0.0 and fix for faster operation
    ft = 0.0
    if (irez != 0):
        # sgp4fix streamline check
        if ((atime == 0.0) or (t * atime <= 0.0) or (np.abs(t) < np.abs(atime))):
            atime = 0.0
            xni = no
            xli = xlamo
        # sgp4fix move check outside loop
        if (t >= 0.0):
            delt = stepp
        else:
            delt = stepn
        iretn = 381
        iret = 0
        while (iretn == 381):

            # /* ------------------- dot terms calculated ------------- */
# /* ----------- near - synchronous resonance terms ------- */
            if (irez != 2):
                xndt = (del1 * np.sin(xli - fasx2)
                        + del2 * np.sin(2.0 * (xli - fasx4))
                        + del3 * np.sin(3.0 * (xli - fasx6)))
                xldot = xni + xfact
                xnddt = (del1 * np.cos(xli - fasx2)
                         + 2.0 * del2 * np.cos(2.0 * (xli - fasx4))
                         + 3.0 * del3 * np.cos(3.0 * (xli - fasx6)))
                xnddt = xnddt * xldot
            else:
                # /* --------- near - half-day resonance terms -------- */
                xomi = argpo + argpdot * atime
                x2omi = xomi + xomi
                x2li = xli + xli
                xndt = (d2201 * np.sin(x2omi + xli - g22)
                        + d2211 * np.sin(xli - g22)
                        + d3210 * np.sin(xomi + xli - g32)
                        + d3222 * np.sin(- xomi + xli - g32)
                        + d4410 * np.sin(x2omi + x2li - g44)
                        + d4422 * np.sin(x2li - g44)
                        + d5220 * np.sin(xomi + xli - g52)
                        + d5232 * np.sin(- xomi + xli - g52)
                        + d5421 * np.sin(xomi + x2li - g54)
                        + d5433 * np.sin(- xomi + x2li - g54))
                xldot = xni + xfact
                xnddt = (d2201 * np.cos(x2omi + xli - g22)
                         + d2211 * np.cos(xli - g22)
                         + d3210 * np.cos(xomi + xli - g32)
                         + d3222 * np.cos(- xomi + xli - g32)
                         + d5220 * np.cos(xomi + xli - g52)
                         + d5232 * np.cos(- xomi + xli - g52)
                         + 2.0 * (d4410 * np.cos(x2omi + x2li - g44)
                                  + d4422 * np.cos(x2li - g44)
                                  + d5421 * np.cos(xomi + x2li - g54)
                                  + d5433 * np.cos(- xomi + x2li - g54)))
                xnddt = xnddt * xldot
            # /* ----------------------- integrator ------------------- */
# sgp4fix move end checks to end of routine
            if (np.abs(t - atime) >= stepp):
                iret = 0
                iretn = 381
            else:
                ft = t - atime
                iretn = 0
            if (iretn == 381):
                xli = xli + xldot * delt + xndt * step2
                xni = xni + xndt * delt + xnddt * step2
                atime = atime + delt

        nm = xni + xndt * ft + xnddt * ft * ft * 0.5
        xl = xli + xldot * ft + xndt * ft * ft * 0.5
        if (irez != 1):
            mm = xl - 2.0 * nodem + 2.0 * theta
            dndt = nm - no
        else:
            mm = xl - nodem - argpm + theta
            dndt = nm - no
        nm = no + dndt

    return atime, em, argpm, inclm, xli, mm, xni, nodem, dndt, nm


#  -----------------------------------------------------------------------------
#
#                            procedure dpper
#
#   this procedure provides deep space long period periodic contributions
#     to the mean elements.  by design, these periodics are zero at epoch.
#     this used to be dscom which included initialization, but it's really a
#     recurring function.
#
# Author:
#   Jeff Beck
#   beckja@alumni.lehigh.edu
#   1.0 (aug 7, 2006) - update for paper dav
# original comments from Vallado C++ version:
#   author        : david vallado                  719-573-2600   28 jun 2005
#
#   inputs        :
#     e3          -
#     ee2         -
#     peo         -
#     pgho        -
#     pho         -
#     pinco       -
#     plo         -
#     se2 , se3 , Sgh2, Sgh3, Sgh4, Sh2, Sh3, Si2, Si3, Sl2, Sl3, Sl4 -
#     t           -
#     xh2, xh3, xi2, xi3, xl2, xl3, xl4 -
#     zmol        -
#     zmos        -
#     ep          - eccentricity                           0.0 - 1.0
#     inclo       - inclination - needed for lyddane modification
#     nodep       - right ascension of ascending node
#     argpp       - argument of perigee
#     mp          - mean anomaly
#
#   outputs       :
#     ep          - eccentricity                           0.0 - 1.0
#     inclp       - inclination
#     nodep       - right ascension of ascending node
#     argpp       - argument of perigee
#     mp          - mean anomaly
#
#   locals        :
#     alfdp       -
#     betdp       -
#     cosip  , sinip  , cosop  , sinop  ,
#     dalf        -
#     dbet        -
#     dls         -
#     f2, f3      -
#     pe          -
#     pgh         -
#     ph          -
#     pinc        -
#     pl          -
#     sel   , ses   , sghl  , sghs  , shl   , shs   , sil   , sinzf , sis   ,
#     sll   , sls
#     xls         -
#     xnoh        -
#     zf          -
#     zm          -
#
#   coupling      :
#     none.
#
#   references    :
#     hoots, roehrich, norad spacetrack report #3 1980
#     hoots, norad spacetrack report #6 1986
#     hoots, schumacher and glover 2004
#     vallado, crawford, hujsak, kelso  2006
#  ----------------------------------------------------------------------------*/

def dpper(e3=None, ee2=None, peo=None, pgho=None, pho=None,
          pinco=None, plo=None, se2=None, se3=None, sgh2=None,
          sgh3=None, sgh4=None, sh2=None, sh3=None, si2=None,
          si3=None, sl2=None, sl3=None, sl4=None, t=None,
          xgh2=None, xgh3=None, xgh4=None, xh2=None, xh3=None,
          xi2=None, xi3=None, xl2=None, xl3=None, xl4=None,
          zmol=None, zmos=None, inclo=None, init=None, ep=None,
          inclp=None, nodep=None, argpp=None, mp=None, opsmode=None):
    # change to variable passed in
    # /* ---------------------- constants ----------------------------- */
    zns = 1.19459e-05
    zes = 0.01675
    znl = 0.00015835218
    zel = 0.0549
    # /* --------------- calculate time varying periodics ----------- */
    zm = zmos + zns * t
    # // be sure that the initial call has time set to zero
    if (init == 'y'):
        zm = zmos

    zf = zm + 2.0 * zes * np.sin(zm)
    sinzf = np.sin(zf)
    f2 = 0.5 * sinzf * sinzf - 0.25
    f3 = - 0.5 * sinzf * np.cos(zf)
    ses = se2 * f2 + se3 * f3
    sis = si2 * f2 + si3 * f3
    sls = sl2 * f2 + sl3 * f3 + sl4 * sinzf
    sghs = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf
    shs = sh2 * f2 + sh3 * f3
    zm = zmol + znl * t
    if (init == 'y'):
        zm = zmol

    zf = zm + 2.0 * zel * np.sin(zm)
    sinzf = np.sin(zf)
    f2 = 0.5 * sinzf * sinzf - 0.25
    f3 = - 0.5 * sinzf * np.cos(zf)
    sel = ee2 * f2 + e3 * f3
    sil = xi2 * f2 + xi3 * f3
    sll = xl2 * f2 + xl3 * f3 + xl4 * sinzf
    sghl = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf
    shll = xh2 * f2 + xh3 * f3
    pe = ses + sel
    pinc = sis + sil
    pl = sls + sll
    pgh = sghs + sghl
    ph = shs + shll
    if (init == 'n'):
        #   //  0.2 rad = 11.45916 deg
        pe = pe - peo
        pinc = pinc - pinco
        pl = pl - plo
        pgh = pgh - pgho
        ph = ph - pho
        inclp = inclp + pinc
        ep = ep + pe
        sinip = np.sin(inclp)
        cosip = np.cos(inclp)
        # /* ----------------- apply periodics directly ------------ */
#   //  sgp4fix for lyddane choice
#   //  strn3 used original inclination - this is technically feasible
#   //  gsfc used perturbed inclination - also technically feasible
#   //  probably best to readjust the 0.2 limit value and limit discontinuity
#   //  use next line for original strn3 approach and original inclination
#   //  if (inclo >= 0.2)
#   //  use next line for gsfc version and perturbed inclination
        if (inclp >= 0.2):
            ph = ph / sinip
            pgh = pgh - cosip * ph
            argpp = argpp + pgh
            nodep = nodep + ph
            mp = mp + pl
        else:
            # /* ---- apply periodics with lyddane modification ---- */
            sinop = np.sin(nodep)
            cosop = np.cos(nodep)
            alfdp = sinip * sinop
            betdp = sinip * cosop
            dalf = ph * cosop + pinc * cosip * sinop
            dbet = - ph * sinop + pinc * cosip * cosop
            alfdp = alfdp + dalf
            betdp = betdp + dbet
            nodep = np.fmod(nodep, twopi)
            # sgp4fix for afspc written intrinsic functions
# nodep used without a trigonometric function ahead
            if ((nodep < 0.0) and (opsmode == 'a')):
                nodep = nodep + twopi
            xls = mp + argpp + cosip * nodep
            dls = pl + pgh - pinc * nodep * sinip
            xls = xls + dls
            xnoh = nodep
            nodep = math.atan2(alfdp, betdp)
            # sgp4fix for afspc written intrinsic functions
# nodep used without a trigonometric function ahead
            if ((nodep < 0.0) and (opsmode == 'a')):
                nodep = nodep + twopi
            if (np.abs(xnoh - nodep) > np.pi):
                if (nodep < xnoh):
                    nodep = nodep + twopi
                else:
                    nodep = nodep - twopi
            mp = mp + pl
            argpp = xls - mp - cosip * nodep

    return ep, inclp, nodep, argpp, mp

# ------------------------------------------------------------------------------
#
#                           function sunill
#
#  this function calculates the illumination due to the sun.
#
#  author        : david vallado                  719-573-2600    9 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    jd          - julian date                    days
#    lat         - location latitude              rad
#    lon         - location longitdue             rad
#    sunaz       - sun azimuth                    rad
#    sunel       - sun elevation                  rad
#
#  outputs       :
#    sunillum    - sun illumination
#
#  locals        :
#                -
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2001, 295-297, eq 5-9
#
# [sunillum] = sunill   (jd, lat, lon, sunaz, sunel)
# ------------------------------------------------------------------------------

def sunill(jd=None, lat=None, lon=None, sunaz=None, sunel=None):
    # -------------------------  implementation   -----------------
    rsun, srtasc, sdecl = sun(jd)

    lst, gst = stu.lstime(lon, jd)
    lha = lst - srtasc
    sunel = (np.arcsin(np.sin(sdecl) * np.sin(lat)
                       + np.cos(sdecl) * np.cos(lat) * np.cos(lha)))
    sinv = (- np.sin(lha) * np.cos(sdecl)
            * np.cos(lat) / (np.cos(sunel) * np.cos(lat)))
    cosv = ((np.sin(sdecl) - np.sin(sunel) * np.sin(lat))
            / (np.cos(sunel) * np.cos(lat)))
    sunaz = math.atan2(sinv, cosv)
    sunel = sunel / deg2rad
    if (sunel > - 18.01):
        x = sunel / 90.0
        if (sunel >= 20):
            l0 = 3.74
            l1 = 3.97
            l2 = - 4.07
            l3 = 1.47
        else:
            if (((sunel >= 5.0) and (sunel < 20.0))):
                l0 = 3.05
                l1 = 13.28
                l2 = - 45.98
                l3 = 64.33
            else:
                if (((sunel >= - 0.8) and (sunel < 5.0))):
                    l0 = 2.88
                    l1 = 22.26
                    l2 = - 207.64
                    l3 = 1034.3
                else:
                    if (((sunel >= - 5.0)and (sunel < - 0.8))):
                        l0 = 2.88
                        l1 = 21.81
                        l2 = - 258.11
                        l3 = - 858.36
                    else:
                        if (((sunel >= - 12.0) and (sunel < - 5.0))):
                            l0 = 2.7
                            l1 = 12.17
                            l2 = - 431.69
                            l3 = - 1899.83
                        else:
                            if (((sunel >= - 18.0) and (sunel < - 12.0))):
                                l0 = 13.84
                                l1 = 262.72
                                l2 = 1447.42
                                l3 = 2797.93
                            else:
                                l0 = 0.0
                                l1 = 0.0
                                l2 = 0.0
                                l3 = 0.0
        l1 = l0 + l1 * x + l2 * x * x + l3 * x * x * x
        sunillum = 10.0 ** l1
        if (((sunillum < - 1e+36) or (sunillum > 999.999))):
            sunillum = 0.0
    else:
        sunillum = 0.0

    return sunillum



# ------------------------------------------------------------------------------
#
#                           function target
#
#  this function accomplishes the targeting problem using kepler/pkepler &
#    lambert.
#
#  author        : david vallado                  719-573-2600    1 mar 2001
#
#  inputs          description                    range / units
#    rint        - initial position vector of int km
#    vint        - initial velocity vector of int km/s
#    rtgt        - initial position vector of tgt km
#    vtgt        - initial velocity vector of tgt km/s
#    dm          - direction of motion for gauss  'l', 's'
#    kind        - type of propagator             'k', 'p'
#    dtsec       - time of flight to the int      s
#
#  outputs       :
#    v1t         - initial transfer velocity vec  km/s
#    v2t         - final transfer velocity vec    km/s
#    dv1         - initial change velocity vec    km/s
#    dv2         - final change velocity vec      km/s
#    error       - error flag from gauss          'ok', ...
#
#  locals        :
#    transnormal - cross product of trans orbit   km
#    intnormal   - cross product of int orbit     km
#    r1tgt       - position vector after dt, tgt  km
#    v1tgt       - velocity vector after dt, tgt  km/s
#    rirt        - rint(4) * r1tgt(4)
#    cosdeltanu  - cosine of deltanu              rad
#    sindeltanu  - sine of deltanu                rad
#    deltanu     - angl between position vectors rad
#    i           - index
#
#  coupling      :
#    kepler      - find r and v at future time
#    lambertu - find velocity vectors at each end of transfer
#    lncom2      - linear combination of two vectors and constants  ###untrue
#
#  references    :
#    vallado       2001, 468-474, alg 58
#
# [v1t, v2t, dv1, dv2] = target (rint, vint, rtgt, vtgt, dm, kind, dtsec, ndot, nddot)
# ------------------------------------------------------------------------------


def target (rint, vint, rtgt, vtgt, dm, kind, dtsec, ndot, nddot):

# -------------------------  implementation   -------------------------

    # ----------- propogate target forward by time ----------------
    if (kind =='k'):
        kepler
        r1tgt, v1tgt, errmsg = kepler(rtgt, vtgt, 0.0)
        print("kepler returned: ", errmsg)
    elif (kind =='p'):
        r1tgt, v1tgt = pkepler(rtgt, vtgt, dtsec, ndot, nddot)
    else:
        print("Invalid target input for kind: ", kind)
        return

    # ----------- calculate transfer orbit between r's ------------
    v1t, v2t, error1 = lambertu(rint, r1tgt, dm, 'n', dtsec)

    if (error1 == '      ok'):
        dv1 = -vint + v1t
        dv2 = v1tgt - v2t
    else:
        print(' err %s ', error1)
        v1t = np.zeros((3))
        v2t = np.zeros((3))
        dvt = np.zeros((3))
        dv2 = np.zeros((3))

    return v1t, v2t, dv1, dv2


# -----------------------------------------------------------------------------
#
#                           procedure findatwbatwa
#
# this procedure finds the a and b matrices for the differential correction
#   problem.  remember that it isn't critical for the propagations to use
#   the highest fidelity techniques because we're only trying to find the
#   "slope". k is an index that allows us to do multiple rows at once. it's
#   used for both the b and a matrix calculations.
#
#  algorithm     : find the a and b matrices by accumulation to reduce matrix sizes
#                  calculate the matrix combinations
#                  atw is found without matrix operations to avoid large matrices
#
#  author        : david vallado                  719-573-2600    6 aug 2008
#
#  inputs          description                                range / units
#    firstob     - number of observations
#    lastob      - number of observations
#    statesize   - size of state                                6 , 7
#    percentchg  - amount to modify the vectors
#                  by in finite differencing
#    deltaamtchg - tolerance for small value in
#                  finite differencing                         0.0000001
#    whichconst  - parameter for sgp4 constants                wgs72, wgs721, wgs84
#    satrec      - structure of satellite parameters for TLE
#    obsrecfile    - array of records containing:
#                  senum, jd, rsvec, obstype,
#                  rng, az, el, drng, daz, del,
#                  trtasc, tdecl data
#    statetype   - type of elements (equinoctial, etc)         'e', 't'
#    scalef      - scale factor to limit the size of the state vector
#    xnom        - state vector                                varied
#
#  outputs       :
#    atwa        - atwa matrix
#    atwb        - atwb matrix
#    atw         - atw matrix
#    b           - b matrix, the residuals
#    drng2       - range residual squared
#    daz2        - azimuth residual squared
#    del2        - elevation residual squared
#    ddrng2      - range rate residual squared
#    ddaz2       - azimuth rate residual squared
#    ddel2       - elevation rate residual squared
#    dtrtasc2    - topocentric right ascension residual squared
#    dtdecl2     - topocentric declination residual squared
#    dx2         - x position residual squared
#    dy2         - y position residual squared
#    dz2         - z position residual squared
#    dxdot2      - xdot position residual squared
#    dydot2      - ydot position residual squared
#    dzdot2      - zdot position residual squared
#
#  locals        :
#    rnom        - nom position vector at epoch                km
#    vnom        - nom velocity vector at epoch                km/s
#    a           - a matrix
#    indobs      -
#    at          -
#    w1          -
#    w2          -
#    w3          -
#    lst         -
#    gst         -
#    dtsec        -
#    deltaamt    -
#    rngpert     -
#    azpert      - modified azimuth                            -2pi to 2pi
#    elpert      - modified azimuth                            -pi/2 to pi/2
#    drng        -
#    daz         -
#    del         -
#    error       -
#
#  coupling      :
#    findsenptr  - find sensor data
#    rv_razel    - find r and v given range, az, el, and rates
#    rv_tradec   - find r and v given topocentric rtasc and decl
#
#  references    :
#    vallado       2013, 753-765
# --------------------------------------------------------------------------- */

def findatwaatwb(firstobs=None, lastobs=None, obsrecarr=None,
                 statesize=None, percentchg=None, deltaamtchg=None,
                 xnom=None):
    # --------------------- initialize parameters ------------------
    if obsrecarr[0]['obstype'] == 0:
        indobs = 1
    else:
        if obsrecarr[0]['obstype'] == 2:
            indobs = 3
        else:
            indobs = 2

    drng2 = 0.0
    daz2 = 0.0
    del2 = 0.0
    dtrtasc2 = 0.0
    dtdecl2 = 0.0
    atwaacc = np.array([[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])
    atwbacc = np.array([[0], [0], [0], [0], [0], [0]])
    # ------------- reset these since they will accumulate ---------
# zero out matrices
    atwa = np.zeros((statesize, statesize))
    atwb = np.zeros((statesize, 1))
    rnom = np.zeros(3)
    vnom = np.zeros(3)
    a = np.zeros((indobs, statesize))
    b = np.zeros((indobs, 1))
    atw = np.zeros((statesize, indobs))
    # ------------------- loop through all the observations ------------------
    for obsktr in range(firstobs, lastobs + 1):
        currobsrec = obsrecarr[obsktr]
        #printf("ob #2i rsecef0 #11.5f jd #11.5f dtmin #8.3f hr #3i rng #11.7f \n",
#         i, currobsrec['rsecef'][0], currobsrec['jd'], currobsrec['dtmin'], currobsrec['hr'], currobsrec['rng'])  # ritrf
        # ----------------- determine sensor characteristics -----------------
#         rs[0] = currobsrec['rsecef'][0]
#         rs(2) = currobsrec['rsecef'](2)
#         rs(3) = currobsrec['rsecef'](3)
        # temporary sensor for now
#  getsensorparams(currobsrec.sennum, currsenrec)
        # --------- propagate the nominal vector to the epoch time -----------
#         sgp4 (whichconst, satrec,  currobsrec.dtmin, rteme, vteme)
        dtsec = (currobsrec['time'] + currobsrec['timef'] - obsrecarr[0]['time'] - obsrecarr[0]['timef']) * 86400.0
        rnom[0] = xnom[0, 0]
        rnom[1] = xnom[1, 0]
        rnom[2] = xnom[2, 0]
        vnom[0] = xnom[3, 0]
        vnom[1] = xnom[4, 0]
        vnom[2] = xnom[5, 0]
        # [reci1, veci1] = kepler (rnom, vnom, dtsec)
        reci1, veci1 = pkepler(rnom, vnom, dtsec, 0, 0)
        # ------------------------- find b matrix ----------------------------
        if currobsrec['obstype'] != 3:
            rngnom, aznom, elnom, drngnom, daznom, delnom = sc.rv2razel(reci1.T, veci1.T, currobsrec['latgd'], currobsrec['lon'], currobsrec['alt'], currobsrec['ttt'], currobsrec['jdut1'], 0.0, currobsrec['xp'], currobsrec['yp'], 2, 0.0, 0.0)
        else:
            rngnom, trtascnom, tdeclnom, drngnom, dtrtascnom, dtdeclnom = sc.rv2tradc(reci1.T, veci1.T, currobsrec['latgd'], currobsrec['lon'], currobsrec['alt'], currobsrec['ttt'], currobsrec['jdut1'], 0.0, currobsrec['xp'], currobsrec['yp'], 2, 0.0, 0.0)
        if 0 == (currobsrec['obstype']):
            b[0, 0] = currobsrec['rng'] - rngnom
        elif 1 == (currobsrec['obstype']):
            b[0, 0] = currobsrec['az'] - aznom
            #fix for 0-360...
            if (np.abs(b[0, 0]) > np.pi):
                b[0, 0] = b[0, 0] - math.sgn(b[0, 0]) * 2.0 * np.pi
            b[1, 0] = currobsrec['el'] - elnom
        elif 2 == (currobsrec['obstype']):
            b[0, 0] = currobsrec['rng'] - rngnom
            b[1, 0] = currobsrec['az'] - aznom
            # fix for 0-360...
            if np.abs(b[1, 0]) > np.pi:
                b[1, 0] = b[1, 0] - np.sign(b[1, 0]) * 2.0 * np.pi
            b[2, 0] = currobsrec['el'] - elnom
        elif 3 == (currobsrec['obstype']):
            b[0, 0] = currobsrec['trtasc'] - trtascnom
            # fix for 0-360...
            if (np.abs(b[0, 0]) > np.pi):
                b[0, 0] = b[0, 0] - np.sign(b[0, 0]) * 2.0 * np.pi
            b[1, 0] = currobsrec['tdecl'] - tdeclnom
        #printf("rnom #11.5f #11.5f #11.5f #8.3f #8.3f #8.3f #8.3f \n",
#              rteme[0], rteme[1], rteme[2], rngnom, aznom * rad2deg, elnom * rad2deg, currobsrec['rng'])  # ritrf
        # ------------------------ find a matrix -----------------------------
# ------------- reset the perturbed vector to the nominal ------------
#  satrecp = satrec
        xnomp = np.copy(xnom)
        rnomp = np.zeros(3)
        vnomp = np.zeros(3)
        # ----- perturb each element in the state (elements or vectors) ------
        for j in range(statesize):
            deltaamt, xnomp = smu.finitediff(j, percentchg, deltaamtchg, xnom)
            #   sgp4 (whichconst, satrecp,  currobsrec['dtmin'], r3, v3)
            dtsec = (currobsrec['time'] + currobsrec['timef']
                     - obsrecarr[0]['time'] - obsrecarr[0]['timef']) * 86400
            rnomp[0] = xnomp[0, 0]
            rnomp[1] = xnomp[1, 0]
            rnomp[2] = xnomp[2, 0]
            vnomp[0] = xnomp[3, 0]
            vnomp[1] = xnomp[4, 0]
            vnomp[2] = xnomp[5, 0]
            #  [reci3, veci3] = kepler (rnomp, vnomp, dtsec)
            reci3, veci3 = pkepler(rnomp, vnomp, dtsec, 0, 0)
            # teme to itrf if observation type
#  if (currobsrec['obstype'] ~= 4)
#      mfme = currobsrec['hr'](j) * 60 + currobsrec['min'](j) + currobsrec['sec'](j)/60.0
#      findeopparam (currobsrec['jd'], mfme, interp, eoparr, jdeopstart,
#                     dut1, dat, lod, xp, yp, ddpsi, ddeps,
#                     iaudx, dy, icrsx, y, s, deltapsi, deltaeps)
#      [ut1, tut1, jdut1, utc, tai, tt, ttt, jdtt, tdb, ttdb, jdtdb ] ...
#       = convtime (year, mon, day, hr, min, sec, timezone, dut1, dat)
#      iau76fk5_itrf_teme(ritrf, vitrf, aitrf, eFrom, r3, v3, ateme, ttt, xp, yp, jdut1, lod, trans)
#  end # if obstype
            if (currobsrec['obstype'] == 3):
                trrpert, trtascpert, tdeclpert, tdrrpert, tdrtascpert, tddeclpert =\
                      sc.rv2tradc(reci3.T, veci3.T, currobsrec['latgd'],
                                  currobsrec['lon'], currobsrec['alt'],
                                  currobsrec['ttt'], currobsrec['jdut1'], 0.0,
                                  currobsrec['xp'], currobsrec['yp'], 2, 0.0, 0.0)
            elif (currobsrec['obstype'] == 4):
                bstarpert = currobsrec['bstar'] * (1.0 + percentchg)
            else:
                rngpert, azpert, elpert, drhopert, dazpert, delpert =\
                    sc.rv2razel(reci3.T, veci3.T, currobsrec['latgd'],
                                currobsrec['lon'], currobsrec['alt'],
                                currobsrec['ttt'], currobsrec['jdut1'], 0.0,
                                currobsrec['xp'], currobsrec['yp'], 2, 0.0, 0.0)
                    #        fprintf(1, 'rnom  #16.8f #16.8f #16.8f #16.8f #16.8f #16.8f #16.8f km \n', reci1, rngnom, aznom, elnom, deltaamt)
#        fprintf(1, 'rpert #16.8f #16.8f #16.8f #16.8f #16.8f #16.8f #16.8f km \n', reci3, rngpert, azpert, elpert, dtsec)

            if 0 == currobsrec['obstype']:
                a[0, j] = (rngpert - rngnom) / deltaamt
            elif 1 == currobsrec['obstype']:
                a[0, j] = (azpert - aznom) / deltaamt
                a[1, j] = (elpert - elnom) / deltaamt
            elif 2 == currobsrec['obstype']:
                a[0, j] = (rngpert - rngnom) / deltaamt
                a[1, j] = (azpert - aznom) / deltaamt
                a[2, j] = (elpert - elnom) / deltaamt
            elif 3 == currobsrec['obstype']:
                a[0, j] = (trtascpert - trtascnom) / deltaamt
                a[1, j] = (tdeclpert - tdeclnom) / deltaamt
            #printf("rpert #11.5f #11.5f #11.5f #8.3f #8.3f #8.3f \n",
#              r3[0], r3[0], r3[2], rngpert, azpert * rad2deg, elpert * rad2deg)  # ritrf
            # ----------------- reset the modified vector --------------------
#satrecp = satrec
            xnomp = np.copy(xnom)
        # ----------------- now form the matrix combinations -----------------
        at = a.T
        # ------------------------- assign weights ---------------------------
        if 0 == (currobsrec['obstype']):
            w1 = 1.0 / (currobsrec['noiserng'] ** 2)
            rng2 = rng2 + b[0, 0] * b[0, 0] * w1
        elif 1 == (currobsrec['obstype']):
            w1 = 1.0 / (currobsrec['noiseaz'] ** 2)
            w2 = 1.0 / (currobsrec['noiseel'] ** 2)
            daz2 = daz2 + b[0, 0] * b[0, 0] * w1
            del2 = del2 + b[1, 0] * b[1, 0] * w2
        elif 2 == (currobsrec['obstype']):
            w1 = 1.0 / (currobsrec['noiserng'] * currobsrec['noiserng'])
            w2 = 1.0 / (currobsrec['noiseaz'] * currobsrec['noiseaz'])
            w3 = 1.0 / (currobsrec['noiseel'] * currobsrec['noiseel'])
            drng2 = drng2 + b[0, 0] * b[0, 0] * w1
            daz2 = daz2 + b[1, 0] * b[1, 0] * w2
            del2 = del2 + b[2, 0] * b[2, 0] * w3
        elif 3 == (currobsrec['obstype']):
            w1 = 1.0 / (currobsrec['noisetrtasc'] * currobsrec['noisetrtasc'])
            w2 = 1.0 / (currobsrec['noisetdecl'] * currobsrec['noisetdecl'])
            dtrtasc2 = dtrtasc2 + b[0, 0] * b[0, 0] * w1
            dtdecl2 = dtdecl2 + b[1, 0] * b[1, 0] * w2

        # for rowc in range(statesize):
        #     for colc in range(indobs):
        #         weight = 1.0
        #         if 0 == (colc):
        #             weight = w1
        #         elif 1 == (colc):
        #             weight = w2
        #         elif 2 == (colc):
        #             weight = w3
        #         # w4-w7 never assigned!
        #         # elif 3 == (colc):
        #         #     weight = w4
        #         # elif 4 == (colc):
        #         #     weight = w5
        #         # elif 5 == (colc):
        #         #     weight = w6
        #         # elif 6 == (colc):
        #         #     weight = w7
        #         atw[rowc, colc] = at[rowc, colc] * weight

        w = np.array([[w1, 0, 0], [0, w2, 0], [0, 0, w3]])
        atw = at @ w

        # ----------------- find the atwa / atwb matrices --------------------
        atwaacc = atw @ a
        #matmult(atw, a, atwaacc, statesize, indobs, statesize)
        atwbacc = atw @ b
        #matmult(atw, b, atwbacc, statesize, indobs, 1)
        # ------------------- accumulate the matricies -----------------------
        for r in range(statesize):
            for c in range(statesize):
                atwa[r, c] = atwaacc[r, c] + atwa[r, c]
        c = 0
        for r in range(statesize):
            atwb[r, c] = atwbacc[r, c] + atwb[r, c]
        #writeexpmat("atwa ", atwa, statesize, statesize)
#writeexpmat("atwb ", atwb, statesize, 1)
#writeexpmat("a ", a, indobs, statesize)
#writemat("b ", b, indobs, 1)

    return atwa, atwb, atw, b, drng2, daz2, del2




##############################################################################################################
##############################################################################################################
##############################################################################################################

if __name__ == '__main__':

  jd = 60206

  sunillum = sunill(jd, 30.0, 45.0, 10.0, 20.0)
  print("sunill returned: ", sunillum)

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
  raan = 80.0
  rp = 1.0

  x, y, z = makeorbitrv(2455545.0, 'j',
                      [5003.400903511, -3817.812007872, 4720.200666830],
                      [5.489294908, 3.005055561, -3.39013016])
  print("makeorbitrv returned ", x, y, z)

  rsun, rtasc, decl = sun(jd)
  print(rsun, rtasc, decl)
  Een, Eex = ShadowEntryExit(rsun, rp, a, ecc, incl, raan, argp, nu, mu)

  r1 = np.array([4e6, 5e6, 6e6])
  r2 = np.array([1e6, 2e6, 3e6])
  v1 = np.array([5.0, -15.0,  -2.0])
  v2 = np.array([4.0, -10.0,  -3.0])

  # these seem to have old versions of kepler and lambertu and dnw -jmb
  #r, v = target(r1, v1, r2, v2, 'L', 'k', 30, 2, 1)
  #print("target returned: ", r, v)
  #r, v = target(r1, v1, r2, v2, 'L', 'p', 30, 2, 1)
  #print("target returned: ", r, v)

  v1t = np.array([1e3, 1e3, 99e3])
  v2t = np.array([1e3, 50e3, 50e3])
  altpad = 10.

  hitearth, hitearthstr = checkhitearth (altpad, r1, v1t, r2, v2t, 3)
  print("checkhitearth returned ", hitearth, hitearthstr)

  axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti = iau06in()
  print(axs0, a0xi, ays0, a0yi, ass0, a0si, apn, apni, appl, appli, agst, agsti)


  st = iau06era (jd)
  print(st)

  ttt = (jd - 2451545.0)/ 36525.0  #0.34698738576
  ddpsi = math.pi
  ddeps = math.pi

  #prec, psia, wa, ea, xa = precess (ttt, "50")
  #print("precess : " , prec, psia, wa, ea, xa)

  #deltapsi, trueeps, meaneps, omega, nut = nutation (ttt, ddpsi, ddeps)
  #print("nutation: " , deltapsi, trueeps, meaneps, omega, nut)

  #los = sight(r1, r2, 's')
  #print("los is ", los)

  #lit = light(r1, jd, 's')
  #print("lit is ", lit)




























