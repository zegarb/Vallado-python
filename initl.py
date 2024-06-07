# -----------------------------------------------------------------------------
#
#                            procedure initl
#
#   this procedure initializes the spg4 propagator. all the initialization is
#     consolidated here instead of having multiple loops inside other routines.
#
# Author:
#   Jeff Beck
#   beckja@alumni.lehigh.edu
#   1.0 (aug 7, 2006) - update for paper dav
#   1.1 nov 16, 2007 - update for better compliance
# original comments from Vallado C++ version:
#   author        : david vallado                  719-573-2600   28 jun 2005
#
#   inputs        :
#     ecco        - eccentricity                           0.0 - 1.0
#     epoch       - epoch time in days from jan 0, 1950. 0 hr
#     inclo       - inclination of satellite
#     no          - mean motion of satellite
#     satn        - satellite number
#
#   outputs       :
#     ainv        - 1.0 / a
#     ao          - semi major axis
#     con41       -
#     con42       - 1.0 - 5.0 cos(i)
#     cosio       - cosine of inclination
#     cosio2      - cosio squared
#     einv        - 1.0 / e
#     eccsq       - eccentricity squared
#     method      - flag for deep space                    'd', 'n'
#     omeosq      - 1.0 - ecco * ecco
#     posq        - semi-parameter squared
#     rp          - radius of perigee
#     rteosq      - square root of (1.0 - ecco*ecco)
#     sinio       - sine of inclination
#     gsto        - gst at time of observation               rad
#     no          - mean motion of satellite
#
#   locals        :
#     ak          -
#     d1          -
#     del         -
#     adel        -
#     po          -
#
#   coupling      :
#     gstime      - find greenwich sidereal time from the julian date
#
#   references    :
#     hoots, roehrich, norad spacetrack report #3 1980
#     hoots, norad spacetrack report #6 1986
#     hoots, schumacher and glover 2004
#     vallado, crawford, hujsak, kelso  2006
#  ----------------------------------------------------------------------------*/

import numpy as np
import spacetime_utils as stu
from space_constants import twopi


def initl(xke = None,j2 = None,ecco = None,epoch = None,inclo = None,no_kozai = None,opsmode = None):
    # /* -------------------- wgs-72 earth constants ----------------- */
#     // sgp4fix identify constants and allow alternate values
# global tumin mu radiusearthkm xke j2 j3 j4 j3oj2
    x2o3 = 2.0 / 3.0
    #   global opsmode

    # /* ------------- calculate auxillary epoch quantities ---------- */
    eccsq = ecco * ecco
    omeosq = 1.0 - eccsq
    rteosq = np.sqrt(omeosq)
    cosio = np.cos(inclo)
    cosio2 = cosio * cosio
    # /* ------------------ un-kozai the mean motion ----------------- */
    ak = (xke / no_kozai) ** x2o3
    d1 = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq)
    del_ = d1 / (ak * ak)
    adel = ak * (1.0 - del_ * del_ - del_ * (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0))
    del_ = d1 / (adel * adel)
    no_unkozai = no_kozai / (1.0 + del_)
    ao = (xke / no_unkozai) ** x2o3
    sinio = np.sin(inclo)
    po = ao * omeosq
    con42 = 1.0 - 5.0 * cosio2
    con41 = - con42 - cosio2 - cosio2
    ainv = 1.0 / ao
    einv = 1.0 / ecco
    posq = po * po
    rp = ao * (1.0 - ecco)
    method = 'n'
    # sgp4fix modern approach to finding sidereal time
    if (opsmode != 'a'):
        gsto = stu.gstime(epoch + 2433281.5)
    else:
        # sgp4fix use old way of finding gst
# count integer number of days from 0 jan 1970
        ts70 = epoch - 7305.0
        ids70 = int(np.floor(ts70 + 1e-08))
        tfrac = ts70 - ids70
        # find greenwich location at epoch
        c1 = 0.017202791694070362
        thgr70 = 1.7321343856509375
        fk5r = 5.075514194322695e-15
        c1p2p = c1 + twopi
        gsto = np.fmod(thgr70 + c1 * ids70 + c1p2p * tfrac + ts70 * ts70 * fk5r,twopi)

    if (gsto < 0.0):
        gsto = gsto + twopi

    return method,ainv,ao,con41,con42,cosio,cosio2,eccsq,omeosq,posq,rp,rteosq,sinio,gsto,no_unkozai
