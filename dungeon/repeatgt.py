# ------------------------------------------------------------------------------
#
#                           function repeatgt
#
#  these subroutine calculates repeat ground tracks
#
#  author        : david vallado                  719-573-2600   10 feb 2006
#
#  revisions
#
#  inputs          description                    range / units
#    r           - ecef position vector           km
#
#  outputs       :
#    latgc       - geocentric latitude            -pi to pi rad
#
#  locals        :
#    temp        - diff between geocentric/
#
#  coupling      :
#
#  references    :
#    vallado       2001, 174-179, alg 12 and alg 13, ex 3-3
#
# [a] = repeatgt ( b );
# ------------------------------------------------------------------------------

import numpy as np
from space_constants import *

# def repeatgt(r = None):

#     a = 6400.0
#     ecc = 0.0001
#     incl = 98.0 * deg2rad
#     nanom = np.sqrt(mu / (a * a * a))
#     nanom
#     #        nnodal = 1;

#     p = a * (1.0 - ecc * ecc)
#     # -------------------------  implementation   -----------------
#     raanrate = - 1.5 * j2 * nanom * np.cos(incl) / (p * p)
#     revs2rep = 107
#     revspday = 16.0
#     days2rpt = revs2rep / revspday
#     days2rpt
#     lonshift = days2rpt * 2 * np.pi / revs2rep
#     lonshift * rad2deg
#     nnodal = oearth * revspday
#     nnodal
#     # old a422 way
#     rp = 160.0
#     period2b = 2.0 * np.pi * np.sqrt(mu / (a * a * a))
#     periodnew = period2b + period2b * raanrate / oearth
#     anew = (mu * (periodnew / oearth) ** 2) ** (1.0 / 3.0)
#     enew = 1 - rp / anew

#     # ------------- iterate to find geodetic latitude -------------
#     i = 1
#     olddelta = latgd + 10.0
#     while ((np.abs(olddelta - latgd) >= small) and (i < 10)):

#         olddelta = latgd
#         sintemp = np.sin(latgd)
#         c = re / (np.sqrt(1.0 - eccearthsqrd * sintemp * sintemp))
#         latgd = np.arctan((r[2] + c * eccearthsqrd * sintemp) / temp)
#         i = i + 1


#     return latgd

def repeatgt(rev2rep,day2rep,ecc,incl):
    revpday = rev2rep/day2rep
    n = revpday * 2 * np.pi / 86400.0
    a = (mu * (1 / n) ** 2) ** 0.33333
    lonshiftrev = (2 * np.pi * day2rep)/rev2rep

    for i in range(10):
        p = a * (1 - ecc**2)
        raandot = -1.5 * j2 * (re/p)**2 *cos(incl)
        lonshiftperiod = (2*np.pi*raandot) / n
        londelta = lonshiftrev + lonshiftperiod
        # Another 2pi? yes or no?
        #n = (4 * np.pi * np.pi / 86400.0) / londelta
        n = (2 * np.pi / 86400.0) / londelta
        a = (mu * (1 / n) ** 2) ** 0.33333

    return a

