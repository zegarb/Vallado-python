# ------------------------------------------------------------------------------
#
#                           function repeatgt
#
#  This subroutine calculates repeat ground tracks.
#
#  author        : david vallado                  719-573-2600   2013
#
#  revisions     : michael courville                             June 2024
#
#  inputs          description                          range / units
#    krev2rep      - revs to repeat (crossing points)   km
#    kday2rep      - days to repeat (frequency)         days
#    ecc           - eccentrictiy                       0.0 to
#    incl          - inclination                        rad
#
#  outputs       :
#    a             - semimajor axis                     km
#
#  references    :
#    vallado       2013, pg 873, alg 71
#
# [a] = repeatgt ( b );
# ------------------------------------------------------------------------------

from space_constants import mu, re, j2
import math

def repeatgt(krev2rep, kday2rep, ecc, incl):
    revpday = krev2rep / kday2rep
    n = revpday * 2 * math.pi / 86400.0
    a = (mu * (1 / n) ** 2) ** 0.33333
    lonshiftrev = (2 * math.pi * kday2rep) / krev2rep

    # iteration of 10 is arbitrary
    for i in range(10):
        p = a * (1 - ecc**2)
        raandot = -1.5 * j2 * (re / p)**2 * math.cos(incl)
        lonshiftperiod = (2 * math.pi*raandot) / n
        londelta = lonshiftrev + lonshiftperiod
        # Another 2pi? yes or no?
        n = (4 * math.pi * math.pi / 86400.0) / londelta
        #n = (2 * math.pi / 86400.0) / londelta
        a = (mu * (1 / n) ** 2) ** 0.33333

    return a


# Old Code
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
