from space_constants import *
import spacemath_utils as smu

class Planet:
    def __init__(self):
        self.abbr = None

    def approxcalc_a(self, ttdb):
        print("Approximating semi-major axis")
        return None

    def approxcalc_ecc(self, ttdb):
        print("Approximating eccentricity")
        return None

    def approxcalc_incl(self, ttdb):
        print("Approximating inclination")
        return None

    def approxcalc_omega(self, ttdb):
        print("Approximating right ascension of ascending node")
        return None

    def approxcalc_argp1(self, ttdb):
        print("Approximating longitude of periapsis")
        return None

    def approxcalc_lonmean(self, ttdb):
        print("Approximating mean longitude")
        return None

    # To approximate heliocentric orbital elements for the mean equator
    # using IAU-76/FK5 (J-2000)
    def approxcalc_coe(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            a = self.approxcalc_a(ttdb)
            ecc = self.approxcalc_ecc(ttdb)
            p = a * (1.0 - ecc ** 2)
            incl = self.approxcalc_incl(ttdb)
            omega = self.approxcalc_omega(ttdb)
            argp1 = self.approxcalc_argp1(ttdb)
            lonmean = self.approxcalc_lonmean(ttdb)

            m = lonmean - argp1
            argp = argp1 - omega
            _ , nu = smu.newtonm(ecc,m)

            return p, ecc, incl, omega, argp, nu

class Jupiter(Planet):
    def __init__(self):
        self.abbr = "j"

    def approxcalc_a(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (5.202603191 + 1.913e-07 * ttdb) * au

    def approxcalc_ecc(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return 0.04849485 + 0.000163244 * ttdb - 4.719e-07 * ttdb ** 2 \
            - 1.97e-9 * ttdb ** 3

    def approxcalc_incl(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (1.30327 - 0.0019872 * ttdb + 3.318e-05 * ttdb ** 2
            + 9.2e-08 * ttdb ** 3) * deg2rad

    def approxcalc_omega(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (100.464441 + 0.1766828 * ttdb + 0.000903877 * ttdb ** 2
            - 7.032e-06 * ttdb ** 3) * deg2rad

    def approxcalc_argp1(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (14.331309 + 0.2155525 * ttdb + 0.00072252 * ttdb ** 2
            - 4.59e-06 * ttdb ** 3) * deg2rad

    def approxcalc_lonmean(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (34.351484 + 3034.9056746 * ttdb - 8.501e-05 * ttdb ** 2
            + 4e-09 * ttdb ** 3) * deg2rad


# Finds semi-major axis
# Only accurate between 1800-2050 AD

def approxcalc_a(self, ttdb):

    print(ttdb)

    if self.abbr == 'me':
        return 0.387098310 * au
    elif self.abbr == 'v':
        return 0.723329820 * au
    elif self.abbr == 'e':
        return 1.000001018 * au
    elif self.abbr == 'ma':
        return 1.523679342 * au
    elif self.abbr == 'j':
        return (5.202603191 + 0.0000001913 * ttdb) * au
    elif self.abbr == 's':
        return (9.554909596 - 0.0000021389 * ttdb) * au
    elif self.abbr == 'u':
        return (19.218466062 - 0.0000000372 * ttdb + \
                0.00000000098 * ttdb ** 2) * au
    elif self.abbr == 'n':
        return (30.110386869 - 0.0000001663 * ttdb + \
                0.00000000069 * ttdb ** 2) * au
    else:
        print("planet unknown")
        return None

