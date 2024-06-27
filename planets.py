from space_constants import *
import spacemath_utils as smu



class Planet:
    def __init__(self):
        self.abbr = None #Planet abbrevation

    # Approximates semi-major axis given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_a(self, ttdb):
        print("Approximating semi-major axis")
        return None

    # Approximates eccentricity given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_ecc(self, ttdb):
        print("Approximating eccentricity")
        return None

    # Approximates inclination given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_incl(self, ttdb):
        print("Approximating inclination")
        return None

    # Approximates omega (right ascension) given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_omega(self, ttdb):
        print("Approximating right ascension of ascending node")
        return None

    # Approximates longitude of periapsis given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_argp1(self, ttdb):
        print("Approximating longitude of periapsis")
        return None

    # Approximates mean longitude given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_lonmean(self, ttdb):
        print("Approximating mean longitude")
        return None

    # To approximate heliocentric orbital elements (coes) referenced to the
    # mean equator and mean equinox of J2000 (IAU-76/FK5)
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

class Mercury(Planet):
    def __init__(self):
        self.abbr = 'me'


    # All approx calculations: From Meeus (1991:202–204), the ecliptic orbital
    # elements of the planets, referenced to the mean equator and mean equinox
    # # of J2000, are as follows.

    def approxcalc_a(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return 0.387098310 * au

    def approxcalc_ecc(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return 0.20563175 + 2.0406e-05 * ttdb - 2.84e-08 * ttdb**2 \
                     - 1.7e-10 * ttdb**3

    def approxcalc_incl(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (7.004986 - 0.0059516 * ttdb + 8.1e-07 * ttdb**2
            + 4.1e-08 * ttdb**3) * deg2rad

    def approxcalc_omega(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (48.330893 - 0.1254229 * ttdb - 8.833e-05 * ttdb**2
            - 1.96e-07 * ttdb**3) * deg2rad

    def approxcalc_argp1(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (77.456119 + 0.1588643 * ttdb - 1.343e-05 * ttdb**2
            + 3.9e-08 * ttdb**3) * deg2rad

    def approxcalc_lonmean(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (252.250906 + 149472.674635 * ttdb - 5.35e-06 * ttdb**2
            + 2.0e-09 * ttdb**3) * deg2rad

class Venus(Planet):
    def __init__(self):
        self.abbr = 'v'


    # All approx calculations: From Meeus (1991:202–204), the ecliptic orbital
    # elements of the planets, referenced to the mean equator and mean equinox
    # # of J2000, are as follows.

    def approxcalc_a(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return 0.723329820 * au

    def approxcalc_ecc(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return 0.00677188 - 4.7766e-05 * ttdb + 9.75e-08 * ttdb**2 \
                     + 4.4e-10 * ttdb**3

    def approxcalc_incl(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (3.394662 - 8.8568e-04 * ttdb - 3.244e-05 * ttdb**2 \
            + 4.4e-10 * ttdb**3) * deg2rad

    def approxcalc_omega(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (76.679920 - 0.278008 * ttdb - 1.38232e-03 * ttdb**2
            - 1.98e-07 * ttdb**3) * deg2rad

    def approxcalc_argp1(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (131.563707 + 0.004864 * ttdb - 0.00138232 * ttdb**2
            - 5.332e-06 * ttdb**3) * deg2rad

    def approxcalc_lonmean(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (181.979801 + 58517.815676 * ttdb + 1.65e-06 * ttdb**2
            - 2.0e-09 * ttdb**3) * deg2rad


class Earth(Planet):
    def __init__(self):
        self.abbr = 'e'


    # All approx calculations: From Meeus (1991:202–204), the ecliptic orbital
    # elements of the planets, referenced to the mean equator and mean equinox
    # # of J2000, are as follows.

    def approxcalc_a(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return 1.000001018 * au

    def approxcalc_ecc(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return 0.01670862 - 0.000042037 * ttdb - 0.0000001236 * ttdb**2 \
                     + 0.00000000004 * ttdb**3

    def approxcalc_incl(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (0.0000000 + 0.0130546 * ttdb - 0.00000931 * ttdb**2 \
            - 0.000000034 * ttdb**3) * deg2rad

    def approxcalc_omega(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (174.873174 - 0.2410908 * ttdb + 0.00004067 * ttdb**2
            - 0.000000034 * ttdb**3) * deg2rad

    def approxcalc_argp1(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (102.937348 + 0.3225557 * ttdb + 0.00015026 * ttdb**2
            + 0.000000478 * ttdb**3) * deg2rad

    def approxcalc_lonmean(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (100.466449 + 35999.372851 * ttdb - 0.00000568 * ttdb**2
            + 0.0e-09 * ttdb**3) * deg2rad


class Mars(Planet):
    def __init__(self):
        self.abbr = 'ma'


    # All approx calculations: From Meeus (1991:202–204), the ecliptic orbital
    # elements of the planets, referenced to the mean equator and mean equinox
    # # of J2000, are as follows.

    def approxcalc_a(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return 1.523679342 * au

    def approxcalc_ecc(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return 0.09340062 + 0.000090483 * ttdb - 9.75e-08 * ttdb**2 \
                     - 3.5e-10 * ttdb**3

    def approxcalc_incl(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (1.849726 - 0.0081479 * ttdb - 0.00002255 * ttdb**2 \
            - 0.000000027 * ttdb**3) * deg2rad

    def approxcalc_omega(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (49.558093 - 0.2949846 * ttdb - 0.00063993 * ttdb**2
            - 0.000002143 * ttdb**3) * deg2rad

    def approxcalc_argp1(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (336.060234 + 0.4438898 * ttdb - 0.00017321 * ttdb**2
            + 0.000000300 * ttdb**3) * deg2rad

    def approxcalc_lonmean(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (355.433275 + 19140.299331 * ttdb + 0.00000261 * ttdb**2
            - 0.000000003 * ttdb**3) * deg2rad

class Jupiter(Planet):
    def __init__(self):
        self.abbr = "j"


    # All approx calculations: From Meeus (1991:202–204), the ecliptic orbital
    # elements of the planets, referenced to the mean equator and mean equinox
    # # of J2000, are as follows.

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
            return 0.04849485 + 0.000163244 * ttdb - 4.719e-07 * ttdb**2 \
            - 1.97e-9 * ttdb**3

    def approxcalc_incl(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (1.30327 - 0.0019872 * ttdb + 3.318e-05 * ttdb**2
            + 9.2e-08 * ttdb**3) * deg2rad

    def approxcalc_omega(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (100.464441 + 0.1766828 * ttdb + 0.000903877 * ttdb**2
            - 7.032e-06 * ttdb**3) * deg2rad

    def approxcalc_argp1(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (14.331309 + 0.2155525 * ttdb + 0.00072252 * ttdb**2
            - 4.59e-06 * ttdb**3) * deg2rad

    def approxcalc_lonmean(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (34.351484 + 3034.9056746 * ttdb - 8.501e-05 * ttdb**2
            + 4e-09 * ttdb**3) * deg2rad


# Not complete (need to update parameters)
# class Saturn(Planet):
#     def __init__(self):
#         self.abbr = 's'


#     # All approx calculations: From Meeus (1991:202–204), the ecliptic orbital
#     # elements of the planets, referenced to the mean equator and mean equinox
#     # # of J2000, are as follows.

#     def approxcalc_a(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return 0.723329820 * au

#     def approxcalc_ecc(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return 0.00677188 - 4.7766e-05 * ttdb + 9.75e-08 * ttdb**2 \
#                      + 4.4e-10 * ttdb**3

#     def approxcalc_incl(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (3.394662 - 8.8568e-04 * ttdb - 3.244e-05 * ttdb**2 \
#             + 4.4e-10 * ttdb**3) * deg2rad

#     def approxcalc_omega(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (76.679920 - 0.278008 * ttdb - 1.38232e-03 * ttdb**2
#             - 1.98e-07 * ttdb**3) * deg2rad

#     def approxcalc_argp1(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (131.563707 + 0.004864 * ttdb - 0.00138232 * ttdb**2
#             - 5.332e-06 * ttdb**3) * deg2rad

#     def approxcalc_lonmean(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (181.979801 + 58517.815676 * ttdb + 1.65e-06 * ttdb**2
#             - 2.0e-09 * ttdb**3) * deg2rad


# # Not complete (need to update parameters)
# class Uranus(Planet):
#     def __init__(self):
#         self.abbr = 'u'


#     # All approx calculations: From Meeus (1991:202–204), the ecliptic orbital
#     # elements of the planets, referenced to the mean equator and mean equinox
#     # # of J2000, are as follows.

#     def approxcalc_a(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return 0.723329820 * au

#     def approxcalc_ecc(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return 0.00677188 - 4.7766e-05 * ttdb + 9.75e-08 * ttdb**2 \
#                      + 4.4e-10 * ttdb**3

#     def approxcalc_incl(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (3.394662 - 8.8568e-04 * ttdb - 3.244e-05 * ttdb**2 \
#             + 4.4e-10 * ttdb**3) * deg2rad

#     def approxcalc_omega(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (76.679920 - 0.278008 * ttdb - 1.38232e-03 * ttdb**2
#             - 1.98e-07 * ttdb**3) * deg2rad

#     def approxcalc_argp1(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (131.563707 + 0.004864 * ttdb - 0.00138232 * ttdb**2
#             - 5.332e-06 * ttdb**3) * deg2rad

#     def approxcalc_lonmean(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (181.979801 + 58517.815676 * ttdb + 1.65e-06 * ttdb**2
#             - 2.0e-09 * ttdb**3) * deg2rad

# # Not complete (need to update parameters)
# class Neptune(Planet):
#     def __init__(self):
#         self.abbr = 'n'


#     # All approx calculations: From Meeus (1991:202–204), the ecliptic orbital
#     # elements of the planets, referenced to the mean equator and mean equinox
#     # # of J2000, are as follows.

#     def approxcalc_a(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return 0.723329820 * au

#     def approxcalc_ecc(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return 0.00677188 - 4.7766e-05 * ttdb + 9.75e-08 * ttdb**2 \
#                      + 4.4e-10 * ttdb**3

#     def approxcalc_incl(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (3.394662 - 8.8568e-04 * ttdb - 3.244e-05 * ttdb**2 \
#             + 4.4e-10 * ttdb**3) * deg2rad

#     def approxcalc_omega(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (76.679920 - 0.278008 * ttdb - 1.38232e-03 * ttdb**2
#             - 1.98e-07 * ttdb**3) * deg2rad

#     def approxcalc_argp1(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (131.563707 + 0.004864 * ttdb - 0.00138232 * ttdb**2
#             - 5.332e-06 * ttdb**3) * deg2rad

#     def approxcalc_lonmean(self, ttdb):
#         if ttdb < -2.0 or ttdb > 0.5:
#             print("Error: Only accurate between 1800 and 2050")
#             return None
#         else:
#             return (181.979801 + 58517.815676 * ttdb + 1.65e-06 * ttdb**2
#             - 2.0e-09 * ttdb**3) * deg2rad