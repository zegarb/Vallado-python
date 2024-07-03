from space_constants import *
import spacemath_utils as smu
import numpy as np
import os


class Planet:
    def __init__(self):
        self.abbr = None #Planet abbrevation
        # Poly coefficients of coe changes per century
        self.a0 = 0.0
        self.a1 = 0.0
        self.a2 = 0.0
        self.ecc0 = 0.0
        self.ecc1 = 0.0
        self.ecc2 = 0.0
        self.ecc3 = 0.0
        self.incl0 = 0.0
        self.incl1 = 0.0
        self.incl2 = 0.0
        self.incl3 = 0.0
        self.omega0 = 0.0
        self.omega1 = 0.0
        self.omega2 = 0.0
        self.omega3 = 0.0
        self.argp0 = 0.0
        self.argp1 = 0.0
        self.argp2 = 0.0
        self.argp3 = 0.0
        self.lonmean0 = 0.0
        self.lonmean1 = 0.0
        self.lonmean2 = 0.0
        self.lonmean3 = 0.0

    def readorbdata1(self):
        datacol = 0
        if self.abbr == 'me':
            datacol= 1
        elif self.abbr == 'v':
            datacol = 2
        elif self.abbr == 'e':
            datacol = 3
        elif self.abbr == 'ma':
            datacol = 4
        elif self.abbr == 'j':
            datacol = 5
        elif self.abbr == 's':
            datacol = 6
        elif self.abbr == 'u':
            datacol = 7
        elif self.abbr == 'n':
            datacol = 8
        elif self.abbr == 'p':
            datacol = 9
        else:
            print('Error: Invalid Planet Type')
            return

        fn = os.path.join(os.path.dirname(__file__), "data", "planetorbdata1.dat")
        orbdata = np.loadtxt(fname=fn,dtype='float',skiprows=2, usecols=datacol)
        self.a0 = orbdata[0]
        self.a1 = orbdata[1]
        self.a2 = orbdata[2]
        self.ecc0 = orbdata[3]
        self.ecc1 = orbdata[4]
        self.ecc2 = orbdata[5]
        self.ecc3 = orbdata[6]
        self.incl0 = orbdata[7]
        self.incl1 = orbdata[8]
        self.incl2 = orbdata[9]
        self.incl3 = orbdata[10]
        self.omega0 = orbdata[11]
        self.omega1 = orbdata[12]
        self.omega2 = orbdata[13]
        self.omega3 = orbdata[14]
        self.argp0 = orbdata[15]
        self.argp1 = orbdata[16]
        self.argp2 = orbdata[17]
        self.argp3 = orbdata[18]
        self.lonmean0 = orbdata[19]
        self.lonmean1 = orbdata[20]
        self.lonmean2 = orbdata[21]
        self.lonmean3 = orbdata[22]

    # All approx calculations: From Meeus (1991:202–204), the ecliptic orbital
    # elements of the planets, referenced to the standard equinox of J2000.

    # Approximates semi-major axis given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_a(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (self.a0 + self.a1 * ttdb + self.a2 * ttdb**2) * au

    # Approximates eccentricity given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_ecc(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (self.ecc0 + self.ecc1 * ttdb + self.a2 * ttdb**2 +
                    self.ecc3 * ttdb**3)

    # Approximates inclination given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_incl(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (self.incl0 + self.incl1 * ttdb + self.incl2 * ttdb**2 +
                    self.incl3 * ttdb**3) * deg2rad

    # Approximates omega (right ascension) given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_omega(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (self.omega0 + self.omega1 * ttdb + self.omega2 * ttdb**2 +
                    self.omega3 * ttdb**3) * deg2rad

    # Approximates longitude of periapsis given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_argp(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (self.argp0 + self.argp1 * ttdb + self.argp2 * ttdb**2 +
                    self.argp3 * ttdb**3) * deg2rad


    # Approximates mean longitude given number of centuries (ttdb) from year 2000 using
    # planetary ephemerides listed in appendix D of Vallado (2013, pgs 1046-1048)
    def approxcalc_lonmean(self, ttdb):
        if ttdb < -2.0 or ttdb > 0.5:
            print("Error: Only accurate between 1800 and 2050")
            return None
        else:
            return (self.lonmean0 + self.lonmean1 * ttdb + self.lonmean2 * ttdb**2 +
                    self.lonmean3 * ttdb**3) * deg2rad

    # To approximate heliocentric orbital elements (coes) referenced to the
    # standard equinox of J2000 (IAU-76/FK5)
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
            argp = self.approxcalc_argp(ttdb)
            lonmean = self.approxcalc_lonmean(ttdb)

            m = lonmean - argp
            argp = argp - omega
            _ , nu = smu.newtonm(ecc,m)

            return p, ecc, incl, omega, argp, nu

class Mercury(Planet):
    def __init__(self, orbdata_option = 1):
        self.abbr = 'me'
    #Can extend this if wanting to use other data sets
        if orbdata_option == 1:
            self.readorbdata1()
        else:
            self.readorbdata1()

class Venus(Planet):
    def __init__(self, orbdata_option = 1):
        self.abbr = 'v'
    #Can extend this if wanting to use other data sets
        if orbdata_option == 1:
            self.readorbdata1()
        else:
            self.readorbdata1()

class Earth(Planet):
    def __init__(self, orbdata_option = 1):
        self.abbr = 'e'
    # Can extend this if wanting to use other data sets
        if orbdata_option == 1:
            self.readorbdata1()
        else:
            self.readorbdata1()

class Mars(Planet):
    def __init__(self, orbdata_option = 1):
        self.abbr = 'ma'
    # Can extend this if wanting to use other data sets
        if orbdata_option == 1:
            self.readorbdata1()
        else:
            self.readorbdata1()

class Jupiter(Planet):
    def __init__(self, orbdata_option = 1):
        self.abbr = "j"
    # Can extend this if wanting to use other data sets
        if orbdata_option == 1:
            self.readorbdata1()
        else:
            self.readorbdata1()

class Saturn(Planet):
    def __init__(self, orbdata_option = 1):
        self.abbr = 's'
    # Can extend this if wanting to use other data sets
        if orbdata_option == 1:
            self.readorbdata1()
        else:
            self.readorbdata1()

class Uranus(Planet):
    def __init__(self, orbdata_option = 1):
        self.abbr = 'u'
    # Can extend this if wanting to use other data sets
        if orbdata_option == 1:
            self.readorbdata1()
        else:
            self.readorbdata1()

class Neptune(Planet):
    def __init__(self, orbdata_option = 1):
        self.abbr = 'n'
    # Can extend this if wanting to use other data sets
        if orbdata_option == 1:
            self.readorbdata1()
        else:
            self.readorbdata1()

class Pluto(Planet):
    def __init__(self, orbdata_option = 1):
        self.abbr = 'p'
    # Can extend this if wanting to use other data sets
        if orbdata_option == 1:
            self.readorbdata1()
        else:
            self.readorbdata1()

if __name__ == '__main__':
    planet = Mercury(1)