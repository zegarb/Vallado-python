import math
# ------------------------------------------------------------------------------
#
#                           function constmath
#
#  this function sets constants for mathematical operations.
#
#  author        : david vallado                  719-573-2600    2 apr 2007
#
#  revisions
#
#  inputs        : description                    range / units
#    none
#
#  outputs       :
#    rad, twopi, halfpi
#    ft2m, mile2m, nm2m, mile2ft, mileph2kmph, nmph2kmph
#
#  locals        :
#                -
#
#  coupling      :
#    none.
#
# constmath
# ------------------------------------------------------------------------------


smalle6 = 1.0e-6
smalle8 = 1.0e-8
small = 1.0e-10
smalle12 = 1.0e-12
smalle18 = 1.0e-18


infinite = math.inf
undefined = 999999.1

# -------------------------  mathematical  --------------------
twopi = 2.0 * math.pi
halfpi = math.pi * 0.5
# -------------------------  conversions  ---------------------
ft2m = 0.3048
mile2m = 1609.344
nm2m = 1852
mile2ft = 5280
mileph2kmph = 0.44704
nmph2kmph = 0.5144444

rad2deg = 180.0 / math.pi
###added
deg2rad = math.pi / 180.0
auer = 149597870.0 /6378.1363
j2 = 0.001826267

# 1" to rad
arcsec2rad = math.pi / (3600.0 * 180)
rad2arcsec = 1 / arcsec2rad

#milli-arcseconds to rad
marcsec2rad = 1e-03 * arcsec2rad

#micro-arcseconds to rad
uarcsec2rad = 1e-06 * arcsec2rad

# Redundant with eccearthsqrd?
# eesqrd = 0.006694385000     # eccentricity of earth sqrd



# ------------------------------------------------------------------------------
#
#                           function constastro
#
#  this function sets constants for various astrodynamic operations.
#
#  author        : david vallado                  719-573-2600    2 apr 2007
#
#  revisions
#
#  inputs        : description                    range / units
#    none
#
#  outputs       :
#    re, flat, earthrot, mu
#    eccearth, eccearthsqrd
#    renm, reft, tusec, tumin, tuday, omegaearthradptu, omegaearthradpmin
#    velkmps, velftps, velradpmin
#    degpsec, radpday
#    speedoflight, au, earth2moon, moonradius, sunradius
#
#  locals        :
#                -
#
#  coupling      :
#    none.
#
# ------------------------------------------------------------------------------


# -----------------------  physical constants  ----------------
# EGM-08 constants used here
re = 6378.1363         # km
flat = 1.0/298.257223563
earthrot = 7.292115e-5     # rad /s  old 7.29211514670698e-05
mu = 398600.4415      # km3 /s2
mum = 3.986004415e14   # m3 /s2

# derived constants from the base values
eccearth = math.sqrt(2.0*flat - flat**2)
eccearthsqrd = eccearth**2

renm = re / nm2m
reft = re * 1000.0 / ft2m

tusec = math.sqrt(re**3/mu)
tumin = tusec / 60.0
tuday = tusec / 86400.0
tudaysid = tusec / 86164.090524

omegaearthradptu = earthrot * tusec
omegaearthradpmin = earthrot * 60.0

velkmps = math.sqrt(mu / re)
velftps = velkmps * 1000.0 / ft2m
velradpmin = velkmps * 60.0 / re
# for afspc
# velkmps1 = velradpmin*6378.135/60.0   7.90537051051763
# mu1 = velkmps*velkmps*6378.135        3.986003602567418e+005
degpsec = (180.0 / math.pi) / tusec
radpday = 2.0 * math.pi * 1.002737909350795

speedoflight = 299792.458 # km/s
au = 149597870.7      # km
earth2moon = 384400.0 # km
moonradius = 1738.0 # km
sunradius = 696000.0 # km

masssun = 1.9891e30
massearth = 5.9742e24
massmoon = 7.3483e22

musun = 1.32712428e11
mumoon = 4902.799

#Umbra and penumbra angle using mean sun-earth distance
angumbearth = 0.004695073367587701 #rad
angpenearth = 0.004609804797371252 #rad



# ------------------------------------------------------------------------------
#
#                           function sethelp
#
#  this function sets help flags to control intermediate output during
#  debugging
#
#  author        : david vallado                  719-573-2600   16 jul 2002
#
#  revisions
#
#  inputs        : description                    range / units
#    none
#
#  outputs       :
#    iauhelp
#
#
# sethelp;
# ------------------------------------------------------------------------------

class sethelp():
  iauhelp = False;
  iaupnhelp = True;
  show = False;

