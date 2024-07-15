import math

#--------------------------------------------------------------
#--------------------------constmath --------------------------
#--------------------------------------------------------------

smalle6 = 1.0e-6
smalle8 = 1.0e-8
small = 1.0e-10
smalle12 = 1.0e-12
smalle18 = 1.0e-18

infinite = math.inf
undefined = float('nan')

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
j2 = 0.001826267 # EGM-08

# 1" to rad
arcsec2rad = math.pi / (3600.0 * 180)
rad2arcsec = 1 / arcsec2rad

#milli-arcseconds to rad
marcsec2rad = 1e-03 * arcsec2rad

#micro-arcseconds to rad
uarcsec2rad = 1e-06 * arcsec2rad

# Redundant with eccearthsqrd?
# eesqrd = 0.006694385000     # eccentricity of earth sqrd



#--------------------------------------------------------------
#---------------------constastro-------------------------------
#--------------------------------------------------------------

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

tusec = math.sqrt(re**3/mu) # EGM-08
tumin = tusec / 60.0 # EGM-08
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


# -----------------------------------------------------------------------------
#
#                           function getgravc
#
#  this function gets Earth gravitional constants for a propagator.
#  note that mu is identified to facilitiate comparisons with newer models.
#
#  author        : david vallado                  719-573-2600   21 jul 2006
#
#  inputs        :
#    whichconst  - which set of constants to use  wgs721, wgs72, wgs84, egm08
#
#  outputs       :
#    gravc dictionary
#                 mu          - earth gravitational parameter
#                 re          - radius of the earth in km
#                 j2, j3, j4  - un-normalized zonal harmonic values
#                 j3oj2       - j3 divided by j2
#                 xke         - reciprocal of tumin
#                 tumin       - minutes in one time unit
#
#  locals        :
#
#  coupling      :
#
#  references    :
#    norad spacetrack report #3
#    vallado, crawford, hujsak, kelso  2006
# [tumin, mu, re, xke, j2, j3, j4, j3oj2] = getgravc(whichconst)
#  --------------------------------------------------------------------------- */

def getgravc(whichconst=None):
  """this function gets Earth gravitional constants for a propagator.
  note that mu is identified to facilitiate comparisons with newer models.

  Parameters
  ----------
  whichconst : str
      which set of constants to use: 'wgs721', 'wgs72', 'wgs84', 'egm08'

  Returns
  -------
  gravc : Dict[str, float]
      gravitional constants dictionary
              mu - earth gravitational parameter: km^3/s^2
              re - radius of the earth: km
              j2, j3, j4 - un-normalized zonal harmonic values
              j3oj2 - j3 divided by j2
              xke - reciprocal of tumin
              tumin - reciprocal of xke (minutes in one time unit)


  """
  #global tumin mu re xke j2 j3 j4 j3oj2
  gravc = {}
  if 'wgs721' == whichconst:
      # -- wgs-72 low precision str#3 constants --

    mu = 398600.79964
    re = 6378.135 # radius of earth in km
    j2 = 0.001082616
    j3 = - 2.53881e-06
    j4 = - 1.65597e-06
    j3oj2 = j3 / j2
    xke = 0.0743669161
    tumin = 1.0 / xke

    gravc = {
      'mu': mu,
      're': re,
      'j2': j2,
      'j3': j3,
      'j4': j4,
      'j3oj2': j3oj2,
      'xke': xke,
      'tumin': tumin
    }

  elif 'wgs72' == whichconst:
    # ------------ wgs-72 constants ------------
    mu = 398600.8
    re = 6378.135
    j2 = 0.001082616
    j3 = - 2.53881e-06
    j4 = - 1.65597e-06
    j3oj2 = (j3 / j2)
    xke = 60.0 / math.sqrt(re**3 / mu)
    tumin = 1.0 / xke

    gravc = {
      'mu': mu,
      're': re,
      'j2': j2,
      'j3': j3,
      'j4': j4,
      'j3oj2': j3oj2,
      'xke': xke,
      'tumin': tumin
    }

  elif 'wgs84' == whichconst:
    # ------------ wgs-84 constants ------------

    mu = 398600.5
    re = 6378.137
    j2 = 0.00108262998905
    j3 = - 2.53215306e-06
    j4 = - 1.61098761e-06
    j3oj2 = j3 / j2
    xke = 60.0 / math.sqrt(re**3 / mu)
    tumin = 1.0 / xke

    gravc = {
      'mu': mu,
      're': re,
      'j2': j2,
      'j3': j3,
      'j4': j4,
      'j3oj2': j3oj2,
      'xke': xke,
      'tumin': tumin
    }

  elif 'egm08' == whichconst:
    # ------------ egm-08 constants ------------
    mu = 398600.4415
    re = 6378.1363
    j2 = 0.00108262617
    j3 = - 2.53241052e-06
    j4 = - 1.6198976e-06
    j3oj2 = j3 / j2
    xke = 60.0 / math.sqrt(re**3 / mu)
    tumin = 1.0 / xke

    gravc = {
      'mu': mu,
      're': re,
      'j2': j2,
      'j3': j3,
      'j4': j4,
      'j3oj2': j3oj2,
      'xke': xke,
      'tumin': tumin
    }

  else:
      print('unknown gravity option (%d)\n' % (whichconst))


  return gravc

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
  """this function sets help flags to control intermediate output during
  debugging

  show : bool
      general view print flag
  iauhelp : bool
      flag used in iau06pnb
  iauhelp : bool
      flag used in iau functions

  """
  show = False;
  iauhelp = False;
  iaupnhelp = True;


