#     -----------------------------------------------------------------
#
#                              Ex2_2
#
#                 This file demonstrates Example 2-2.
#                     (4th edition - pgs. 69-70)
#
#                          Companion code for
#             Fundamentals of Astrodynamics and Applications
#                                (2013)
#                           By David Vallado
#
#     (h)               email davallado@gmail.com
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#            16 feb 19  david vallado
#                         update for new constants
#  changes :
#            13 feb 07  david vallado
#                         original baseline
#
#     *****************************************************************
#
#   Given :               Description                   Range / Units
#     p                  - semiparameter                km
#     dt                 - change in time               sec
#     ecc                - eccentricity                 0.0 to
#
#   Intermediate :
#     np                 - parabolic mean motion        rad/sec
#
#   Find :
#     b                  - parabolic anomaly            0.0 to 2pi rad


from spacemath_utils import newtonm, newtone, newtonnu, cot
from space_constants import *

p = 25512
dt = 53.7874 * 60 # convert min to sec
ecc = 1.0

# First find mean motion
np = 2.0 * math.sqrt(mu/(p**3))

# Intermediate subsitions/relations (clever math tricks)
cot2s = 1.5*np*dt # Time equation
s = 0.5 * math.atan(1.0/cot2s)
w = math.atan((math.tan(s)) ** (1/3))

# Result
b = 2.0 * cot(2.0*w)
print('\n-------- ex 2-2 --------\n')
print('np: %11.7f rad/s    cot2s: %11.7f    s: %11.6f deg' \
      '       w: %11.6f  deg\n\nb: %11.6f rad \n\n' \
      %(np, cot2s, s * rad2deg, w * rad2deg, b))

# Function Check (All math above is accounted for within function)
# --- newtonm - Finds eccentric and true anomaly given ecc and mean anomaly.
print('Ex 2-2 Equivalent Function Test:\n')
b, nu = newtonm(ecc,np*dt)
print('               m               e           nu           b   \n')
print('newm  %14.8f %14.8f %14.8f %14.8f  rad \n'
      % ((np*dt), ecc, nu, b))
print('newm  %14.8f %14.8f %14.8f %14.8f  deg \n'
      % ((np*dt) * rad2deg, ecc, nu * rad2deg, b * rad2deg ))


# -------------------- Related Tests --------------------
print('Related tests:\n')

# newtone - Finds true and mean anomaly given ecc and eccentric anomaly.
ecc= 0.4
e0 = 334.566986 * deg2rad
m, nu = newtone(ecc, e0)
print('               e0              e            m           nu   \n')
print('newe  %14.8f %14.8f %14.8f %14.8f  deg \n'
      % (e0 * rad2deg, ecc, m * rad2deg, nu * rad2deg))

# --- newtonm - Finds eccentric and true anomaly given ecc and mean anomaly.
ecc= 0.34
m = 235.4 * deg2rad
e0, nu = newtonm(ecc, m)
print('               m               e           nu           e0   \n')
print('newm  %14.8f %14.8f %14.8f %14.8f  deg \n'
      % (m * rad2deg, ecc, nu * rad2deg, e0 * rad2deg ))

# --- newtonnu - Finds eccentric and mean anomaly given ecc and true anomaly.
ecc= 0.34
nu = 134.567001 * deg2rad
e0, m = newtonnu(ecc, nu)
print('               nu              e           e0            m      \n')
print('newnu %14.8f %14.8f %14.8f %14.8f  deg \n'
      % (nu * rad2deg, ecc, e0 * rad2deg, m * rad2deg ))



