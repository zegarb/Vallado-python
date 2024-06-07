#     -----------------------------------------------------------------
#
#                              Ex2_3.py 
#                                 
#                 This file demonstrates Example 2-3.
#                     (4th edition - pgs. 71-72)
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
#     ecc                - eccentricity                 0.0 to 
#     m                  - mean anomaly                 0.0 to 2pi rad
#
#
#   Find : 
#     h                  - hyperbolic anomaly           0.0 to 2pi rad
#     nu (extra)         - true anomaly                 0.0 to 2pi rad  


from spacemath_utils import newtonm
from space_constants import *

# --- newtonm - Finds eccentric and true anomaly given ecc and mean anomaly.
ecc = 2.4
m = 235.4 * deg2rad

# Function performs iteration as shown in table 2-3 (pg 72)
h, nu = newtonm(ecc, m)

print('\n-------- ex 2-3 --------\n')
print('            m                e              nu             h   \n')
print('newm  %14.8f %14.8f %14.8f %14.8f  rad \n'
      % (m, ecc, nu, h))
print('newm  %14.8f %14.8f %14.8f %14.8f  deg \n'
      % (m * rad2deg, ecc, nu * rad2deg, h * rad2deg))

