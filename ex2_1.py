#     -----------------------------------------------------------------
#
#                              Ex2_1.py 
#                                 
#                 This file demonstrates Example 2-1.
#                     (4th edition - pgs. 66-67)
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
#   Find :
#     e0                 - eccentric anomaly            0.0 to 2pi rad
#     nu (extra)         - true anomaly                 0.0 to 2pi rad  


from spacemath_utils import newtonm
from space_constants import *

# --- newtonm - Finds eccentric and true anomaly given ecc and mean anomaly.
ecc = 0.4
m = 235.4 * deg2rad
e0, nu = newtonm(ecc, m) # iterative function

print('\n-------- ex 2-1 --------\n')
print('              m                e           nu           e0   \n')
print('newm  %14.8f %14.8f %14.8f %14.8f  rad \n' % (m, ecc, nu, e0))
print('newm  %14.8f %14.8f %14.8f %14.8f  deg \n'
      % (m * rad2deg, ecc, nu * rad2deg, e0 * rad2deg))
