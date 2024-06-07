#     -----------------------------------------------------------------
#
#                              Ex1_1.py
#
#                 This file demonstrates Example 1-1.
#                     (4th edition - pgs. 32-33)
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
#   Given :            Description                   Range / Units
#      periodsid       - one orbit period (sidereal)   sec
#
#   Constants :
#      mu              - gravitional parameter
#
#   Find :
#      a               - semimajor axis                km^3/s^2

from space_constants import *

# A geosynchronous orbit (approximate Earth revolution time)
# Orbit Period = 24 sidereal hours ~ 86,164.09 secs (NOT 24 solar hours)
periodsid = 86164.090524

print('\n-------- ex 1-1 --------\n')

# Semimajor axis, a
a = (mu * (periodsid / twopi) ** 2) ** (1.0 / 3.0)

print('Results:\n')
# ER = Earth Radius (Earth's equatorial radius) (Sec. 3.8, pgs 237-238)
print('a        %16.8f km   %18.10f er \n' % (a, a / re))

print('Provided:\n')
print('re      %16.8f km           mu   %18.10f km^3/s^2\n' % (re, mu))

print('periodsid  %16.10f sec\
      periodsid alt   %16.10f sec\n' % (periodsid, (86400 / 1.002737909350795)))

# TU = Time Unit (Sec. 3.8, pg. 238)
print('tuday      %16.14f           1/tuday %18.12f  \
      \n' % (tuday, 1.0 / tuday))
print('tudaysid   %16.14f           1/tuday %18.12f  \
      \n' % (tudaysid, 1.0 / tudaysid))

