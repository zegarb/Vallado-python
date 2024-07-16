import numpy as np

from space_conversions import coe2rv
from space_constants import *
from orbit_utils import kepler, keplercoe



#     -----------------------------------------------------------------
#
#                              Ex2_4
#
#  this file demonstrates example 2-4.
#
#                          companion code for
#             fundamentals of astrodynamics and applications
#                                 2013
#                            by david vallado
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



print('\n-------- kepler  ex 2-4, pg 102 --------- \n')

# initial vectors in km and km/s
ro = np.array([ 1131.340,  -2282.343,  6672.423])
vo = np.array([ -5.64305,  4.30333,  2.42879 ])
print('input: \n')
print('ro %16.8f %16.8f %16.8f km \n' % (ro[0], ro[1], ro[2]))
print('vo %16.8f %16.8f %16.8f km/s \n' % (vo[0], vo[1], vo[2]))

# convert 40 minutes to seconds
dtsec = 40.0*60.0
print('dt %16.8f sec \n' % dtsec)
print('intermediate values: \n')

r1, v1, errk =  kepler(ro, vo, dtsec)

# answer in km and km/s
print('kepler errk: %s\n' % errk)
print('r1 %16.8f %16.8f %16.8f er \n' % (r1[0]/re, r1[1]/re, r1[2]/re))
print('r1 %16.8f %16.8f %16.8f km \n' % (r1[0], r1[1], r1[2]))
print('v1 %16.8f %16.8f %16.8f er/tu \n' % (v1[0]/velkmps, v1[1]/velkmps,
                                            v1[2]/velkmps))
print('v1 %16.8f %16.8f %16.8f km/s \n\n\n' % (v1[0], v1[1], v1[2]))

r1, v1 = keplercoe(ro, vo, dtsec)

print('keplercoe test:')
print('r1 %16.8f %16.8f %16.8f er \n' % (r1[0]/re, r1[1]/re, r1[2]/re))
print('r1 %16.8f %16.8f %16.8f km \n' % (r1[0], r1[1], r1[2]))
print('v1 %16.8f %16.8f %16.8f er/tu \n' % (v1[0]/velkmps, v1[1]/velkmps,
                                            v1[2]/velkmps))
print('v1 %16.8f %16.8f %16.8f km/s \n' % (v1[0], v1[1], v1[2]))


# alt tests
# initial coes with more than one period = 6281.815597 sec
ro, vo = coe2rv(7358.39, 0.0, 28.5 * deg2rad, 0.0 * deg2rad, 30.0 * deg2rad,
                0.0 * deg2rad, 0.0, 0.0, 0.0)
print('input: \n')
print('ro %16.8f %16.8f %16.8f km \n'% (ro[0], ro[1], ro[2]))
print('vo %16.8f %16.8f %16.8f km/s \n'% (vo[0], vo[1], vo[2]))

# convert 40 minutes to seconds
dtsec = 4000.0*60.0
dtsec = 1.291007302335531e+03
dtsec = 6281.815597
print('dt %16.8f sec \n' % dtsec)
print('intermediate values: \n')

r1, v1, errk =  kepler(ro, vo, dtsec)

# answer in km and km/s
print('kepler errk: %s\n' % errk)
print('r1 %16.8f %16.8f %16.8f er \n' % (r1[0]/re, r1[1]/re, r1[2]/re))
print('r1 %16.8f %16.8f %16.8f km \n' % (r1[0], r1[1], r1[2]))
print('v1 %16.8f %16.8f %16.8f er/tu \n' % (v1[0]/velkmps, v1[1]/velkmps,
                                            v1[2]/velkmps))
print('v1 %16.8f %16.8f %16.8f km/s \n' % (v1[0], v1[1], v1[2]))


ro=np.array([-3244.01178958993, 5561.5015207476, 3181.63137126354])
ro = ro.T
vo=np.array([-0.311911476329513, 3.55766787343696, -6.53796978233233])
vo = vo.T
dtsec = 240.0
print('dt %16.8f sec \n' % dtsec)
print('intermediate values: \n')

r1, v1, errk =  kepler(ro, vo, dtsec)

# answer in km and km/s
print('kepler errk: %s\n' % errk)
print('r1 %16.8f %16.8f %16.8f er \n' % (r1[0]/re, r1[1]/re, r1[2]/re))
print('r1 %16.8f %16.8f %16.8f km \n' % (r1[0], r1[1], r1[2]))
print('v1 %16.8f %16.8f %16.8f er/tu \n' % (v1[0]/velkmps, v1[1]/velkmps,
                                            v1[2]/velkmps))
print('v1 %16.8f %16.8f %16.8f km/s \n' % (v1[0], v1[1], v1[2]))

r1, v1 = keplercoe(ro, vo, dtsec)

print('keplercoe test:')
print('r1 %16.8f %16.8f %16.8f er \n' % (r1[0]/re, r1[1]/re, r1[2]/re))
print('r1 %16.8f %16.8f %16.8f km \n' % (r1[0], r1[1], r1[2]))
print('v1 %16.8f %16.8f %16.8f er/tu \n' % (v1[0]/velkmps, v1[1]/velkmps,
                                            v1[2]/velkmps))
print('v1 %16.8f %16.8f %16.8f km/s \n' % (v1[0], v1[1], v1[2]))

