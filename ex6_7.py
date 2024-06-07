#     -----------------------------------------------------------------
#
#                              ex6_7.m
#
#  this file demonstrates example 6-7.
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


import numpy as np
import orbit_utils as obu
from space_constants import au, rad2deg, deg2rad, re, velkmps, tumin, mu
import space_conversions as sc
import spacemath_utils as smu


print('-------------------- problem ex 6-7 \n')
rinit = (re + 191.0) / re
rfinal = (re + 35780.0) / re
print('from radius %11.5f to %11.5f ---------------------- \n'
      % (rinit * re, rfinal * re))
einit = 0.0
efinal = 0.0
iinit = 28.5 * deg2rad
ifinal = 0.0 * deg2rad
print('from inclination %11.5f to %11.5f \n' % (iinit * rad2deg, ifinal * rad2deg))
deltai = - ifinal + iinit
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2 = obu.combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai)
print('combined maneuver \n' % ())
print(' deltava  %11.7f %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltai1  %11.7f  %11.7f  \n' % (deltai1 * rad2deg, deltai2 * rad2deg))
print(' gam1  %11.7f \n' % (gam1 * rad2deg))
print(' gam2  %11.7f \n' % (gam2 * rad2deg))
#End of Book Example


atran = (rinit + rfinal) * re / 2.0
etran = re * (rfinal - rinit) / (2 * atran)
ptran = atran * (1.0 - etran ** 2)
#                   a      e       i               node         argp    nu      ci arglat ce truelon ee lonper
r1, v1 = sc.coe2rv(rinit * re, 0.0, 28.5 * deg2rad, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0 * deg2rad, 0.0, 0.0)
r2, v2 = sc.coe2rv(ptran, etran, 28.5 * deg2rad - deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0, 0.0, 0.0)
r3, v3 = sc.coe2rv(ptran, etran, 28.5 * deg2rad - deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  180.0 * deg2rad, 0.0, 0.0, 0.0)
r4, v4 = sc.coe2rv(rfinal * re, 0.0, ifinal, 0.0 * deg2rad, 0.0 * deg2rad, 0.0 * deg2rad,
                  0.0, 225.0 * deg2rad, 0.0)
vx = v2 - v1
print(' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n' % (vx[0], vx[1], vx[2], smu.mag(vx)))
vy = v4 - v3
print(' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n' % (vy[0], vy[1], vy[2], smu.mag(vy)))
#pause
# other tests
# ===================================================================
rinit = (re + 222.0) / re
rfinal = (re + 35786.0) / re
einit = 0.0
efinal = 0.0
iinit = 28.5 * deg2rad
ifinal = 89.0 * deg2rad
deltai = ifinal - iinit
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
print('from radius %11.5f to %11.5f GEO Polar ---------------------- \n' % (rinit * re, rfinal * re))
print('from inclination %11.5f to %11.5f \n' % (iinit * rad2deg, ifinal * rad2deg))
deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2 = \
    obu.combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai)
print('combined maneuver \n')
print(' deltava  %11.7f %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltai1  %11.7f  %11.7f  \n' % (deltai1 * rad2deg, deltai2 * rad2deg))
print(' gam1  %11.7f \n' % (gam1 * rad2deg))
print(' gam2  %11.7f \n' % (gam2 * rad2deg))
atran = (rinit + rfinal) * re / 2.0
etran = re * (rfinal - rinit) / (2 * atran)
ptran = atran * (1.0 - etran ** 2)
#                   a      e       i               node         argp    nu      ci arglat ce truelon ee lonper
r1, v1 = sc.coe2rv(rinit * re, 0.0, 28.5 * deg2rad, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0 * deg2rad, 0.0, 0.0)
r2, v2 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0, 0.0, 0.0)
r3, v3 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  180.0 * deg2rad, 0.0, 0.0, 0.0)
r4, v4 = sc.coe2rv(rfinal * re, 0.0, ifinal, 45.0 * deg2rad, 0.0 * deg2rad, 0.0 * deg2rad,
                  180.0, 0.0, 0.0)
vx= v2 - v1
print(' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vx[0], vx[1], vx[2], smu.mag(vx)))
vy= v4 - v3
print(' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vy[0], vy[1], vy[2], smu.mag(vy)))

# ===================================================================
rinit = (re + 222.0) / re
rfinal = (re + 35786.0) / re
einit = 0.0
efinal = 0.0
iinit = 28.5 * deg2rad
ifinal = 0.0 * deg2rad
deltai = ifinal - iinit
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
print('from radius %11.5f to %11.5f GEO ---------------------- \n'
      % (rinit * re, rfinal * re))
print('from inclination %11.5f to %11.5f \n' % (iinit * rad2deg, ifinal * rad2deg))
deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2 \
    = obu.combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai)
print('combined maneuver \n')
print(' deltava  %11.7f %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltai1  %11.7f  %11.7f  \n' % (deltai1 * rad2deg, deltai2 * rad2deg))
print(' gam1  %11.7f \n' % (gam1 * rad2deg))
print(' gam2  %11.7f \n' % (gam2 * rad2deg))
atran = (rinit + rfinal) * re / 2.0
etran = re * (rfinal - rinit) / (2 * atran)
ptran = atran * (1.0 - etran ** 2)
#                   a      e       i               node         argp    nu       arglat  truelon lonper
r1, v1 = sc.coe2rv(rinit * re, 0.0, 28.5 * deg2rad, 45.0 * deg2rad, 0.0 * deg2rad,
                   0.0 * deg2rad, 0.0 * deg2rad, 0.0, 0.0)
r2, v2 = sc.coe2rv(ptran, etran, 28.5 * deg2rad - deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                   0.0 * deg2rad, 0.0, 0.0, 0.0)
r3, v3 = sc.coe2rv(ptran, etran, 28.5 * deg2rad - deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                   180.0 * deg2rad, 0.0, 0.0, 0.0)
r4, v4 = sc.coe2rv(rfinal * re, 0.0, ifinal, 0.0 * deg2rad, 0.0 * deg2rad, 0.0 * deg2rad,
                   0.0, 225.0 * deg2rad, 0.0)
vx = v2 - v1
print(' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vx[0], vx[1], vx[2], smu.mag(vx)))
vy= v4-v3
print(' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vy[0], vy[1], vy[2], smu.mag(vy)))

# ===================================================================
rinit = (re + 222.0) / re
rfinal = (re + 15392.0) / re
einit = 0.0
efinal = 0.64887
iinit = 28.5 * deg2rad
ifinal = 63.4 * deg2rad
deltai = ifinal - iinit
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
print('from radius %11.5f to %11.5f Mod HEO ---------------------- \n'
      % (rinit * re, rfinal * re))
print('from inclination %11.5f to %11.5f \n' % (iinit * rad2deg, ifinal * rad2deg))
deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2 = \
    obu.combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai)
print('combined maneuver \n' % ())
print(' deltava  %11.7f %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltai1  %11.7f  %11.7f  \n' % (deltai1 * rad2deg, deltai2 * rad2deg))
print(' gam1  %11.7f \n' % (gam1 * rad2deg))
print(' gam2  %11.7f \n' % (gam2 * rad2deg))
atran = (rinit + rfinal) * re / 2.0
etran = re * (rfinal - rinit) / (2 * atran)
ptran = atran * (1.0 - etran ** 2)
#                   a      e       i               node         argp    nu      ci arglat ce truelon ee lonper
r1, v1 = sc.coe2rv(rinit * re, 0.0, 28.5 * deg2rad, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0 * deg2rad, 0.0, 0.0)
r2, v2 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0, 0.0, 0.0)
r3, v3 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  180.0 * deg2rad, 0.0, 0.0, 0.0)
r4, v4 = sc.coe2rv(rfinal * re, efinal, ifinal, 45.0 * deg2rad, 270.0 * deg2rad,
                  180.0 * deg2rad, 0.0, 0.0, 0.0)
vx = v2 - v1
print(' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vx[0], vx[1], vx[2], smu.mag(vx)))
vy = v4 - v3
print(' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vy[0], vy[1], vy[2], smu.mag(vy)))

# ===================================================================
rinit = (re + 222.0) / re
rfinal = (re + 35786.0) / re
einit = 0.0
efinal = 0.0
iinit = 28.5 * deg2rad
ifinal = 49.3 * deg2rad
deltai = ifinal - iinit
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
print('from radius %11.5f to %11.5f Quasi-Zen ---------------------- \n'
      % (rinit * re, rfinal * re))
print('from inclination %11.5f to %11.5f \n' % (iinit * rad2deg, ifinal * rad2deg))
deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2 = \
    obu.combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai)
print('combined maneuver \n')
print(' deltava  %11.7f %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltai1  %11.7f  %11.7f  \n' % (deltai1 * rad2deg, deltai2 * rad2deg))
print(' gam1  %11.7f \n' % (gam1 * rad2deg))
print(' gam2  %11.7f \n' % (gam2 * rad2deg))
atran = (rinit + rfinal) * re / 2.0
etran = re * (rfinal - rinit) / (2 * atran)
ptran = atran * (1.0 - etran ** 2)
#                   a      e       i               node         argp    nu       arglat  truelon lonper
r1, v1 = sc.coe2rv(rinit * re, 0.0, 28.5 * deg2rad, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0 * deg2rad, 0.0, 0.0)
r2, v2 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0, 0.0, 0.0)
r3, v3 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  180.0 * deg2rad, 0.0, 0.0, 0.0)
r4, v4 = sc.coe2rv(rfinal * re, 0.0, ifinal, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 180.0, 0.0, 0.0)
vx = v2 - v1
print(' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vx[0], vx[1], vx[2], smu.mag(vx)))
vy = v4 - v3
print(' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vy[0], vy[1], vy[2], smu.mag(vy)))

# ===================================================================
rinit = (re + 222.0) / re
rfinal = (re + 37028.0) / re
einit = 0.0
efinal = 0.0
iinit = 28.5 * deg2rad
ifinal = 0.0 * deg2rad
deltai = ifinal - iinit
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
print('from radius %11.5f to %11.5f Super GEO ---------------------- \n'
      % (rinit * re, rfinal * re))
print('from inclination %11.5f to %11.5f \n' % (iinit * rad2deg, ifinal * rad2deg))
deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2 = \
    obu.combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai)
print('combined maneuver \n')
print(' deltava  %11.7f %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltai1  %11.7f  %11.7f  \n' % (deltai1 * rad2deg, deltai2 * rad2deg))
print(' gam1  %11.7f \n' % (gam1 * rad2deg))
print(' gam2  %11.7f \n' % (gam2 * rad2deg))
atran = (rinit + rfinal) * re / 2.0
etran = re * (rfinal - rinit) / (2 * atran)
ptran = atran * (1.0 - etran ** 2)
#                   a      e       i               node         argp    nu       arglat  truelon lonper
r1, v1 = sc.coe2rv(rinit * re, 0.0, 28.5 * deg2rad, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0 * deg2rad, 0.0, 0.0)
r2, v2 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0, 0.0, 0.0)
r3, v3 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  180.0 * deg2rad, 0.0, 0.0, 0.0)
r4, v4 = sc.coe2rv(rfinal * re, 0.0, ifinal, 0.0 * deg2rad, 0.0 * deg2rad, 0.0 * deg2rad,
                  0.0, 225.0 * deg2rad, 0.0)
vx = v2 - v1
print(' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vx[0], vx[1], vx[2], smu.mag(vx)))
vy = v4 - v3
print(' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vy[0], vy[1], vy[2], smu.mag(vy)))

# ===================================================================
rinit = (re + 222.0) / re
rfinal = (re + 35786.0) / re
einit = 0.0
efinal = 0.268
iinit = 28.5 * deg2rad
ifinal = 63.4 * deg2rad
deltai = ifinal - iinit
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
print('from radius %11.5f to %11.5f Sirius ---------------------- \n'
      % (rinit * re, rfinal * re))
print('from inclination %11.5f to %11.5f \n' % (iinit * rad2deg, ifinal * rad2deg))
deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2 = \
    obu.combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai)
print('combined maneuver \n')
print(' deltava  %11.7f %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltai1  %11.7f  %11.7f  \n' % (deltai1 * rad2deg, deltai2 * rad2deg))
print(' gam1  %11.7f \n' % (gam1 * rad2deg))
print(' gam2  %11.7f \n' % (gam2 * rad2deg))
atran = (rinit + rfinal) * re / 2.0
etran = re * (rfinal - rinit) / (2 * atran)
ptran = atran * (1.0 - etran ** 2)
#                   a      e       i               node         argp    nu       arglat  truelon lonper
r1, v1 = sc.coe2rv(rinit * re, 0.0, 28.5 * deg2rad, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0 * deg2rad, 0.0, 0.0)
r2, v2 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  0.0 * deg2rad, 0.0, 0.0, 0.0)
r3, v3 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                  180.0 * deg2rad, 0.0, 0.0, 0.0)
r4, v4 = sc.coe2rv(rfinal * re, efinal, ifinal, 45.0 * deg2rad, 270.0 * deg2rad,
                  180.0 * deg2rad, 0.0, 0.0, 0.0)
vx = v2 - v1
print(' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vx[0], vx[1], vx[2], smu.mag(vx)))
vy = v4 - v3
print(' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n'
      % (vy[0], vy[1], vy[2], smu.mag(vy)))

# ===================================================================
rinit = (re + 222.0) / re
rfinal = (re + 35786.0) / re
einit = 0.0
efinal = 0.700381
iinit = 28.5 * deg2rad
ifinal = 63.4 * deg2rad
deltai = ifinal - iinit
nuinit = 0.0 * deg2rad
nufinal = 180.0 * deg2rad
print('from radius %11.5f to %11.5f Tundra-like ---------------------- \n'
      % (rinit * re, rfinal * re))
print('from inclination %11.5f to %11.5f \n' % (iinit * rad2deg, ifinal * rad2deg))
deltai1, deltai2, deltava, deltavb, dttu, gam1, gam2 = \
    obu.combined(rinit, rfinal, einit, efinal, nuinit, nufinal, deltai)
print('combined maneuver \n')
print(' deltava  %11.7f %11.7f km/s \n' % (deltava, deltava * velkmps))
print(' deltavb  %11.7f %11.7f km/s \n' % (deltavb, deltavb * velkmps))
print(' deltai1  %11.7f  %11.7f  \n' % (deltai1 * rad2deg, deltai2 * rad2deg))
print(' gam1  %11.7f \n' % (gam1 * rad2deg))
print(' gam2  %11.7f \n' % (gam2 * rad2deg))
atran = (rinit + rfinal) * re / 2.0
etran = re * (rfinal - rinit) / (2 * atran)
ptran = atran * (1.0 - etran ** 2)
#                   a      e       i               node         argp    nu       arglat  truelon lonper
r1, v1 = sc.coe2rv(rinit * re, 0.0, 28.5 * deg2rad, 45.0 * deg2rad, 0.0 * deg2rad,
                   0.0 * deg2rad, 0.0 * deg2rad, 0.0, 0.0)
r2, v2 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                   0.0 * deg2rad, 0.0, 0.0, 0.0)
r3, v3 = sc.coe2rv(ptran, etran, 28.5 * deg2rad + deltai1, 45.0 * deg2rad, 0.0 * deg2rad,
                   180.0 * deg2rad, 0.0, 0.0, 0.0)
r4, v4 = sc.coe2rv(rfinal * re, efinal, ifinal, 45.0 * deg2rad, 270.0 * deg2rad,
                   180.0 * deg2rad, 0.0, 0.0, 0.0)
vx = v2 - v1
print(' manv1 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n' % (vx[0], vx[1], vx[2], smu.mag(vx)))
vy = v4 - v3
print(' manv2 %11.7f %11.7f %11.7f km/s in icrf %11.7f \n' % (vy[0], vy[1], vy[2], smu.mag(vy)))

