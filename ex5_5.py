#     -----------------------------------------------------------------
#
#                              Ex5_5
#
#  this file demonstrates example 5-5.
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

import spacetime_utils as stu
import space_conversions as sc
from space_constants import au, rad2deg, deg2rad
import spacemath_utils as smu
import orbit_utils as obu

year = 1994
mon = 5
day = 20
hr = 20
min = 0
sec = 0.0
dut1 = 0.0
dat = 0
timezone = 0


#, tcg, jdtcg, jdtcgfrac, tcb, jdtcb, jdtcbfrac
ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, \
jdtdb, jdtdbfrac = stu.convtime(year, mon, day, hr, min, sec, timezone, dut1,
                                dat)


# find coes for Jupiter
a = 5.202603191 + 1.913e-07 * ttdb

ecc = 0.04849485 + 0.000163244 * ttdb - 4.719e-07 * ttdb ** 2
incl = (1.30327 - 0.0019872 * ttdb + 3.318e-05 * ttdb ** 2
        + 9.2e-08 * ttdb ** 3) * deg2rad
omega = (100.464441 + 0.1766828 * ttdb + 0.000903877 * ttdb ** 2
         - 7.032e-06 * ttdb ** 3) * deg2rad
argp = (14.331309 + 0.2155525 * ttdb + 0.00072252 * ttdb ** 2
         - 4.59e-06 * ttdb ** 3) * deg2rad
lonmean = (34.351484 + 3034.9056746 * ttdb - 8.501e-05 * ttdb ** 2
           + 4e-09 * ttdb ** 3) * deg2rad
argp * rad2deg
lonmean * rad2deg
mean = lonmean - argp
argp = argp - omega
eccanom, nu = smu.newtonm(ecc, mean)
au = 149597870.7

p = a * (1.0 - ecc * ecc) * au

musun = 132712428000.0

# answer in km/s
r, v = sc.coe2rv(p, ecc, incl, omega, argp, nu, 0.0, 0.0, 0.0, musun)
# r in au
r = r / au
# v in au/day
v = (v / au) * 86400
print('\nEx 5.5: Jupiter rv\n')
print('r  %11.6f %11.6f %11.6f  AU \n' % (r[0], r[1], r[2]))
print('v  %11.6f %11.6f %11.6f  AU/day \n' % (v[0], v[1], v[2]))
print('coes %11.4f %11.4f %11.7f %11.5f %11.5f ' % (p / au, a, ecc, incl * rad2deg,
                                                    omega * rad2deg))
print('%11.5f %11.5f %11.5f\n' % (argp * rad2deg, nu * rad2deg, mean * rad2deg))
# Less precise eps (why values slightly differ from book to function)
eps = 23.440021 * deg2rad
reci = smu.rot1(r, - eps)
veci = smu.rot1(v, - eps)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci[0], reci[1], reci[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci[0], veci[1], veci[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci[0] * au, reci[1] * au,
                                             reci[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci[0] * au / 86400,
                                               veci[1] * au / 86400,
                                               veci[2] * au / 86400))
#mag(veci * au / 86400)

#Modularized the code above and made classes/function for planets

print("\nFunction Test: Mercury\n")
reci2, veci2 = obu.planetrv('me', jdut1 + jdut1frac)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci2[0], reci2[1], reci2[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci2[0], veci2[1], veci2[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci2[0] * au, reci2[1] * au,
                                             reci2[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci2[0] * au / 86400,
                                               veci2[1] * au / 86400,
                                               veci2[2] * au / 86400))

print("Function Test: Venus\n")
reci2, veci2 = obu.planetrv('v', jdut1 + jdut1frac)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci2[0], reci2[1], reci2[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci2[0], veci2[1], veci2[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci2[0] * au, reci2[1] * au,
                                             reci2[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci2[0] * au / 86400,
                                               veci2[1] * au / 86400,
                                               veci2[2] * au / 86400))

print("Function Test: Earth\n")
reci2, veci2 = obu.planetrv('e', jdut1 + jdut1frac)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci2[0], reci2[1], reci2[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci2[0], veci2[1], veci2[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci2[0] * au, reci2[1] * au,
                                             reci2[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci2[0] * au / 86400,
                                               veci2[1] * au / 86400,
                                               veci2[2] * au / 86400))

print("Function Test: Mars\n")
reci2, veci2 = obu.planetrv('ma', jdut1 + jdut1frac)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci2[0], reci2[1], reci2[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci2[0], veci2[1], veci2[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci2[0] * au, reci2[1] * au,
                                             reci2[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci2[0] * au / 86400,
                                               veci2[1] * au / 86400,
                                               veci2[2] * au / 86400))

print("Function Test: Jupiter\n")
reci2, veci2 = obu.planetrv('j', jdut1 + jdut1frac)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci2[0], reci2[1], reci2[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci2[0], veci2[1], veci2[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci2[0] * au, reci2[1] * au,
                                             reci2[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci2[0] * au / 86400,
                                               veci2[1] * au / 86400,
                                               veci2[2] * au / 86400))

print("Function Test: Saturn\n")
reci2, veci2 = obu.planetrv('s', jdut1 + jdut1frac)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci2[0], reci2[1], reci2[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci2[0], veci2[1], veci2[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci2[0] * au, reci2[1] * au,
                                             reci2[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci2[0] * au / 86400,
                                               veci2[1] * au / 86400,
                                               veci2[2] * au / 86400))

print("Function Test: Uranus\n")
reci2, veci2 = obu.planetrv('u', jdut1 + jdut1frac)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci2[0], reci2[1], reci2[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci2[0], veci2[1], veci2[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci2[0] * au, reci2[1] * au,
                                             reci2[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci2[0] * au / 86400,
                                               veci2[1] * au / 86400,
                                               veci2[2] * au / 86400))


print("Function Test: Neptune\n")
reci2, veci2 = obu.planetrv('n', jdut1 + jdut1frac)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci2[0], reci2[1], reci2[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci2[0], veci2[1], veci2[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci2[0] * au, reci2[1] * au,
                                             reci2[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci2[0] * au / 86400,
                                               veci2[1] * au / 86400,
                                               veci2[2] * au / 86400))


print("Function Test: Pluto\n")
reci2, veci2 = obu.planetrv('p', jdut1 + jdut1frac)
print('reci  %11.6f %11.6f %11.6f  AU \n' % (reci2[0], reci2[1], reci2[2]))
print('veci  %11.6f %11.6f %11.6f  AU/day \n' % (veci2[0], veci2[1], veci2[2]))
# now in km and km/s
print('reci  %11.1f %11.1f %11.1f  km \n' % (reci2[0] * au, reci2[1] * au,
                                             reci2[2] * au))
print('veci  %11.5f %11.5f %11.5f  km/s \n' % (veci2[0] * au / 86400,
                                               veci2[1] * au / 86400,
                                               veci2[2] * au / 86400))