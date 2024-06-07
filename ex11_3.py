#     -----------------------------------------------------------------
#
#                              Ex11_3.m
#
#  this file demonstrates example 11-3.
#
#                          companion code for
#             fundamentals of astrodynamics and applications
#                                 2007
#                            by david vallado
#
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#            30 dec 12  david vallado
#                         original
#  changes :
#            30 dec 12 david vallado
#                         original baseline
#
#     *****************************************************************

import numpy as np
from space_constants import *
from spacemath_utils import mag
from spacetime_utils import lstime


print("--------- Example 11-3 ----------------------")

j2 = 0.00108263
# --------  repeat gt calculations
alt = 160.0

ecc = None
incl = 45.0 * deg2rad
latgd = 41.52 * deg2rad
lon = 12.37 * deg2rad
salt = 0.152

etafov = 20.0 * deg2rad
etactr = 15.0 * deg2rad
sinlat = np.sin(latgd)
# ------  find rdel and rk components of site vector  ---------

# In book rdel = 4782.104 and rk = 4207.153; slightly off from
# answer here; doesn't effect future math
cearth = re / np.sqrt(1.0 - (eccearthsqrd * sinlat * sinlat))
rdel = (cearth + salt) * np.cos(latgd)
rk = ((1.0 - eccearthsqrd) * cearth + salt) * sinlat
# ---------------  find site position vector  -----------------
rs = np.empty(3)
rs[0] = rdel * np.cos(lon)
rs[1] = rdel * np.sin(lon)
rs[2] = rk
r = mag(rs)
rp = r + alt
print('rdel %11.7f  rk %11.7f r %11.7f  rpsite %11.7f km \n'
      % (rdel, rk, r, rp))
r = re + alt
fovmin = 0.5 * etafov + etactr
gamma = np.pi - np.arcsin(r * np.sin(fovmin) / re)

print('gamma %11.7f gamma %11.7f rho %11.7f fovmin %11.7f \n'
      % (gamma * rad2deg, (np.pi - gamma) * rad2deg, r, fovmin))
rho = re * np.cos(gamma) + r * np.cos(fovmin)
print('fovmin %11.7f gamma %11.7f gamma %11.7f rho %11.7f  \n'
      % (fovmin * rad2deg, gamma * rad2deg, (np.pi - gamma) * rad2deg, rho))
lambda_ = np.arcsin(rho * np.sin(fovmin) / re)
print('lambda %11.7f  %11.7f km  \n' % (lambda_ * rad2deg, lambda_ * re))
revpday = 16.4
n = revpday * 2 * np.pi / 86400.0

a = (mu * (1 / n) ** 2) ** 0.33333
ecc = (a - rp) / a
print('n %11.7f r/d %11.7f rad/s  a %11.7f km ecc %11.7f  \n'
      % (revpday, n, a, ecc))
print('\n\n now start the design process letting ecc = 0.001 \n\n')
ecc = 0.001
a = rp / (1.0 - ecc)
n = np.sqrt(mu / (a * a * a))
omega = n / earthrot
print('n %11.7f rad/s  a %11.7f km ecc %11.7f omega  %11.7f \n'
      % (n, a, ecc, omega))
a = 6535.4713
ecc = 0.001
pii = 16.387
qi = 1
p = a * (1.0 - ecc * ecc)
n = np.sqrt(mu / (a * a * a))
raandot = - 1.5 * j2 * (re / p) ** 2 * n * np.cos(incl)
nn = pii / qi * (earthrot - raandot) * (1 - 1.5 * j2 * (re / a) ** 2 *
                                        (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
kepperiod = (2.0 * np.pi) / n
nodalperiodG = (2.0 * np.pi) / (earthrot - raandot)
nodalperiod = kepperiod / (1 - 1.5 * j2 * (re / a) ** 2
                           * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
print('p %3.3f q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n'
      % (pii, qi, nodalperiodG / 60, nodalperiod / 60, kepperiod / 60, nn))
print('p %3.3f q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/TU  \n'
      % (pii, qi, nodalperiodG / tusec, nodalperiod / tusec, kepperiod / tusec,
         nn * tusec))

pii = 16.0
qi = 1
p = a * (1.0 - ecc * ecc)
n = np.sqrt(mu / (a * a * a))
raandot = - 1.5 * j2 * (re / p) ** 2 * n * np.cos(incl)
nn = pii / qi * (earthrot - raandot) * (1 - 1.5 * j2 * (re / a) ** 2 *
                                        (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
kepperiod = (2.0 * np.pi) / n
nodalperiodG = (2.0 * np.pi) / (earthrot - raandot)
nodalperiod = kepperiod / (1 - 1.5 * j2 * (re / a) ** 2
                           * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
# nodalperiod = kepperiod * (1 + 1.5 * j2 * (re / a) ** 2
#                            * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n'
      % (pii, qi, nodalperiodG / 60, nodalperiod / 60, kepperiod / 60, nn))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/TU  \n'
      % (pii, qi, nodalperiodG / tusec, nodalperiod / tusec, kepperiod / tusec,
         nn * tusec))
a1 = (mu * (1 / nn) ** 2) ** 0.33333
ecc1 = (a1  - rp) / a1
# Needs second iteration, get a and ecc for final answer

p = a1 * (1.0 - ecc1 * ecc1)
n = np.sqrt(mu / (a1 * a1 * a1))
raandot = - 1.5 * j2 * (re / p) ** 2 * n * np.cos(incl)
nn = pii / qi * (earthrot - raandot) * (1 - 1.5 * j2 * (re / a1) ** 2 *
                                        (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
kepperiod = (2.0 * np.pi) / n
nodalperiodG = (2.0 * np.pi) / (earthrot - raandot)
nodalperiod = kepperiod / (1 - 1.5 * j2 * (re / a1) ** 2
                           * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
# nodalperiod = kepperiod * (1 + 1.5 * j2 * (re / a) ** 2
#                            * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
print('iteration 1 p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n'
      % (pii, qi, nodalperiodG / 60, nodalperiod / 60, kepperiod / 60, nn))
print('iteration 1 p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/TU  \n'
      % (pii, qi, nodalperiodG / tusec, nodalperiod / tusec, kepperiod / tusec,
         nn * tusec))
a2 = (mu * (1 / nn) ** 2) ** 0.33333
ecc2 = (a2  - rp) / a2

pii = 33
qi = 2
p = a * (1.0 - ecc * ecc)
n = np.sqrt(mu / (a * a * a))
raandot = - 1.5 * j2 * (re / p) ** 2 * n * np.cos(incl)
nn = pii / qi * (earthrot - raandot) * (1 - 1.5 * j2 * (re / a) ** 2 *
                                        (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
kepperiod = (2.0 * np.pi) / n
nodalperiodG = (2.0 * np.pi) / (earthrot - raandot)
nodalperiod = kepperiod / (1 - 1.5 * j2 * (re / a) ** 2
                           * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n'
      % (pii, qi, nodalperiodG / 60, nodalperiod / 60, kepperiod / 60, nn))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/TU  \n'
      % (pii, qi, nodalperiodG / tusec, nodalperiod / tusec, kepperiod / tusec,
         nn * tusec))

pii = 49
qi = 3
p = a * (1.0 - ecc * ecc)
n = np.sqrt(mu / (a * a * a))
raandot = - 1.5 * j2 * (re / p) ** 2 * n * np.cos(incl)
nn = pii / qi * (earthrot - raandot) * (1 - 1.5 * j2 * (re / a) ** 2 *
                                        (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
kepperiod = (2.0 * np.pi) / n
nodalperiodG = (2.0 * np.pi) / (earthrot - raandot)
nodalperiod = kepperiod / (1 - 1.5 * j2 * (re / a) ** 2
                           * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n'
      % (pii, qi, nodalperiodG / 60, nodalperiod / 60, kepperiod / 60, nn))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/TU  \n'
      % (pii, qi, nodalperiodG / tusec, nodalperiod / tusec, kepperiod / tusec,
         nn * tusec))

pii = 82
qi = 5
p = a * (1.0 - ecc * ecc)
n = np.sqrt(mu / (a * a * a))
raandot = - 1.5 * j2 * (re / p) ** 2 * n * np.cos(incl)
nn = pii / qi * (earthrot - raandot) * (1 - 1.5 * j2 * (re / a) ** 2 *
                                        (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
kepperiod = (2.0 * np.pi) / n
nodalperiodG = (2.0 * np.pi) / (earthrot - raandot)
nodalperiod = kepperiod / (1 - 1.5 * j2 * (re / a) ** 2
                           * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n'
      % (pii, qi, nodalperiodG / 60, nodalperiod / 60, kepperiod / 60, nn))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/TU  \n'
      % (pii, qi, nodalperiodG / tusec, nodalperiod / tusec, kepperiod / tusec,
         nn * tusec))

pii = 213
qi = 13
p = a * (1.0 - ecc * ecc)
n = np.sqrt(mu / (a * a * a))
raandot = - 1.5 * j2 * (re / p) ** 2 * n * np.cos(incl)
nn = pii / qi * (earthrot - raandot) * (1 - 1.5 * j2 * (re / a) ** 2
                                        * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
kepperiod = (2.0 * np.pi) / n
nodalperiodG = (2.0 * np.pi) / (earthrot - raandot)
nodalperiod = kepperiod / (1 - 1.5 * j2 * (re / a) ** 2
                           * (3.0 - 4.0 * np.sin(incl) * np.sin(incl)))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/s  \n'
      % (pii, qi, nodalperiodG / 60, nodalperiod / 60, kepperiod / 60, nn))
print('p %3i q %3i  pg %11.7f  %11.7f  %11.7f min nn %11.7f rad/TU  \n'
      % (pii, qi, nodalperiodG / tusec, nodalperiod / tusec, kepperiod / tusec,
         nn * tusec))

az = math.asin(math.cos(incl) /  math.cos(latgd))
lambdau = math.acos(math.cos(az) / math.sin(incl))
print(f'az = {az * rad2deg}, lambdau = {lambdau * rad2deg}')
jd = 2451637.0
_, thetagmst = lstime(lon, jd)
lonasc = thetagmst - lambdau + lon
argper = math.asin(math.sin(latgd) / math.sin(incl))
print(f'a = {a2}, ecc = {ecc2}, incl = {incl * rad2deg},'
      f'lonacs = {lonasc * rad2deg} argper = {argper * rad2deg}, M = 0.0000')

