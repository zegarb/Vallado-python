#     -----------------------------------------------------------------
#
#                              Ex11_4
#
#  this file demonstrates example 11-4.
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
#            30 aug 11  david vallado
#                         original
#  changes :
#            22 aug 11 david vallado
#                         original baseline
#
#     *****************************************************************

import numpy as np
from space_constants import *

print("------ Example 11-4 --------")
gravc = getgravc('egm08')
j2 = gravc['j2']
mu = gravc['mu']
re = gravc['re']
# --------  repeat gt calculations
a = 6570.3358

ecc = 0.006301
incl = 45.0 * deg2rad
p = a * (1.0 - ecc * ecc)
nanom = np.sqrt(mu / (a * a * a))
print('p %11.7f  %11.7f \n' % (p, p / re))
n = nanom
print('n %11.7f  %11.7f \n' % (n, n * deg2rad))
raanrate = - 1.5 * j2 * re ** 2 * n * np.cos(incl) / (p * p)
print('raanrate %11.7f  %11.7f \n' % (raanrate, raanrate * 180 * 86400 / np.pi))

deltatoa = (3.0 * np.pi / (n * a)) * (1.0 + 0.5 * j2 * (re / a) ** 2 * (4.0 * (np.cos(incl)) ** 2 - 1.0))
print('deltatoa %11.7f  %11.7f \n' % (deltatoa, deltatoa * re / tusec))
deltaOdotoa = - 3.5 * raanrate / a
print('deltaOdotoa %d  %11.11f \n' % (deltaOdotoa, deltaOdotoa * re / tusec))
deltapoi = 12.0 * np.pi / n * j2 * (re / a) ** 2 * np.sin(2.0 * incl)
print('deltapoi %11.7f  %11.7f \n' % (deltapoi, deltapoi / tusec))
deltaraanoi = - raanrate * np.tan(incl)
print('deltaraanoi %11.7e  %11.7e \n' % (deltaraanoi, deltaraanoi / tusec))
pnodal = 2.0 * np.pi / n
print('pnodal %11.7f s %11.7f min \n' % (pnodal, pnodal / 60.0))
anomper = pnodal * (1.0 / (1.0 + 0.75 * j2 * (re / p) ** 2 * (np.sqrt(1.0 - ecc * ecc) * (2.0 - 3.0 * np.sin(incl) ** 2) + (4.0 - 5.0 * np.sin(incl) ** 2))))
print('anomper %11.7f s %11.7f min \n' % (anomper, anomper / 60.0))
dellon = (earthrot - raanrate) * anomper
print('dellon %11.7f  %11.7f \n' % (dellon, dellon * re))
dlpa = re * (earthrot - raanrate) * deltatoa - deltaOdotoa * re * anomper
print('dlpa %11.7f  %11.7f \n' % (dlpa, dlpa))
dlpi = re * (earthrot - raanrate) * deltapoi - deltaraanoi * re * anomper
print('dlpi %11.7f  %11.7f \n' % (dlpi, dlpi))
# osculating values from stk
#2 Apr 2000 12:00:00.000  6570.3400000 0.00630100  45.00000  321.21600  69.629    0.000   5300.206
#2 Apr 2000 13:28:20.000  6565.1860071 0.00579447  44.99776  320.82382  70.494  359.550   5293.971
# mean values from stk

dadt = - 0.174562 / anomper
didt = 0.000132 * np.pi / (180 * anomper)

#        dadt = -1.59082 / anomper;
#        didt =  0.01188 * pi / (180 * anomper);

print('dadt %11.7f km/rev %e km/s \n' % (dadt * anomper, dadt))
print('didt %11.7f deg/rev %e deg/sec \n' % (didt * anomper * 180 / np.pi, didt))
k2 = 1.0 / anomper * (dlpa * dadt + dlpi * didt)
k1 = np.sqrt(2.0 * k2 * (- 50 - 50))
print('k2 %e  %11.7e \n' % (k2, k2 * tusec / re))
print('k1 %11.7f  %11.7f \n' % (k1, k1 * tusec / re))
da = k1 * anomper / dlpa
print('da %11.7f  %11.7f \n' % (da, da))
tdrift = - k1 / k2
print('tdrift %11.7f sec  %11.7f min %11.7f day \n' % (tdrift, tdrift / 60.0, tdrift / 86400.0))
deltav = n * 0.5 * da
print('deltav %11.7f km/s %11.7f m/s \n' % (deltav, deltav * 1000))
print('change frequency deltav %11.7f m/s timing  %11.7f min, %11.7f days \n'
      % (2 * deltav * 1000, 2 * tdrift / 60, 2 * tdrift / 86400))
