#     -----------------------------------------------------------------
#
#                              Ex11_1
#
#  this file demonstrates example 11-1.
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
#             7 jun 07  david vallado
#                         original
#  changes :
#            13 feb 07  david vallado
#                         original baseline
#
#     *****************************************************************

import numpy as np
import spacetime_utils as stu
import orbit_utils as obu
from space_constants import *
import space_conversions as sc

print('---------------Example 11-1 --------------------------\n')

# --------  satfov calculations
incl = 60.0 * deg2rad
az = 40.0 * deg2rad
slatgd = 50.0 * deg2rad
slon = 40.0 * deg2rad
salt = 800.0

tfov = 25.0 * deg2rad
etactr = 0.0 * deg2rad
rhomin, rhomax, lambda_ = obu.satfov(incl, az, slatgd, slon, salt, tfov, etactr)
print(f'{rhomax = } km, {rhomin = } km, {lambda_ = } rad')
print(f'lambda = {lambda_ * rad2deg}')
print('\noff axis test')
# run twice with +- tfov angles to get the max and min from nadir
# location
incl = 60.0 * deg2rad
az = 40.0 * deg2rad
slatgd = 50.0 * deg2rad
slon = 40.0 * deg2rad
salt = 800.0

tfov = 25.0 * deg2rad
etactr = 40.0 * deg2rad
rhomin, rhomax, lambda_ = obu.satfov(incl, az, slatgd, slon, salt, tfov, etactr)
print(f'{rhomax = } km, {rhomin = } km, {lambda_ = } rad')
print(f'lambda = {lambda_ * rad2deg}')
print('\noff axis test 2nd half \n')
incl = 60.0 * deg2rad
az = 40.0 * deg2rad
slatgd = 50.0 * deg2rad
slon = 40.0 * deg2rad
salt = 800.0

tfov = -25.0 * deg2rad
etactr = 40.0 * deg2rad
rhomin, rhomax, lambda_ = obu.satfov(incl, az, slatgd, slon, salt, tfov, etactr)
print('\noff axis test circle calcs \n')
incl = 60.0 * deg2rad
az = 140.0 * deg2rad
slatgd = 50.0 * deg2rad
slon = 40.0 * deg2rad
salt = 800.0

tfov = 25.0 * deg2rad
etactr = 27.5 * deg2rad
# location sat fov sensor if on satellite direction of motion az
az = np.arcsin(np.cos(incl) / np.cos(slatgd))
u = np.arcsin(np.sin(slatgd) / np.sin(incl))
if (np.cos(u) < 0.0):
    az = np.pi - az

print('az %11.7f %11.7f u %11.7f  %11.7f \n\n' % (az * rad2deg,
                                                  (np.pi - az) * rad2deg,
                                                  u * rad2deg,
                                                  (np.pi - u) * rad2deg))
# pick az = 140 deg
az = 140.0 * deg2rad
rhomin, rhomax, lambda_ = obu.satfov(incl, az, slatgd, slon, salt, tfov, etactr)
# ----- loop around the new circle with the sensor range ------
for i in range(0, 19):
    az = i * 20.0 * deg2rad
    tgtlat, tgtlon = obu.pathm(slatgd, slon, rhomin * 0.5 / re, az)
    if (i == 0):
        maxlat = tgtlat
    if (i == 9):
        minlat = tgtlat
    print('az %11.7f lat %11.7f lon %11.7f\n' % (az * rad2deg,
                                                 tgtlat * rad2deg,
                                                 tgtlon * rad2deg))
