import numpy as np
from space_constants import *
import space_conversions as sc
import spacetime_utils as stu
import orbit_utils as obu
import spacemath_utils as smu


def riset(rsat_eci_initial, vsat_eci_initial, latgd, lon, hellp, jd, dut1, dat, lod, xp, yp, ddpsi, ddeps, rholim, azlim, ellim):
    rsite_ecef = obu.site(latgd,lon,hellp)


    year, mon, day, hr, minute, sec = stu.invjday (jd)
    timezone = 0
    eqeterms  = 2
    ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, \
    jdtdb, jdtdbfrac = stu.convtime(year, mon, day, hr, minute,
                                    sec, timezone, dut1, dat)

    sinlat = math.sin(latgd)
    coslat = math.cos(latgd)
    sinlon = math.sin(lon)
    coslon = math.cos(lon)
    sez2ecef = np.array([[sinlat*coslon, -sinlon, coslat*coslon],
                         [sinlat*sinlon, coslon, coslat*sinlon],
                         [-coslat, 0, sinlat]])
    rsite_sez = sez2ecef @ rsite_ecef

    # Size: 1440 min in a day / 3 min (180sec) intervals
    size = 1440/3
    rsat_ecef_array = np.zeros(size)
    asat_eci = np.zeros(3)
    rsat_ecef_initial, _ = sc.eci2ecef(rsat_eci_initial, vsat_eci_initial, asat_eci, ttt, jdut1 + jdut1frac, lod, xp, yp, eqeterms, ddpsi, ddeps)
    rsat_ecef_array[0] = rsat_ecef_initial
    for i in range(size - 1):
        dtsec = (i+1) * 180
        jdut1new = jdut1 + jdut1frac + dtsec/86400
        rsat_eci_new, vrat_eci_new = obu.pkepler(rsat_eci_initial,vsat_eci_initial,dtsec, 0.0, 0.0)
        rsat_ecef_new, _ = sc.eci2ecef(rsat_eci_new, vrat_eci_new, asat_eci, ttt, jdut1new, lod, xp, yp, eqeterms, ddpsi, ddeps)
        rsat_ecef_array[i+1] = rsat_ecef_new


if __name__ == '__main__':

    # Vernal Equinox of year 2000
    timezone = 0
    year = 2000
    mon = 3
    day = 20
    hr = 7 #UTC
    min = 35
    sec = 0.0
    utc = sec

    # source: https://hpiers.obspm.fr/iers/series/opa/eopc04
    xp = 0.073758
    yp = 0.353317
    dut1 = 0.2847777
    ut1 = utc + dut1
    lod = 0.0010450
    ddpsi = -0.048346
    ddeps = -0.005581
    dat = 32

    jdut1, jdut1frac = stu.jday(year, mon, day, hr, min, ut1)

    # Site: U.S. Air Force Academy
    latgd = 39.007 * deg2rad
    lon = -104.883 * deg2rad
    hellp = 2.918
    rholim = 50000 #km
    azlim = (120,240)
    ellim = (-35, 35)

    # Sat Object 6 data
    n = 16.09769232 * revday2radsec
    ecc = 0.0078742
    incl = 82.8709
    omega = 0.0
    argp = 0.0
    m = 0.0
    nu = m
    a = (mu / n ** 2) ** (1/3)
    p = a * (1 - ecc**2)

    rsat_eci, vsat_eci = sc.coe2rv(p, ecc, incl, omega, argp, nu, 0.0, 0.0, 0.0)

    riset(rsat_eci, vsat_eci, latgd, lon, hellp, jdut1 + jdut1frac, dut1, dat, lod, xp, yp, ddpsi, ddeps, rholim, azlim, ellim)
