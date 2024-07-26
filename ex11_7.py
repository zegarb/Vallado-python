import numpy as np
from space_constants import *
import space_conversions as sc
import spacetime_utils as stu
import orbit_utils as obu
import spacemath_utils as smu

# Make sure time is in UTC
def riset(rsat_eci_initial, vsat_eci_initial, latgd, lon, hellp, jd, dut1, dat, lod, xp, yp, ddpsi, ddeps, rholim, azlim, ellim):

    rsite_ecef, _ = obu.site(latgd,lon,hellp)

    year, mon, day, hr, minute, sec = stu.invjday(jd)
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
    jdut1_initial = jdut1 + jdut1frac
    size = int(1440/3)
    rsat_ecef_array = np.zeros((size,3))
    rhosat_sez_array = np.zeros((size,3))
    frange = np.zeros(size)
    faz_high = np.zeros(size)
    faz_low = np.zeros(size)
    fel_high = np.zeros(size)
    fel_low = np.zeros(size)
    asat_eci = np.zeros(3)
    for i in range(size):
        dtsec = (i) * 180
        jdut1new = jdut1_initial + dtsec/86400
        rsat_eci_new, vrat_eci_new = obu.pkepler(rsat_eci_initial,vsat_eci_initial,dtsec, 0.0, 0.0)
        rsat_ecef_new, _ ,_ = sc.eci2ecef(rsat_eci_new, vrat_eci_new, asat_eci, ttt, jdut1new, lod, xp, yp, eqeterms, ddpsi, ddeps)
        rsat_ecef_array[i] = rsat_ecef_new
        rhosat_ecef = rsat_ecef_new - rsite_ecef
        rhosat_sez_array[i] = sez2ecef @ rhosat_ecef

        rhosatmag = smu.mag(rhosat_sez_array[i])
        frange[i] = rhosatmag - rholim
        faz_high[i] = rhosat_sez_array[i][1] + rhosat_sez_array[i][0] * math.tan(max(azlim))
        faz_low[i] = rhosat_sez_array[i][1] + rhosat_sez_array[i][0] * math.tan(min(azlim))

        rsatmag = smu.mag(rsat_ecef_new)
        el = np.arcsin(rhosat_sez_array[i][2] / rhosatmag)
        fel_high[i] = np.arccos(np.cos(max(ellim))/rsatmag) - max(ellim) - np.arccos(np.cos(el)/rsatmag) + el
        fel_low[i] = np.arccos(np.cos(min(ellim))/rsatmag) - min(ellim) - np.arccos(np.cos(el)/rsatmag) + el

    for i in range(0,480,6):
        print(i)
        minfound1, rootf1, _ = smu.quartbln(frange[i],frange[i+1],frange[i+2],frange[i+3],frange[i+4],frange[i+5])
        if minfound1 == 'y':
            print('root 1: ', rootf1)
        minfound2, rootf2, _ = smu.quartbln(faz_high[i],faz_high[i+1],faz_high[i+2],faz_high[i+3],faz_high[i+4],faz_high[i+5])
        if minfound2 == 'y':
            print('root 2: ', rootf2)
        minfound3, rootf3, _ = smu.quartbln(faz_low[i],faz_low[i+1],faz_low[i+2],faz_low[i+3],faz_low[i+4],faz_low[i+5])
        if minfound3 == 'y':
            print('root 3: ', rootf3)
        minfound4, rootf4, _ = smu.quartbln(fel_high[i],fel_high[i+1],fel_high[i+2],fel_high[i+3],fel_high[i+4],fel_high[i+5])
        if minfound4 == 'y':
            print('root 4: ', rootf4)
        minfound5, rootf5, _ = smu.quartbln(fel_low[i],fel_low[i+1],fel_low[i+2],fel_low[i+3],fel_low[i+4],fel_low[i+5])
        if minfound5 == 'y':
            print('root 5: ', rootf5)

        if minfound1 == 'y' or minfound2 == 'y' or minfound3 == 'y' or minfound4 == 'y' or minfound5 == 'y':
            range_crossing = smu.recovqt(frange[i],frange[i+1],frange[i+2],frange[i+3],frange[i+4],frange[i+5], rootf1)
            az_crossing_low = smu.recovqt(faz_high[i],faz_high[i+1],faz_high[i+2],faz_high[i+3],faz_high[i+4],faz_high[i+5], rootf2)
            az_crossing_high = smu.recovqt(faz_low[i],faz_low[i+1],faz_low[i+2],faz_low[i+3],faz_low[i+4],faz_low[i+5], rootf3)
            el_crossing_low = smu.recovqt(fel_high[i],fel_high[i+1],fel_high[i+2],fel_high[i+3],fel_high[i+4],fel_high[i+5], rootf4)
            el_crossing_high = smu.recovqt(fel_low[i],fel_low[i+1],fel_low[i+2],fel_low[i+3],fel_low[i+4],fel_low[i+5], rootf5)
            print('range_crossing', range_crossing)
            print('az_crossing_low',  az_crossing_low)
            print('az_crossing_high', az_crossing_high)
            print(' el_crosing_low',  el_crossing_low)
            print(' el_crosing_high',  el_crossing_high)

    return frange, faz_high, faz_low, fel_high, fel_low


if __name__ == '__main__':

    # Vernal Equinox of year 2000
    timezone = 0
    year = 2000
    mon = 3
    day = 20
    hr = 7 #UTC
    minute = 35
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

    jdut1, jdut1frac = stu.jday(year, mon, day, hr, minute, ut1)

    # Site: U.S. Air Force Academy
    latgd = 39.007 * deg2rad
    lon = -104.883 * deg2rad
    hellp = 2.918
    rholim = 12000 #km
    azlim = (0*deg2rad,359*deg2rad)
    ellim = (-90*deg2rad, 90*deg2rad)

    # Sat Object 6 data
    n = 16.09769232 * revday2radsec
    ecc = 0.0078742
    incl = 82.8709 * rad2deg
    omega = 0.0
    argp = 0.0
    m = 0.0
    nu = m
    a = (mu / n ** 2) ** (1/3)
    p = a * (1 - ecc**2)

    rsat_eci, vsat_eci = sc.coe2rv(p, ecc, incl, omega, argp, nu, 0.0, 0.0, 0.0)

    frange, faz_high, faz_low, fel_high, fel_low = riset(rsat_eci, vsat_eci, latgd, lon, hellp, jdut1 + jdut1frac, dut1, dat, lod, xp, yp, ddpsi, ddeps, rholim, azlim, ellim)

