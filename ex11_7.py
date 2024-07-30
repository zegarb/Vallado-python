import numpy as np
from space_constants import *
import space_conversions as sc
import spacetime_utils as stu
import orbit_utils as obu
import spacemath_utils as smu

# Make sure time is in UTC
def riset(rsat_eci_initial, vsat_eci_initial, latgd, lon, hellp, jd, dut1, dat, lod, xp, yp, ddpsi, ddeps, rholim, azlim, ellim):

    riset_time = []
    riset_az = []
    riset_el = []

    timezone = 0
    eqeterms  = 2

    rsite_ecef, _ = obu.site(latgd,lon,hellp)

    year, mon, day, hr, minute, sec = stu.invjday(jd)

    sinlat = math.sin(latgd)
    coslat = math.cos(latgd)
    sinlon = math.sin(lon)
    coslon = math.cos(lon)
    ecef2sez = np.array([[sinlat*coslon, sinlat*sinlon, -coslat],
                         [-sinlon, coslon, 0],
                         [coslat*coslon, coslat*sinlon, sinlat]])

    rsite_sez = ecef2sez @ rsite_ecef

    # Size: 1440 min in a day / 3 min (180sec) intervals
    size = int(1440/3) + 5
    rsat_ecef_array = np.zeros((size,3))
    rhosat_sez_array = np.zeros((size,3))
    drhosat_sez_array = np.zeros((size,3))
    frange = np.zeros(size)
    faz_high = np.zeros(size)
    faz_low = np.zeros(size)
    fel_high = np.zeros(size)
    fel_low = np.zeros(size)
    asat_eci = np.zeros(3)
    for i in range(size):
        dtsec = (i) * 180
        ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, \
        jdtdb, jdtdbfrac = stu.convtime(year, mon, day, hr, minute,
                                    sec + dtsec, timezone, dut1, dat)
        jdut1new = jdut1 + jdut1frac
        rsat_eci_new, vrat_eci_new = obu.pkepler(rsat_eci_initial,vsat_eci_initial,dtsec, 0.0, 0.0)
        rsat_ecef_new, vsat_ecef_new ,_ = sc.eci2ecef(rsat_eci_new, vrat_eci_new, asat_eci, ttt, jdut1new, lod, xp, yp, eqeterms, ddpsi, ddeps)
        rsat_ecef_array[i] = rsat_ecef_new

        rhosat_ecef = rsat_ecef_new - rsite_ecef
        drhosat_ecef = vsat_ecef_new

        rhosat_sez_array[i] = ecef2sez @ rhosat_ecef
        drhosat_sez_array[i] = ecef2sez @ drhosat_ecef

        rho, az, el, _, _, _ = sc.sez2razel(rhosat_sez_array[i], drhosat_sez_array[i])

        frange[i] = rho - rholim
        faz_high[i] = rhosat_sez_array[i][1] - rhosat_sez_array[i][0] * math.tan(max(azlim))
        faz_low[i] = rhosat_sez_array[i][1] - rhosat_sez_array[i][0] * math.tan(min(azlim))

        rsatmag = smu.mag(rsat_ecef_new)
        fel_high[i] = (np.arccos(np.cos(el)/rsatmag) - el) - (np.arccos(np.cos(max(ellim))/rsatmag) - max(ellim))
        fel_low[i] = (np.arccos(np.cos(el)/rsatmag) - el) - (np.arccos(np.cos(min(ellim))/rsatmag) - min(ellim))

    rho, az, el, drho, daz, del_ = sc.sez2razel(rhosat_sez_array[0], drhosat_sez_array[0])
    if rho <= rholim:
        range_inview = True
    else:
        range_inview = False

    if az >= min(azlim) and az <= max(azlim):
        az_inview = True
    else:
        az_inview = False

    if el >= min(ellim) and el <= max(ellim):
        el_inview = True
    else:
        el_inview = False

    view_time = 0.0
    oldview = range_inview and az_inview and el_inview
    if oldview:
        print('in view at (s): ', view_time)
        riset_time.append([view_time, None])
        riset_az.append([az, None])
        riset_el.append([el, None])

    # print('az view: ', az_inview)
    # print('el view: ', el_inview)
    # print('range view: ', range_inview)
    combined_root = -0.0
    minfound = False
    for i in range(0,480):
        # if i > 273 and i < 297:
        print(i)
            # print('frange: ',frange[i],frange[i+1],frange[i+2],frange[i+3],frange[i+4],frange[i+5])
            # print('faz_high: ',faz_high[i],faz_high[i+1],faz_high[i+2],faz_high[i+3],faz_high[i+4],faz_high[i+5])
            # print('faz_low: ',faz_low[i],faz_low[i+1],faz_low[i+2],faz_low[i+3],faz_low[i+4],faz_low[i+5])
            # print('fel_high: ',fel_high[i],fel_high[i+1],fel_high[i+2],fel_high[i+3],fel_high[i+4],fel_high[i+5])
            # print('fel_low: ',fel_low[i],fel_low[i+1],fel_low[i+2],fel_low[i+3],fel_low[i+4],fel_low[i+5])
        minfound1, rootf1, funrate1  = smu.quartbln(frange[i],frange[i+1],frange[i+2],frange[i+3],frange[i+4],frange[i+5])
        if minfound1 == True:
            minfound = True
            combined_root = rootf1
            print('root 1: ', rootf1)
            # print('fundrate 1: ', funrate1)
        minfound2, rootf2, funrate2 = smu.quartbln(faz_high[i],faz_high[i+1],faz_high[i+2],faz_high[i+3],faz_high[i+4],faz_high[i+5])
        if minfound2 == True:
            print('faz_high: ',faz_high[i],faz_high[i+1],faz_high[i+2],faz_high[i+3],faz_high[i+4],faz_high[i+5])
            minfound = True
            combined_root = rootf2
            print('root 2: ', rootf2)
            # print('fundrate 2: ', funrate2)
        minfound3, rootf3, funrate3 = smu.quartbln(faz_low[i],faz_low[i+1],faz_low[i+2],faz_low[i+3],faz_low[i+4],faz_low[i+5])
        if minfound3 == True:
            print('faz_low: ',faz_low[i],faz_low[i+1],faz_low[i+2],faz_low[i+3],faz_low[i+4],faz_low[i+5])
            minfound = True
            combined_root = rootf3
            print('root 3: ', rootf3)
            # print('fundrate 3: ', funrate3)
        minfound4, rootf4, funrate4 = smu.quartbln(fel_high[i],fel_high[i+1],fel_high[i+2],fel_high[i+3],fel_high[i+4],fel_high[i+5])
        if minfound4 == True:
            minfound = True
            combined_root = rootf4
            print('root 4: ', rootf4)
            # print('fundrate 4: ', funrate4)
        minfound5, rootf5, funrate5 = smu.quartbln(fel_low[i],fel_low[i+1],fel_low[i+2],fel_low[i+3],fel_low[i+4],fel_low[i+5])
        if minfound5 == True:
            minfound = True
            combined_root = rootf5
            print('root 5: ', rootf5)
            # print('fundrate 5: ', funrate5)

        if minfound == True:

            view_time, _ = smu.recovqt(i*180, (i+1)*180, (i+2)*180, (i+3)*180, (i+4)*180, (i+5)*180, combined_root)
            rhosat_s, drhosat_s = smu.recovqt(rhosat_sez_array[i][0], rhosat_sez_array[i+1][0], rhosat_sez_array[i+2][0], rhosat_sez_array[i+3][0], rhosat_sez_array[i+4][0], rhosat_sez_array[i+5][0], combined_root)
            rhosat_e, drhosat_e = smu.recovqt(rhosat_sez_array[i][1], rhosat_sez_array[i+1][1], rhosat_sez_array[i+2][1], rhosat_sez_array[i+3][1], rhosat_sez_array[i+4][1], rhosat_sez_array[i+5][1], combined_root)
            rhosat_z, drhosat_z = smu.recovqt(rhosat_sez_array[i][2], rhosat_sez_array[i+1][2], rhosat_sez_array[i+2][2], rhosat_sez_array[i+3][2], rhosat_sez_array[i+4][2], rhosat_sez_array[i+5][2], combined_root)

            rhosat_sez = np.array([rhosat_s,rhosat_e,rhosat_z])
            drhosat_sez = np.array([drhosat_s,drhosat_e,drhosat_z])
            rho,az,el, drho, daz,del_ = sc.sez2razel(rhosat_sez,drhosat_sez)
            if az < 0.0:
                az = math.pi + az


            if minfound1 and drho < 0.0:
                range_inview = True
            elif not(minfound1) and rho < rholim:
                range_inview = True
            else:
                range_inview = False

            #full 360 degree view check
            if min(azlim) == 0.0 and round(max(azlim),4) == round(2*math.pi,4):
                az_inview = True
            elif minfound2 and daz < 0:
                az_inview = True
            elif minfound3 and daz > 0:
                az_inview = True
            elif not(minfound2 or minfound3) and (az > min(azlim) and az < max(azlim)):
                az_inview = True
            else:
                az_inview = False

            if minfound4 and del_ < 0.0:
                el_inview= True
            elif minfound5 and del_ > 0.0:
                el_inview = True
            elif not(minfound4 or minfound5) and (el > min(ellim) and el < min(ellim)):
                el_inview = True
            else:
                el_inview = False


            newview = range_inview and az_inview and el_inview

            print('az view: ', az_inview)
            print('el view: ', el_inview)
            print('range view: ', range_inview)
            print('az: ', az)
            print('el: ', el)
            print('del;', del_)
            #print('range: ', rho)

            if newview == True and newview != oldview:
                print('in view at (s): ', view_time)
                riset_time.append([view_time, None])
                riset_az.append([az, None])
                riset_el.append([el, None])
            elif newview == False and newview != oldview:
                print('out of view at (s): ', view_time)
                time_temp = riset_time.pop()
                time_temp[1] = view_time
                riset_time.append(time_temp)

                az_temp = riset_az.pop()
                az_temp[1] = az
                riset_az.append(az_temp)

                el_temp = riset_el.pop()
                el_temp[1] = el
                riset_el.append(el_temp)

            minfound = False
            oldview = newview

    return riset_time, riset_az, riset_el

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
    rholim = 100000 #km
    # Superimposed circle over actual
    azlim = (0*deg2rad,360*deg2rad)
    ellim = (30*deg2rad, 90*deg2rad)

    # Sat Object 8 data
    n = 12.41552416 * revday2radsec
    ecc = 0.0036498
    incl = 74.0186 * deg2rad
    omega = 0.0
    argp = 0.0
    m = 0.0
    nu = m
    a = (mu / n ** 2) ** (1/3)
    p = a * (1 - ecc**2)

    rsat_eci, vsat_eci = sc.coe2rv(p, ecc, incl, omega, argp, nu, 0.0, 0.0, 0.0)

    riset_time, riset_az, riset_el = riset(rsat_eci, vsat_eci, latgd, lon, hellp, jdut1 + jdut1frac, dut1, dat, lod, xp, yp, ddpsi, ddeps, rholim, azlim, ellim)


    for i in range(len(riset_time)):
        print(f'rise time: {riset_time[i][0]}')
        print(f'rise az  : {riset_az[i][0]}')
        print(f'rise el  : {riset_el[i][0]} \n')
        print(f'set time: {riset_time[i][1]}')
        print(f'set az  : {riset_az[i][1]}')
        print(f'set el  : {riset_el[i][1]}\n\n')