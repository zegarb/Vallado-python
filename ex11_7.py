import numpy as np
from space_constants import *
import space_conversions as sc
import spacetime_utils as stu
import orbit_utils as obu
import spacemath_utils as smu


def azse_signcheck(azlim_low, azlim_high, sez_array, highlow):
    az_southsign = -0.0
    az_eastsign = -0.0
    azlim = -0.0
    valid_root = False
    sez_sign = np.sign(sez_array)


    if azlim_high < 0.0:
        azlim_high + 2 * math.pi
    if azlim_low < 0.0:
        azlim_low + 2 * math.pi

    if highlow == 'high':
        azlim = azlim_high
    elif highlow == 'low':
        azlim = azlim_low
    else:
        print('INVALID OPTION')
        return None

    # North Axis
    if np.abs(0.0 - azlim) < small or np.abs(2*math.pi - azlim) < small:
        az_southsign = -1
        az_eastsign = 0
    # East Axis
    elif np.abs(math.pi/2 - azlim) < small:
        az_southsign = 0
        az_eastsign = 1
    # South Axis
    elif np.abs(math.pi - azlim) < small:
        az_southsign = 1
        az_eastsign = 0
    # West Axis
    elif np.abs(3*math.pi/2 - azlim) < small:
        az_southsign = 0
        az_eastsign = -1
    # Quadrant 1
    elif azlim > 0.0 and azlim < math.pi/2:
        az_southsign = -1
        az_eastsign = 1
    # Quadrant 2
    elif azlim > math.pi/2 and azlim < math.pi:
        az_southsign = 1
        az_eastsign = 1
    # Quadrant 3
    elif azlim > math.pi and azlim < 3*math.pi/2:
        az_southsign = 1
        az_eastsign = -1
    # Quadrant 4
    elif azlim >  3*math.pi/2 and azlim < 2*math.pi:
        az_southsign = -1
        az_eastsign = -1

    # If angle falls on south axis only use east component to check angle validity
    if az_southsign == 0 and (sez_sign[1] + az_eastsign) != 0:
        valid_root = True
    # If angle falls on east axis only use south component to check angle validity
    elif az_eastsign == 0 and (sez_sign[0] + az_southsign) != 0:
        valid_root = True
    elif (sez_sign[0] + az_southsign) != 0 and (sez_sign[1] + az_eastsign) != 0:
        valid_root = True

    return valid_root


# Make sure time is in UTC
def riset(rsat_eci_initial, vsat_eci_initial, latgd, lon, hellp, jd, dut1, dat, lod, xp, yp, ddpsi, ddeps, rholim, azlim_low, azlim_high, ellim_low, ellim_high):

    riset_time = []
    riset_az = []
    riset_el = []

    timezone = 0
    eqeterms  = 2

    year, mon, day, hr, minute, sec = stu.invjday(jd)

    rsite_ecef, vsite_ecef = obu.site(latgd,lon,hellp)

    rhosite_sez, _ = sc.ecef2sez(rsite_ecef, vsite_ecef, latgd, lon)

    # Size: 1440 min in a day / 3 min (180sec) intervals
    size = int(1440/1) + 5
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
        rsat_eci_new, vsat_eci_new = obu.pkepler(rsat_eci_initial,vsat_eci_initial,dtsec, 0.0, 0.0)
        # rsat_eci_new, vsat_eci_new, _ = obu.kepler(rsat_eci_initial,vsat_eci_initial,dtsec)

        rsat_ecef_new, vsat_ecef_new ,_ = sc.eci2ecef(rsat_eci_new, vsat_eci_new, asat_eci, ttt, jdut1new, lod, xp, yp, eqeterms, ddpsi, ddeps)
        rsat_ecef_array[i] = rsat_ecef_new

        rhosat_ecef = rsat_ecef_new - rsite_ecef
        drhosat_ecef = vsat_ecef_new

        rhosat_sez_array[i], drhosat_sez_array[i] = sc.ecef2sez(rhosat_ecef, drhosat_ecef, latgd, lon)

        rho, az, el, _, _, _ = sc.sez2razel(rhosat_sez_array[i], drhosat_sez_array[i])

        frange[i] = rho - rholim
        faz_high[i] = rhosat_sez_array[i][1] + rhosat_sez_array[i][0] * math.tan(azlim_high)
        faz_low[i] = rhosat_sez_array[i][1] + rhosat_sez_array[i][0] * math.tan(azlim_low)

        rsatmag = smu.mag(rsat_ecef_new)
        fel_high[i] = (np.arccos(np.cos(el)/rsatmag) - el) - (np.arccos(np.cos(ellim_high)/rsatmag) - ellim_high)
        fel_low[i] = (np.arccos(np.cos(el)/rsatmag) - el) - (np.arccos(np.cos(ellim_low)/rsatmag) - ellim_low)
        rho, az, el, drho, daz, del_ = sc.sez2razel(rhosat_sez_array[i], drhosat_sez_array[i])

    rho, az, el, drho, daz, del_ = sc.sez2razel(rhosat_sez_array[0], drhosat_sez_array[0])
    if rho <= rholim:
        range_inview = True
    else:
        range_inview = False

    if az >= azlim_high and az <= azlim_low:
        az_inview = True
    else:
        az_inview = False

    if el >= ellim_low and el <= ellim_high:
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
    full_az = False
    full_el = False

    #full 360 degree view check
    if (0.0 - azlim_low) < small and (2*math.pi - azlim_high) < small:
        full_az = True

    if (-math.pi - ellim_low) < small and (math.pi - ellim_high) < small:
        full_el = True

    for i in range(0,480):
        # print(i)
        az_rootvalid = True
        minfound = False
        minfound1 = False
        minfound2 = False
        minfound3 = False
        minfound4 =  False
        minfound5 = False

        minfound1, rootf1, funrate1 = smu.quartbln(frange[i],frange[i+1],frange[i+2],frange[i+3],frange[i+4],frange[i+5])
        if not(full_az):
            minfound2, rootf2, funrate2 = smu.quartbln(faz_high[i],faz_high[i+1],faz_high[i+2],faz_high[i+3],faz_high[i+4],faz_high[i+5])
            minfound3, rootf3, funrate3 = smu.quartbln(faz_low[i],faz_low[i+1],faz_low[i+2],faz_low[i+3],faz_low[i+4],faz_low[i+5])
        if not(full_el):
            minfound4, rootf4, funrate4 = smu.quartbln(fel_high[i],fel_high[i+1],fel_high[i+2],fel_high[i+3],fel_high[i+4],fel_high[i+5])
            minfound5, rootf5, funrate5 = smu.quartbln(fel_low[i],fel_low[i+1],fel_low[i+2],fel_low[i+3],fel_low[i+4],fel_low[i+5])
            # print(fel_low[i],fel_low[i+1],fel_low[i+2],fel_low[i+3],fel_low[i+4],fel_low[i+5])

        if minfound1:
            minfound = True
            combined_root = rootf1
            print('root 1: ', rootf1)
            # print('fundrate 1: ', funrate1)
        if minfound2:
            minfound = True
            combined_root = rootf2
            print('root 2: ', rootf2)
            print('fundrate 2: ', funrate2)
        if minfound3:
            minfound = True
            combined_root = rootf3
            print('root 3: ', rootf3)
            print('fundrate 3: ', funrate3)
        if minfound4:
            minfound = True
            combined_root = rootf4
            print('root 4: ', rootf4)
            # print('fundrate 4: ', funrate4)
        if minfound5:
            minfound = True
            combined_root = rootf5
            print('root 5: ', rootf5)
            # print('fundrate 5: ', funrate5)

        if minfound:
            print(i)
            view_time, _ = smu.recovqt(i*180, (i+1)*180, (i+2)*180, (i+3)*180, (i+4)*180, (i+5)*180, combined_root)
            rhosat_s, drhosat_s = smu.recovqt(rhosat_sez_array[i][0], rhosat_sez_array[i+1][0], rhosat_sez_array[i+2][0], rhosat_sez_array[i+3][0], rhosat_sez_array[i+4][0], rhosat_sez_array[i+5][0], combined_root)
            rhosat_e, drhosat_e = smu.recovqt(rhosat_sez_array[i][1], rhosat_sez_array[i+1][1], rhosat_sez_array[i+2][1], rhosat_sez_array[i+3][1], rhosat_sez_array[i+4][1], rhosat_sez_array[i+5][1], combined_root)
            rhosat_z, drhosat_z = smu.recovqt(rhosat_sez_array[i][2], rhosat_sez_array[i+1][2], rhosat_sez_array[i+2][2], rhosat_sez_array[i+3][2], rhosat_sez_array[i+4][2], rhosat_sez_array[i+5][2], combined_root)
            rhosat_sez = np.array([rhosat_s,rhosat_e,rhosat_z])
            drhosat_sez = np.array([drhosat_s,drhosat_e,drhosat_z])

            # Angles are seperated at 180 degrees.
            if minfound2 and minfound3:
                az_rootvalid = True
                if azse_signcheck(azlim_low, azlim_high, rhosat_sez, 'high'):
                    print('high azimuth')
                    minfound3 = False
                elif azse_signcheck(azlim_low, azlim_high, rhosat_sez, 'low'):
                    print('low azimuth')
                    minfound2 =  False

            elif minfound2:
                az_rootvalid = azse_signcheck(azlim_low, azlim_high, rhosat_sez, 'high')
            elif minfound3:
                az_rootvalid = azse_signcheck(azlim_low, azlim_high, rhosat_sez, 'low')

            if az_rootvalid:
                rho,az,el, drho, daz,del_ = sc.sez2razel(rhosat_sez,drhosat_sez)

                if az < 0.0:
                    az = 2 * math.pi + az

                if minfound1 and drho < 0.0:
                    range_inview = True
                elif not(minfound1) and rho < rholim:
                    range_inview = True
                else:
                    range_inview = False

                if minfound2 and daz < 0:
                    print('here 1')
                    az_inview = True
                elif minfound3 and daz > 0:
                    print('here 2')
                    az_inview = True
                elif not(minfound2 or minfound3) and (az >= azlim_low and az <= azlim_high):
                    print('here 3')
                    az_inview = True
                elif full_az:
                    print('here 4')
                    az_inview = True
                else:
                    az_inview = False

                if minfound4 and del_ < 0.0:
                    el_inview= True
                elif minfound5 and del_ > 0.0:
                    el_inview = True
                elif not(minfound4 or minfound5) and (el >= ellim_low and el <= ellim_high):
                    el_inview = True
                elif full_el:
                    el_inview = False
                else:
                    el_inview = False


                newview = range_inview and az_inview and el_inview

                print('az view: ', az_inview)
                print('el view: ', el_inview)
                print('range view: ', range_inview)
                print('az: ', az)
                print('daz: ', daz)
                print('el: ', el)
                # print('range: ', rho)

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
                oldview = newview

            else:
                rho,az,el, drho, daz,del_ = sc.sez2razel(rhosat_sez,drhosat_sez)
                if az < 0.0:
                    az = 2 * math.pi + az
                print('az: ', az)
                print('Invalid Azimuth Root\n')

            minfound = False

    return riset_time, riset_az, riset_el

if __name__ == '__main__':

    # # Vernal Equinox of year 2000
    # timezone = 0
    # year = 2000
    # mon = 3
    # day = 20
    # hr = 7 #UTC
    # minute = 35
    # sec = 0.0
    # utc = sec

    # # source: https://hpiers.obspm.fr/iers/series/opa/eopc04
    # xp = 0.073758
    # yp = 0.353317
    # dut1 = 0.2847777
    # ut1 = utc + dut1
    # lod = 0.0010450
    # ddpsi = -0.048346
    # ddeps = -0.005581
    # dat = 32

    # Vernal Equinox of year 1992
    timezone = 0
    year = 1992
    mon = 3
    day = 20
    hr = 8 #UTC
    minute = 47
    sec = 0.0
    utc = sec

    # source: https://hpiers.obspm.fr/iers/series/opa/eopc04
    xp = -0.054482
    yp = 0.137482
    dut1 = -0.3234685
    ut1 = utc + dut1
    lod = 0.033059
    ddpsi = -0.009127
    ddeps = -0.004187
    dat = 27

    jdut1, jdut1frac = stu.jday(year, mon, day, hr, minute, ut1)

    # Site: U.S. Air Force Academy
    latgd = 39.007 * deg2rad
    lon = -104.883 * deg2rad
    hellp = 2.918

    # Limits depend on antenna site survey. Values below are arbitrary.
    rholim = 100000 #km
    # Superimposed circle over actual azimuth/elevation limits
    azlim_low = 0.0*deg2rad
    azlim_high = 360*deg2rad
    ellim_low = 30*deg2rad
    ellim_high = 90*deg2rad
    objectnum = 8

    # Sat Objects 1, 2 (?), 4, and 7 are not visible to site
    match objectnum:
        case 1:
            n = 1.00272141 * revday2radsec
            ecc = 0.0000032
            incl = 0.0956 * deg2rad
        case 2:
            n = 8.36589235 * revday2radsec
            ecc = 0.0080158
            incl = 90.0175 * deg2rad
        case 3:
            n = 0.24891961 * revday2radsec
            ecc = 0.9363060
            incl = 64.9874 * deg2rad
        case 4:
            n = 0.21467209 * revday2radsec
            ecc = 0.0668128
            incl = 57.3500 * deg2rad
        case 5:
            n = 13.37659679 * revday2radsec
            ecc = 0.0145072
            incl = 90.2619 * deg2rad
        case 6:
            n = 16.09769232 * revday2radsec
            ecc = 0.0078742
            incl = 82.8709 * deg2rad
        case 7:
            n = 1.00271920 * revday2radsec
            ecc = 0.0003109
            incl = 0.0099 * deg2rad
        case 8:
            n = 12.41552416 * revday2radsec
            ecc = 0.0036498
            incl = 74.0186 * deg2rad
        case 9:
            n = 13.84150848 * revday2radsec
            ecc = 0.0048964
            incl = 144.6414 * deg2rad

    omega = 0.0
    argp = 0.0
    m = 0.0
    nu = m
    a = (mu / n ** 2) ** (1/3)
    p = a * (1 - ecc**2)

    rsat_eci, vsat_eci = sc.coe2rv(p, ecc, incl, omega, argp, nu, 0.0, 0.0, 0.0)

    riset_time, riset_az, riset_el = riset(rsat_eci, vsat_eci, latgd, lon, hellp, jdut1 + jdut1frac, dut1, dat, lod, xp, yp, ddpsi, ddeps, rholim, azlim_low, azlim_high, ellim_low, ellim_high)

    if len(riset_time) == 0:
        print('Satellite is not visible to site')
    else:
        for i in range(len(riset_time)):
            print(f'rise time: {riset_time[i][0]}')
            print(f'rise az  : {riset_az[i][0]}')
            print(f'rise el  : {riset_el[i][0]} \n')
            print(f'set time: {riset_time[i][1]}')
            print(f'set az  : {riset_az[i][1]}')
            print(f'set el  : {riset_el[i][1]}\n\n')