import math
import numpy as np
from pprint import pprint as pp
from space_constants import *

# ------------------------------------------------------------------------------
#
#                           function jdayall
#
#  this function finds the julian date given the year, month, day, and time.
#    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - universal time hour            0 .. 23
#    min         - universal time min             0 .. 59
#    sec         - universal time sec             0.0 .. 59.999
#    whichtype   - julian or gregorian calender   'j' or 'g'
#
#  outputs       :
#    jd          - julian date                    days from 4713 bc
#
#  locals        :
#    b           - var to aid gregorian dates
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2001, 187
#
# [jd, jdfrac] = jdayall(year, mon, day, hr, min, sec, whichtype)
# -----------------------------------------------------------------------------

def jdayall(year=None, mon=None, day=None, hr=None, min=None, sec=None, whichtype=None):
    if mon <= 2:
        year = year - 1
        mon = mon + 12

    if whichtype == 'j':
        # --------- use for julian calender, every 4 years --------
        b = 0.0
    else:
        # ---------------------- use for gregorian ----------------
        b = (2 - int(np.floor(year * 0.01))
             + int(np.floor(int(np.floor(year * 0.01)) * 0.25)))

    jd = (int(np.floor(365.25 * (year + 4716)))
          + int(np.floor(30.6001 * (mon + 1))) + day + b - 1524.5)
    jdfrac = (sec + min * 60.0 + hr * 3600.0) / 86400.0
    # check jdfrac
    if jdfrac > 1.0:
        jd = jd + int(np.floor(jdfrac))
        jdfrac = jdfrac - int(np.floor(jdfrac))

    return jd, jdfrac


#
# -----------------------------------------------------------------------------
#
#                           function getintda
#
#  this function finds the integer equivalent of the 3 character string
#    representation of the day of the week.
#
#  author        : david vallado                  719-573-2600    5 jul 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    daystr      - day name string                'sun', 'mon' ...
#
#  outputs       :
#    dayn        - integer day equivalent         1 .. 7
#
#  locals        :
#    i           - index
#
#  coupling      :
#    none
#
# [dayn] = getintda(daystr)
# -----------------------------------------------------------------------------

def getintda(daystr=None):
    # ------------------------  implementation   --------------------------
    daytitle = ['sun', 'mon', 'tue', 'wed', 'thr', 'fri', 'sat']
    dayn = daytitle.index(daystr) + 1
    return dayn

# -----------------------------------------------------------------------------
#
#                           function dayofwee
#
#  this function finds the day of the week. integers are used for the days,
#    1 = 'sun', 2 = 'mon', ... 7 = 'sat'.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    jd          - julian date of interest        days from 4713 bc
#
#  outputs       :
#    dayofweek   - answer                         1 to 7
#
#  locals        :
#    none.
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 188, eq 3-39
#
# dayofweek = dayofwee(jd)
# -----------------------------------------------------------------------------

def dayofwee(jd=None):
    # ------------------------  implementation   ------------------
# ------- be sure jd is at 0.0d0 h on the day of interest -----
    jd = int(np.floor(jd + 0.5))
    dayofweek = np.rint(jd - 7 * np.rint((jd + 1) / 7) + 2)
    return dayofweek

#
# -----------------------------------------------------------------------------
#
#                           function getintmon
#
#  this function finds the integer equivalent of the 3 character string
#    representation of month.
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    monstr      - month name                     'jan', 'feb' ...
#
#  outputs       :
#    mon         - integer month equivalent       1 .. 12
#
#  locals        :
#    i           - index
#
#  coupling      :
#    none
#
# [mon] = getintmon(monstr)
# -----------------------------------------------------------------------------


def getintmon(monstr=None):
    # ------------------------  implementation   --------------------------
    monthtitle = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
    mon = monthtitle.index(monstr) + 1
    return mon



# -----------------------------------------------------------------------------
#
#                           function hms2sec
#
#  this function converts hours, minutes and seconds into seconds from the
#    beginning of the day.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    hr          - hours                          0 .. 24
#    min         - minutes                        0 .. 59
#    sec         - seconds                        0.0 .. 59.99
#
#   outputs      :
#    utsec       - seconds                        0.0 .. 86400.0
#
#  locals        :
#    temp        - temporary variable
#
#  coupling      :
#    none.
#
# function [utsec ] = hms2sec(hr, min, sec)
# -----------------------------------------------------------------------------


def hms2sec(hr, min, sec):

        # ------------------------  implementation   ------------------
        utsec = hr * 3600.0 + min * 60.0 + sec
        return utsec


# -----------------------------------------------------------------------------
#
#                           function sec2hms
#
#  this function converts seconds from the beginning of the day into hours,
#    minutes and seconds.
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    utsec       - seconds                        0.0 .. 86400.0
#
#  outputs       :
#    hr          - hours                          0 .. 24
#    min         - minutes                        0 .. 59
#    sec         - seconds                        0.0 .. 59.99
#
#  locals        :
#    temp        - temporary variable
#
#  coupling      :
#    none.
#
# [hr, min, sec] = sec2hms(utsec)
# -----------------------------------------------------------------------------


def sec2hms(utsec):

        # ------------------------  implementation   ------------------
        temp = utsec / 3600.0
        hr = np.fix(temp)
        min = np.fix((temp - hr) * 60.0)
        sec = (temp - hr - min/60.0) * 3600.0
        return hr, min, sec



# -----------------------------------------------------------------------------
#
#                           function jday.m
#
#  this function finds the julian date given the year, month, day, and time.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - universal time hour            0 .. 23
#    min         - universal time min             0 .. 59
#    sec         - universal time sec             0.0 .. 59.999
#    whichtype   - julian .or. gregorian calender   'j' .or. 'g'
#
#  outputs       :
#    jd          - julian date                    days from 4713 bc
#    jdfrac      - julian date fraction of a day   0.0 to 1.0
#
#  locals        :
#    none.
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 189, alg 14, ex 3-14
#
# [jd, jdfrac] = jday(yr, mon, day, hr, min, sec)
# -----------------------------------------------------------------------------


def jday(yr, mon, day, hr, min, sec):

        # ------------------------  implementation   ------------------
        jd = 367.0 * yr  \
             - math.floor((7 * (yr + math.floor((mon + 9) / 12.0))) * 0.25)   \
             + math.floor(275 * mon / 9.0) \
             + day + 1721013.5   # use - 678987.0 to go to mjd directly
        jdfrac = (sec + min * 60.0 + hr *3600.0) / 86400.0

        # check jdfrac
        if jdfrac > 1.0:
            jd = jd + math.floor(jdfrac)
            jdfrac = jdfrac - math.floor(jdfrac)

        #  - 0.5 * sign(100.0 * yr + mon - 190002.5) + 0.5
        return jd, jdfrac


# ------------------------------------------------------------------------------
#
#                           function invjday
#
#  this function finds the year, month, day, hour, minute and second
#    given the julian date. tu can be ut1, tdt, tdb, etc.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    jd          - julian date                    days from 4713 bc
#
#  outputs       :
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - hour                           0 .. 23
#    min         - minute                         0 .. 59
#    sec         - second                         0.0 .. 59.999
#
#  locals        :
#    days        - day of year plus fractional
#                  portion of a day               days
#    tu          - julian centuries from 0 h
#                  jan 0, 1900
#    temp        - temporary real values
#    leapyrs     - number of leap years from 1900
#
#  coupling      :
#    days2mdhms  - finds month, day, hour, minute and second given days and year
#
#  references    :
#    vallado       2007, 208, alg 22, ex 3-13
#
# [year, mon, day, hr, min, sec] = invjday (jd, jdfrac)
# -----------------------------------------------------------------------------

def invjday (jd, jdfrac):

     # check jdfrac for multiple days
     if (abs(jdfrac) >= 1.0):
         jd = jd + math.floor(jdfrac)
         jdfrac = jdfrac - math.floor(jdfrac)

     # check for fraction of a day included in the jd
     dt = jd - math.floor(jd) - 0.5
     if (abs(dt) > 0.00000001):
         jd = jd - dt
         jdfrac = jdfrac + dt

     # ----------------- find year and days of the year ---------------
     temp = jd - 2415019.5
     tu = temp / 365.25
     year = 1900 + math.floor(tu)
     leapyrs = math.floor((year-1901)*0.25)
     days = math.floor(temp - ((year-1900)*365.0 + leapyrs))

     # ------------ check for case of beginning of a year -------------
     if (days + jdfrac < 1.0):
         year = year - 1
         leapyrs = math.floor((year-1901)*0.25)
         days = math.floor(temp - ((year-1900)*365.0 + leapyrs))

     # ------------------- find remaining data  -----------------------
     # now add the daily time in to preserve accuracy
     mon, day, hr, min, sec = days2mdh(year, days + jdfrac)
     return year, mon, day, hr, min, sec


# ------------------------------------------------------------------------------
#
#                           function days2mdh
#
#  this function converts the day of the year, days, to the equivalent month
#    day, hour, minute and second.
#
#  author        : david vallado                  719-573-2600   22 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    year        - year                           1900 .. 2100
#    days        - julian day of the year         1.0  .. 366.0
#
#  outputs       :
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - hour                           0 .. 23
#    minute      - minute                         0 .. 59
#    sec         - second                         0.0 .. 59.999
#
#  locals        :
#    dayofyr     - day of year
#    temp        - temporary extended values
#    inttemp     - temporary integer value
#    i           - index
#    lmonth(12)  - integer array containing the number of days per month
#
#  coupling      :
#    none.
#
# [mon, day, hr, minute, sec] = days2mdh (year, days)
# -----------------------------------------------------------------------------

def days2mdh (year, days):

        # --------------- set up array of days in month  --------------
        lmonth = np.zeros(12)
        for i in range(12):
            lmonth[i] = 31
            if i == 1:
                lmonth[i] = 28
            if (i == 3 or i == 5 or i == 8 or i == 10):
                lmonth[i] = 30

        dayofyr = math.floor(days)

        # ----------------- find month and day of month ---------------
        if (np.fmod(year-1900, 4) == 0):
            lmonth[1] = 29

        i = 0
        inttemp = 0
        while (dayofyr > inttemp + lmonth[i]) and (i < 11):
            inttemp = inttemp + lmonth[i]
            i = i+1

        mon = i+1
        day = dayofyr - inttemp

        # ----------------- find hours minutes and seconds ------------
        temp = (days - dayofyr)*24.0
        hr = np.fix(temp)
        temp = (temp-hr) * 60.0
        minute = np.fix(temp)
        sec = (temp-minute) * 60.0

        return mon, day, hr, minute, sec


# -----------------------------------------------------------------------------
#
#                           function lstime
#
#  this function finds the local sidereal time at a given location.
#
#  author        : david vallado                  719-573-2600   27 may 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    lon         - site longitude (west -)        -2pi to 2pi rad
#    jd          - julian date                    days from 4713 bc
#
#  outputs       :
#    lst         - local sidereal time            0.0 to 2pi rad
#    gst         - greenwich sidereal time        0.0 to 2pi rad
#
#  locals        :
#    none.
#
#  coupling      :
#    gstime        finds the greenwich sidereal time
#
#  references    :
#    vallado       2007, 194, alg 15, ex 3-5
#
# [lst, gst] = lstime (lon, jd)
# -----------------------------------------------------------------------------

def lstime (lon, jd):
        # ------------------------  implementation   ------------------
        gst = gstime(jd)
        lst = lon + gst

        # ----------------------- check quadrants ---------------------
        lst = np.fmod(lst, twopi)
        if (lst < 0.0):
            lst = lst + twopi
        return lst, gst




# ------------------------------------------------------------------------------
#
#                           function convtime
#
#  this function finds the time parameters and julian century values for inputs
#    of utc or ut1. numerous outputs are found as shown in the local variables.
#    because calucations are in utc, you must include timezone if (you enter a
#    local time, otherwise it should be zero.
#
#  author        : david vallado                  719-573-2600    4 jun 2002
#
#  revisions
#    vallado     - add tcg, tcb, etc                              6 oct 2005
#    vallado     - fix documentation for dut1                     8 oct 2002
#
#  inputs          description                    range / units
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - universal time hour            0 .. 23
#    min         - universal time min             0 .. 59
#    sec         - universal time sec (utc)            0.0  .. 59.999
#    timezone    - offset to utc from local site  0 .. 23 hr
#    dut1        - delta of ut1 - utc             sec
#    dat         - delta of tai - utc             sec
#
#  outputs       :
#    ut1         - universal time                 sec
#    tut1        - julian centuries of ut1
#    jdut1       - julian date (days only)           days from 4713 bc
#    jdut1Frac   - julian date (fraction of a day)   days from 0 hr of the day
#    utc         - coordinated universal time     sec
#    tai         - atomic time                    sec
#    tdt         - terrestrial dynamical time     sec
#    ttdt        - julian centuries of tdt
#    jdtt        - julian date (days only)           days from 4713 bc
#    jdttFrac    - julian date (fraction of a day)   days from 0 hr of the day
#    tdb         - terrestrial barycentric time   sec
#    ttdb        - julian centuries of tdb
#    jdtdb       - julian date of tdb             days from 4713 bc
#    tcb         - celestial barycentric time     sec
#    tcg         - celestial geocentric time      sec
#    jdtdb       - julian date (days only)           days from 4713 bc
#    jdtdbFrac   - julian date (fraction of a day)   days from 0 hr of the day
#
#  locals        :
#    hrtemp      - temporary hours                hr
#    mintemp     - temporary minutes              min
#    sectemp     - temporary seconds              sec
#    localhr     - difference to local time       hr
#    jd          - julian date of request         days from 4713 bc
#    me          - mean anomaly of the earth      rad
#
#  coupling      :
#    hms_2_sec   - conversion between hr-min-sec .and. seconds
#    jday        - find the julian date
#
#  references    :
#    vallado       2007, 201, alg 16, ex 3-7
#
#  leave out for now...
#  , tcg, jdtcg, jdtcgfrac, tcb, jdtcb, jdtcbfrac
# [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, ...
#  tdb, ttdb, jdtdb, jdtdbfrac] ...
# = convtime (year, mon, day, hr, min, sec, timezone, dut1, dat)
# ------------------------------------------------------------------------------


def convtime(year, mon, day, hr, min, sec, timezone, dut1, dat):

        # ------------------------  implementation   ------------------
        jd, jdfrac = jday(year, mon, day, hr + timezone, min, sec)
        mjd = jd + jdfrac - 2400000.5
        mfme = hr*60.0 + min + sec/60.0

        # ------------------ start if (ut1 is known ------------------
        localhr = timezone + hr
        utc = hms2sec(localhr, min, sec)

        ut1 = utc + dut1
        hrtemp, mintemp, sectemp = sec2hms(ut1)
        jdut1, jdut1frac = jday(year, mon, day, hrtemp, mintemp, sectemp)
        tut1 = (jdut1 + jdut1frac - 2451545.0) / 36525.0

        tai = utc + dat
        hrtemp, mintemp, sectemp = sec2hms(tai)
        jdtai, jdtaifrac = jday(year, mon, day, hrtemp, mintemp, sectemp)

        tt = tai + 32.184   # sec
        hrtemp, mintemp, sectemp = sec2hms(tt)
        jdtt, jdttfrac = jday(year, mon, day, hrtemp, mintemp, sectemp)
        ttt = (jdtt+jdttfrac - 2451545.0) / 36525.0

#%%%%%%%%%%%%%%%%%%% tdb
        '''
        % vallado approach (extra digits)
        %         me = 357.5277233  + 35999.05034 *ttt
        %         me = mod(me, 360.0 )
        %         me = me * deg2rad
        %         tdb = tt + 0.001658  * sin(me) + 0.00001385 *sin(2.0 *me)
        %         [hrtemp, mintemp, sectemp] = sec2hms(tdb)
        %         [jdtdb, jdtdbfrac] = jday(year, mon, day, hrtemp, mintemp, sectemp)
        %         ttdb = (jdtdb + jdtdbfrac - 2451545.0 )/ 36525.0
        %         fprintf(1, 'book tdb %8.6f ttdb  %16.12f jdtdb  %18.11f %18.11f \n', tdb, ttdb, jdtdb, jdtdbfrac)
        % std approach (digits)
        %         me = 357.53  + 0.9856003 * (jdtt - 2451545.0)
        %         me = mod(me, 360.0 )
        %         me = me * deg2rad
        %         tdb1 = tt + 0.001658  * sin(me) + 0.000014 *sin(2.0 *me)
        %         [hrtemp, mintemp, sectemp] = sec2hms(tdb1)
        %         [jdtdb1, jdtdb1frac] = jday(year, mon, day, hrtemp, mintemp, sectemp)
        %         ttdb1 = (jdtdb1 + jdtdb1frac - 2451545.0 )/ 36525.0
        %         fprintf(1, 'std  tdb %8.6f ttdb  %16.12f jdtdb  %18.11f %18.11f \n', tdb1, ttdb1, jdtdb1, jdtdb1frac)
        % ast alm approach (2012) bradley email
        '''


        me = 357.53  + 0.98560028 * (jdtt - 2451545.0)
        me = np.remainder(me, 360.0 ) ###not quite equivalent to matlab mod command, but close
        me = me * deg2rad
        dlje = 246.11 + 0.90251792*(jdtt - 2451545.0)
        tdb2 = tt + 0.001657  * math.sin(me) + 0.000022 *math.sin(dlje)
        [hrtemp, mintemp, sectemp] = sec2hms(tdb2)
        [jdtdb2, jdtdb2frac] = jday(year, mon, day, hrtemp, mintemp, sectemp)
        ttdb2 = (jdtdb2 + jdtdb2frac - 2451545.0 )/ 36525.0
 #       fprintf(1, 'asta tdb %8.6f ttdb  %16.12f jdtdb  %18.11f %18.11f \n', tdb2, ttdb2, jdtdb2, jdtdb2frac)
# usno circular approach
        tdb = tt + 0.001657*math.sin(628.3076*ttt+6.2401) \
               + 0.000022*math.sin(575.3385*ttt+4.2970) \
               + 0.000014*math.sin(1256.6152*ttt+6.1969) \
               + 0.000005*math.sin(606.9777*ttt+4.0212) \
               + 0.000005*math.sin(52.9691*ttt+0.4444) \
               + 0.000002*math.sin(21.3299*ttt+5.5431) \
               + 0.000010*ttt*math.sin(628.3076*ttt+4.2490)  # USNO circ (14)
        [hrtemp, mintemp, sectemp] = sec2hms(tdb)
        [jdtdb, jdtdbfrac] = jday(year, mon, day, hrtemp, mintemp, sectemp)
        ttdb = (jdtdb + jdtdbfrac - 2451545.0 )/ 36525.0

#        fprintf(1, 'usno tdb %8.6f ttdb  %16.12f jdtdb  %18.11f %18.11f \n', tdb, ttdb, jdtdb, jdtdbfrac)
        [h, m, s] = sec2hms(tdb)
#        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s)

        #
#%%%%%%%%%%%%%%%%%%% tcg
# approx with tai
        tcg = tt + 6.969290134e-10*(jdtai - 2443144.5003725)*86400.0  # AAS 05-352 (10) and IERS TN (104)
        [hrtemp, mintemp, sectemp] = sec2hms(tcg)
        [jdtcg, jdtcgfrac] = jday(year, mon, day, hrtemp, mintemp, sectemp)
        tt2 = tcg-6.969290134e-10*(jdtcg+jdtcgfrac-2443144.5003725)*86400.0

  #      fprintf(1, 'tcg %8.6f jdtcg  %18.11f ', tcg, jdtcg)
        [h, m, s] = sec2hms(tcg)
#        fprintf(1, 'hms %3i %3i %8.6f \n', h, m, s)

        '''
        % binomial approach with days
        %        lg = 6.969290134e-10*86400.0
        %        tcg1 = tt + (jdtt - 2443144.5003725)*(lg + lg*lg + lg*lg*lg)
        % days from 77
        %        jdttx = jday(year, mon, day, 0, 0, 0.0)
        %        ttx = tt/86400.0 + jdttx-2443144.5003725  % days from the 1977 epoch
        %        tcg2 = (jdttx - 6.969290134e-10*2443144.5003725) / (1.0 - 6.969290134e-10) % days
        %        tcg2 = (tcg2 - jdttx)*86400*86400
        % sec from 77
        %        ttx = tt + (jdttx-2443144.5003725)*86400.0  % s from the 1977 epoch
        %        tcg3 = ttx / (1.0 - 6.969290134e-10) % s
        %        tcg3 = tcg3 -(jdttx-2443144.5003725)*86400.0
        % check with tcg
        %        tcg4 = tt + 6.969290134e-10*(jdtcg - 2443144.5003725)*86400.0  % AAS 05-352 (10) and IERS TN (104)
        %        [hrtemp, mintemp, sectemp] = sec2hms(tcg4)
        %        jdtcg4 = jday(year, mon, day, hrtemp, mintemp, sectemp)
        %        tt2 = tcg4-6.969290134e-10*(jdtcg4-2443144.5003725)*86400.0
        %        difchk = tt2-tt
        '''


        tcbmtdb = -1.55051976772e-8*(jdtai+jdtaifrac - 2443144.5003725)*86400.0 - 6.55e-5  # sec, value for de405 AAS 05-352 (10) and IERS TN (104)?
        tcb = tdb + tcbmtdb
        [hrtemp, mintemp, sectemp] = sec2hms(tcb)
        [jdtcb, jdtcbfrac] = jday(year, mon, day, hrtemp, mintemp, sectemp)
        ttcb = (jdtcb + jdtcbfrac - 2451545.0 )/ 36525.0
#        fprintf(1, '     tcb %8.6f ttcb  %16.12f jdtcb  %18.11f %18.11f \n', tcb, ttcb, jdtcb, jdtcbfrac)


        return ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, \
          tdb, ttdb, jdtdb, jdtdbfrac




# -----------------------------------------------------------------------------
#
#                           function gstime
#
#  this function finds the greenwich sidereal time (iau-82).
#
#  author        : david vallado                  719-573-2600    7 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    jdut1       - julian date of ut1             days from 4713 bc
#
#  outputs       :
#    gst         - greenwich sidereal time        0 to 2pi rad
#
#  locals        :
#    temp        - temporary variable for reals   rad
#    tut1        - julian centuries from the
#                  jan 1, 2000 12 h epoch (ut1)
#
#  coupling      :
#
#  references    :
#    vallado       2007, 193, Eq 3-43
#
# gst = gstime(jdut1)
# -----------------------------------------------------------------------------


def gstime(jdut1):

  # ------------------------  implementation   ------------------
  tut1 = (jdut1 - 2451545.0) / 36525.0

  temp = - 6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1  \
         + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 67310.54841

  # 360/86400 = 1/240, to deg, to rad
  temp = np.fmod(temp*deg2rad/240.0, twopi)

  # ------------------------ check quadrants --------------------
  if (temp < 0.0):
      temp = temp + twopi

  gst = temp

  return gst


# -----------------------------------------------------------------------------
#
#                           function gstime0
#
#  this function finds the greenwich sidereal time at the beginning of a year.
#    this formula is derived from the astronomical almanac and is good only
#    0 hr ut1, jan 1 of a year.
#
#  author        : david vallado                  719-573-2600   22 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    year        - year                           1998, 1999, etc.
#
#  outputs       :
#    gst  0      - greenwich sidereal time        0 to 2pi rad
#
#  locals        :
#    jd          - julian date                    days from 4713 bc
#    temp        - temporary variable for reals   rad
#    tut1        - julian centuries from the
#                  jan 1, 2000 12 h epoch (ut1)
#
#  coupling      :
#
#  references    :
#    vallado        2007, 195, Eq 3-46
#
# gst0 = gstime0(year)
# -----------------------------------------------------------------------------

def gstime0(year):
    # ------------------------  implementation   ------------------
    jd = 367.0 * year  \
         - math.floor((7 * (year + math.floor(10 / 12.0))) * 0.25)   \
         + math.floor(275 / 9.0) \
         + 1721014.5

    tut1 = (jd - 2451545.0) / 36525.0

    temp = - 6.2e-6 * tut1 * tut1 * tut1   \
         + 0.093104 * tut1 * tut1  \
         + (876600.0  * 3600.0 + 8640184.812866) * tut1  \
         + 67310.54841

    # ------------------------ check quadrants --------------------
    temp = np.fmod(temp, twopi)
    if (temp < 0.0):
       temp = temp + twopi

    gst0 = temp
    return gst0





# ----------------------------------------------------------------------------
#
#                           function sidereal
#
#  this function calulates the transformation matrix that accounts for the
#    effects of sidereal time. Notice that deltaspi should not be mod'ed to a
#    positive number because it is multiplied rather than used in a
#    trigonometric argument.
#
#  author        : david vallado                  719-573-2600   25 jun 2002
#
#  revisions
#    vallado     - fix units on kinematic terms                   5 sep 2002
#    vallado     - add terms                                     30 sep 2002
#    vallado     - consolidate with iau 2000                     14 feb 2005
#
#  inputs          description                    range / units
#    jdut1       - julian centuries of ut1        days
#    deltapsi    - nutation angle                 rad
#    meaneps     - mean obliquity of the ecliptic rad
#    omega       - long of asc node of moon       rad
#    lod         - length of day                  sec
#    eqeterms    - terms for ast calculation      0, 2
#
#  outputs       :
#    st          - transformation matrix for pef - tod
#    stdot       - transformation matrix for pef - tod rate
#
#  locals        :
#    gmst        - mean greenwich sidereal time   0 to 2pi rad
#    ast         - apparent gmst                   0 to 2pi rad
#    hr          - hour                           hr
#    min         - minutes                        min
#    sec         - seconds                        sec
#    temp        - temporary vector
#    tempval     - temporary variable
#
#  coupling      :
#
#  references    :
#    vallado       2013, 223-224
#
# [st, stdot] = sidereal (jdut1, deltapsi, meaneps, omega, lod, eqeterms)
# ----------------------------------------------------------------------------


def sidereal (jdut1, deltapsi, meaneps, omega, lod, eqeterms):

        # ------------------------ find gmst --------------------------
        gmst = gstime(jdut1)

        # ------------------------ find mean ast ----------------------
        # after 1997, kinematic terms apply as well as gemoetric in eqe
        if (jdut1 > 2450449.5) and (eqeterms > 0):
            ast = gmst + deltapsi* math.cos(meaneps) \
                + 0.00264*math.pi /(3600*180)*math.sin(omega) \
                + 0.000063*math.pi /(3600*180)*math.sin(2.0 *omega)
        else:
            ast = gmst + deltapsi* math.cos(meaneps)

        ast = np.fmod (ast, 2.0*math.pi)
        thetasa = earthrot * (1.0  - lod/86400.0)
        omegaearth = thetasa

#print('st gmst %11.8f ast %11.8f ome  %11.8f \n', gmst*180/math.pi, ast*180/math.pi, omegaearth*180/math.pi)

        st = np.zeros((3, 3))
        st[0, 0] = math.cos(ast)
        st[0, 1] = -math.sin(ast)
        st[0, 2] = 0.0
        st[1, 0] = math.sin(ast)
        st[1, 1] = math.cos(ast)
        st[1, 2] = 0.0
        st[2, 0] = 0.0
        st[2, 1] = 0.0
        st[2, 2] = 1.0

        # compute sidereal time rate matrix
        stdot = np.zeros((3, 3))
        stdot[0, 0] = -omegaearth * math.sin(ast)
        stdot[0, 1] = -omegaearth * math.cos(ast)
        stdot[0, 2] = 0.0
        stdot[1, 0] = omegaearth * math.cos(ast)
        stdot[1, 1] = -omegaearth * math.sin(ast)
        stdot[1, 2] = 0.0
        stdot[2, 0] = 0.0
        stdot[2, 1] = 0.0
        stdot[2, 2] = 0.0

        return st, stdot


# -----------------------------------------------------------------------------
#
#                           function jd2sse.m
#
#  this function finds the seconds since epoch (1 Jan 2000) given the julian date
#
#  author        : david vallado                  719-573-2600   12 dec 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    jd          - julian date                    days from 4713 bc
#
#  outputs       :
#    sse         - seconds since epoch 1 jan 2000
#
#  locals        :
#    none.
#
#  coupling      :
#    none.
#
#  references    :
#    none.
#
# sse = jd2sse(jd)
# -----------------------------------------------------------------------------

def jd2sse(jd=None):
    # ------------------------  implementation   ------------------
    sse = (jd - 2451544.5) * 86400.0
    return sse

# ------------------------------------------------------------------------------
#
#                           function finddays
#
#  this function finds the fractional days through a year given the year,
#    month, day, hour, minute and second.
#
#  author        : david vallado                  719-573-2600   22 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    year        - year                           1900 .. 2100
#    mon         - month                          1 .. 12
#    day         - day                            1 .. 28, 29, 30, 31
#    hr          - hour                           0 .. 23
#    min         - minute                         0 .. 59
#    sec         - second                         0.0 .. 59.999
#
#  outputs       :
#    days        - day of year plus fraction of a
#                    day                          days
#
#  locals        :
#    lmonth      - length of months of year
#    i           - index
#
#  coupling      :
#    none.
#
#  references    :
#    vallado       2007, 207, ex 3-12
#
# [days] = finddays (year, month, day, hr, min, sec)
# -----------------------------------------------------------------------------

def finddays(year=None, month=None, day=None, hr=None, min=None, sec=None):
    lmonth = np.zeros(12)
    for i in range(12):
        lmonth[i] = 31
        if i == 1:
            lmonth[i] = 28
        if (i == 3) or (i == 5) or (i == 8) or (i == 10):
            lmonth[i] = 30

    if (np.fmod(year, 4) == 0):
        lmonth[1] = 29
        if ((np.fmod(year, 100) == 0) and (np.fmod(year, 400) != 0)):
            lmonth[1] = 28

    i = 1
    days = 0.0
    while ((i < month) and (i < 12)):

        days = days + lmonth(i)
        i = i + 1


    days = days + day + hr / 24.0 + min / 1440.0 + sec / 86400.0
    return days

##############################################################################################################
##############################################################################################################
##############################################################################################################

if __name__ == '__main__':

    year = 2004
    mon = 4
    day = 6
    hr = 7
    min = 51
    sec = 28.386009

    dut1 = -0.4399619  # sec
    dat = 32         # sec
    timezone =0
    opt = 'c' # specify the iau00 cio approach

    print('input data \n\n')
    print(' year %5i '% year)
    print(' mon %4i '% mon)
    print(' day %3i '% day)
    print(' %3i:%2i:%8.6f\n'% (hr, min, sec))
    print(' dut1 %8.6f s'% dut1)
    print(' dat %3i s'% dat)
    print('units are km and km/s and km/s2\n')

    # -------- convtime    - convert time from utc to all the others
    print('convtime results\n')
    [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, tdb, ttdb, jdtdb, jdtdbfrac ] \
        = convtime (year, mon, day, hr, min, sec, timezone, dut1, dat)
    print('ut1 %8.6f tut1 %16.12f jdut1 %18.11f '% (ut1, tut1, jdut1+jdut1frac))


    jd, jdfrac = jdayall(year, mon, day, hr, min, sec, 'j')
    print("jdayall returned: ", jd, jdfrac)


    jd, jdfrac = jdayall(year, mon, day, hr, min, sec, 'g')
    print("jdayall returned: ", jd, jdfrac)


    jd, jdfrac = jday(year, mon, day, hr, min, sec)
    print("jday returned: ", jd, jdfrac)

    jdut1 = jd
    deltapsi = math.pi
    meaneps = math.pi
    omega = math.pi
    lod = 3600*24
    eqeterms = 1

    st, stdot = sidereal (jdut1, deltapsi, meaneps, omega, lod, eqeterms)
    print("sidereal returned:")
    pp(st)
    pp(stdot)




























