import math
import os
import numpy as np
from pprint import pprint as pp
from space_constants import *
import spacemath_utils as smu


# -----------------------------------------------------------------------------
#
#                           function getintday
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

def getintday(daystr: str):
    """this function finds the integer equivalent of the 3 character string
    representation of the day of the week.

    Parameters
    ----------
    daystr : str
        name of day: 'sun', 'mon' etc.

    Returns
    -------
    dayn : int
        day integer equivalent: 1 to 7
    """
    daytitle = ['sun', 'mon', 'tue', 'wed', 'thr', 'fri', 'sat']
    dayn = daytitle.index(daystr) + 1
    return dayn

# -----------------------------------------------------------------------------
#
#                           function dayofweek
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

def dayofweek(jd: float):
    """this function finds the day of the week. integers are used for the days,
    1 = 'sun', 2 = 'mon', ... 7 = 'sat'.

    Parameters
    ----------
    jd : float
        julian date: days from 4713 bc

    Returns
    -------
    dayofweek: int
        day of the week: 1 to 7
    """
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


def getintmon(monstr: str):
    """this function finds the integer equivalent of the 3 character string
    representation of month.

    Parameters
    ----------
    monstr : str
        3 character month name: 'jan', 'feb', etc.

    Returns
    -------
    mon: int
        integer month equivalent: 1 to 12
    """
    monthtitle = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug',
                  'sep', 'oct', 'nov', 'dec']
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
#    hr          - hours                          0 .. 23
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


def hms2sec(hr: int, min: int, sec: float):
    """this function converts hours, minutes and seconds into seconds from the
    beginning of the day.

    Parameters
    ----------
    hr : int
        hours: 0 to 23
    min : int
        minutes: 0 to 59
    sec : float
        seconds: 0 to 59.99

    Returns
    -------
    utsec: float
        seconds from start of day: 0 to 86400
    """
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
#    hr          - hours                          0 .. 23
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


def sec2hms(utsec: float):
    """this function converts seconds from the beginning of the day into hours,
    minutes and seconds.

    Parameters
    ----------
    utsec : float
        seconds from beginning of day: 0 to 86400

    Returns
    -------
    hr: int
        hours: 0 to 23
    min: int
        minutes: 0 to 59
    sec: float
        seconds: 0 to 59.99
    """
    temp = utsec / 3600.0
    hr = np.fix(temp)
    min = np.fix((temp - hr) * 60.0)
    sec = (temp - hr - min/60.0) * 3600.0
    return hr, min, sec



# -----------------------------------------------------------------------------
#
#                           function jday
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


def jday(yr: int, mon: int, day: int, hr: int, min: int, sec: float):
    """this function finds the julian date given the year, month, day, and
    time.

    Parameters
    ----------
    yr : int
        years: 1900 to 2100
    mon : int
        months: 1 to 12
    day : int
        days: 1 to 31
    hr : int
        hours: 0 to 23
    min : int
        minutes: 0 to 59
    sec : float
        seconds: 0 to 59.99

    Returns
    -------
    jd : float
        julian date: days from 4713 bc
    jdfrac: float
        julian date fraction of a day: 0.0 to 1.0
    """

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
#    jdfrac      - julian date fraction of a day  0.0 to 1.0
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

def invjday (jd: float, jdfrac: float = 0):
    """this function finds the year, month, day, hour, minute and second
    given the julian date. tu can be ut1, tdt, tdb, etc.

    Parameters
    ----------
    jd : float
        julian date: days from 4713 bc
    jdfrac: float, optional
        julian date fraction of a day: 0.0 to 1.0 , default 0

    Returns
    -------
    yr : int
        years: 1900 to 2100
    mon : int
        months: 1 to 12
    day : int
        days: 1 to 31
    hr : int
        hours: 0 to 23
    min : int
        minutes: 0 to 59
    sec : float
        seconds: 0 to 59.99
    """

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

def jdayall(year: int, mon: int, day: int, hr: int, min: int, sec: float,
            whichtype: str):
    """this function finds the julian date given the year, month, day, and time.
    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.

    Parameters
    ----------
    year : int
        year: 1900 to 2100
    mon : int
        month: 1 to 12
    day : int
        day: 1 to 31
    hr : int
        hour: 0 to 23
    min : int
        minutes: 0 to 59
    sec : float
        seconds: 0 to 59.99
    whichtype : str
        Julian or gegorian: 'j' or 'g'

    Returns
    -------
    jd: float
        julian date: days from 4713 bc
    """

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
#    min         - minute                         0 .. 59
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

def days2mdh (year: int, days: float):
    """this function converts the day of the year, days, to the equivalent month
    day, hour, minute and second.

    Parameters
    ----------
    year : int
        year: 1900 to 2100
    days : float
        julian day of the year: 0 to 366.0

    Returns
    -------
    mon : int
        month: 1 to 12
    day : int
        day: 1 to 31
    hr : int
        hour: 0 to 23
    min : int
        minutes: 0 to 59
    sec: float
        seconds: 0 to 59.99
    """

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
    min = np.fix(temp)
    sec = (temp-min) * 60.0
    sec = round(sec, 3)
    if (sec == 60.0):
        sec = 0.0
        min = min + 1


    return mon, day, hr, min, sec


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

def lstime(lon: float, jd: float):
    """this function finds the local sidereal time at a given location.

    Parameters
    ----------
    lon : float
        site longitude (west -): -2pi to 2pi rad
    jd : float
        julian date: days from 4713

    Returns
    -------
    lst : float
        local sidereal time: 0 to 2pi rad
    gst : float
        greenwich sidereal time: 0 to 2pi rad
    """
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
#  inputs          description                              range / units
#    year        - year                                     1900 .. 2100
#    mon         - month                                    1 .. 12
#    day         - day                                      1 .. 28, 29, 30, 31
#    hr          - universal time hour                      0 .. 23
#    min         - universal time min                       0 .. 59
#    sec         - universal time sec (utc)                 0.0  .. 59.999
#    timezone    - offset to utc from local site            0 .. 23 hr
#    dut1        - delta of ut1 - utc                       sec
#    dat         - delta of tai - utc                       sec
#
#  outputs       :
#    ut1         - universal time                           sec
#    tut1        - julian centuries of ut1
#    jdut1       - julian date (days only)                  days from 4713 bc
#    jdut1frac   - julian date (fraction of a day)          days from 0 hr of the day
#    utc         - coordinated universal time               sec
#    tai         - atomic time                              sec
#    tt          - terrestrial time                         sec
#    ttt         - julian centuries of tt
#    jdtt        - julian date (days only)                  days from 4713 bc
#    jdttfrac    - julian date (fraction of a day)          days from 0 hr of the day
#    tdb         - terrestrial barycentric dynamical time   sec
#    ttdb        - julian centuries of tdb
#    jdtdb       - julian date of tdb                       days from 4713 bc
#    jdtdb       - julian date (days only)                  days from 4713 bc
#    jdtdbfrac   - julian date (fraction of a day)          days from 0 hr of the day
#
#  locals        :
#    hrtemp      - temporary hours                          hr
#    mintemp     - temporary minutes                        min
#    sectemp     - temporary seconds                        sec
#    localhr     - difference to local time                 hr
#    jd          - julian date of request                   days from 4713 bc
#    me          - mean anomaly of the earth                rad
#
#  coupling      :
#    hms_2_sec   - conversion between hr-min-sec .and. seconds
#    jday        - find the julian date
#
#  references    :
#    vallado       2007, 201, alg 16, ex 3-7
#
#  leave out for now...
#    tcb         - celestial barycentric time               sec
#    tcg         - celestial geocentric time                sec
#    (as well as jdtcg, jdtcgfrac, jdtcb, jdtcbfrac)
#
# [ut1, tut1, jdut1, jdut1frac, utc, tai, tt, ttt, jdtt, jdttfrac, ...
#  tdb, ttdb, jdtdb, jdtdbfrac] ...
# = convtime (year, mon, day, hr, min, sec, timezone, dut1, dat)
# ------------------------------------------------------------------------------


def convtime(year: int, mon: int, day: int, hr: int, min: int, sec: float,
             timezone: int, dut1: float, dat: float):
    """this function finds the time parameters and julian century values for
    inputs of utc or ut1. numerous outputs are found as shown in the local
    variables. because calucations are in utc, you must include timezone if
    (you enter a local time, otherwise it should be zero.

    Parameters
    ----------
    year : int
        year: 1900 to 2100
    mon : int
        month: 1 to 12
    day : int
        day: 1 to 31
    hr : int
        hour: 0 to 23
    min : int
        minutes: 0 to 59
    sec : float
        seconds: 0 to 59.99
    timezone : int
        offset to utc from local site: 0 to 23 hr
    dut1 : float
        delta of ut1 - utc: sec
    dat : float
        delta of tai - utc: sec

    Returns
    --------
    ut1 : float
        universal time: sec
    tut1 : float
        julian centuries of ut1
    jdut1 : float
        julian date \(days only): days from 4713 bc
    jdut1frac : float
        julian date \(fraction of a day): days from 0 hr of the day
    utc : float
        coordinated universal time: sec
    tai : float
        atomic time: sec
    tt : float
        terrestrial time: sec
    ttt : float
        julian centuries of tt
    jdtt : float
        julian date \(days only): days from 4713 bc
    jdttfrac : float
        julian date \(fraction of a day): days from 0 hr of the day
    tdb : float
        terrestrial barycentric dynamical time: sec
    ttdb : float
        julian centuries of tdb
    jdtdb : float
        julian date \(days only): days from 4713 bc
    jdtdbfrac : float
        julian date \(fraction of a day): days from 0 hr of the day
    """

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


    # me = 357.53  + 0.98560028 * (jdtt - 2451545.0)
    # me = np.remainder(me, 360.0 ) ###not quite equivalent to matlab mod command, but close
    # me = me * deg2rad
    # dlje = 246.11 + 0.90251792*(jdtt - 2451545.0)
    # tdb2 = tt + 0.001657  * math.sin(me) + 0.000022 *math.sin(dlje)
    # [hrtemp, mintemp, sectemp] = sec2hms(tdb2)
    # [jdtdb2, jdtdb2frac] = jday(year, mon, day, hrtemp, mintemp, sectemp)
    # ttdb2 = (jdtdb2 + jdtdb2frac - 2451545.0 )/ 36525.0
    # fprintf(1, 'asta tdb %8.6f ttdb  %16.12f jdtdb  %18.11f %18.11f \n', tdb2, ttdb2, jdtdb2, jdtdb2frac)

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


    tcbmtdb = (-1.55051976772e-8 * (jdtai + jdtaifrac - 2443144.5003725)
               * 86400.0 - 6.55e-5)  # sec, value for de405 AAS 05-352 (10) and IERS TN (104)?
    tcb = tdb + tcbmtdb
    [hrtemp, mintemp, sectemp] = sec2hms(tcb)
    [jdtcb, jdtcbfrac] = jday(year, mon, day, hrtemp, mintemp, sectemp)
    ttcb = (jdtcb + jdtcbfrac - 2451545.0 ) / 36525.0
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

def gstime(jdut1: float):
    """this function finds the greenwich sidereal time (iau-82).

    Parameters
    ----------
    jdut1 : float
        julian date of ut1: days from 4713 bc

    Returns
    -------
    gst : float
        greenwich sidreal time: rad
    """

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
#    gst0        - greenwich sidereal time        0 to 2pi rad
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

def gstime0(year: int):
    """this function finds the greenwich sidereal time at the beginning of a
    year. this formula is derived from the astronomical almanac and is good
    only 0 hr ut1, jan 1 of a year.

    Parameters
    ----------
    year : int
        year: 1998, 1999 etc

    Returns
    -------
    gst0 : float
        greenwich sidereal time: 0 to 2pi rad
    """

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


def sidereal(jdut1: float, deltapsi: float, meaneps: float, omega: float,
             lod: float, eqeterms: int):
    """this function calulates the transformation matrix that accounts for the
    effects of sidereal time. Notice that deltaspi should not be mod'ed to a
    positive number because it is multiplied rather than used in a
    trigonometric argument.

    Parameters
    ----------
    jdut1 : float
        julian days of ut1: days from 4713 bc
    deltapsi : float
        nutation angle: rad
    meaneps : float
        mean obliquity of the ecliptic: rad
    omega : float
        longitude of ascending node of moon: rad
    lod : float
        length of day: sec
    eqeterms : int
        terms for ast calculation: 0, 2
    """


    # ------------------------ find gmst --------------------------
    gmst = gstime(jdut1)

    # ------------------------ find mean ast ----------------------
    # after 1997, kinematic terms apply as well as geometric in eqe
    if (jdut1 > 2450449.5) and (eqeterms > 0):
        ast = gmst + deltapsi* math.cos(meaneps) \
            + 0.00264 * arcsec2rad * math.sin(omega) \
            + 0.000063 * arcsec2rad * math.sin(2.0 *omega)
    else:
        ast = gmst + deltapsi* math.cos(meaneps)

    ast = np.fmod(ast, 2.0*math.pi)
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
#                           function jd2sse
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

def jd2sse(jd: float):
    """this function finds the seconds since epoch (1 Jan 2000) given the
    julian date

    Parameters
    ----------
    jd : float
        julian date: days from 4713 bc

    Returns
    -------
    sse: float
        seconds since epoch 1 jan 2000
    """
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
#  inputs          description                              range / units
#    year        - year                                     1900 .. 2100
#    mon         - month                                    1 .. 12
#    day         - day                                      1 .. 28, 29, 30, 31
#    hr          - hour                                     0 .. 23
#    min         - minute                                   0 .. 59
#    sec         - second                                   0.0 .. 59.999
#
#  outputs       :
#    days        - day of year plus fraction of a day       days
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

def finddays(year: int, month: int, day: int, hr: int, min: int, sec: float):
    """this function finds the fractional days through a year given the year,
    month, day, hour, minute and second.

    Parameters
    ----------
    year : int
        years: 1900 to 2100
    month : int
        months: 1 to 12
    day : int
        days: 1 t 31
    hr : int
        hours: 0 to 23
    min : int
        minutes: 0 to 59
    sec : float
        seconds: 0 to 59.99

    Returns
    -------
    days: float
        day of the year + fraction: days
    """
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


def initjplde(infilename: str):
    """this function initializes the jpl planetary ephemeris data. the input
       data files are from processing the ascii files into a text file of sun
       and moon positions.

    Parameters
    ----------
    infilename : str
        name of the data file to use

    Returns
    -------
    jpldearr : dictionary of arrays
        arrays of jplde data records
    jdjpldestart : float
        julian date of the start of the jpldearr data
    """
    #double jdtdb, jdtdbf;
    #string pattern;
    #Int32 i;
    #jdjpldestart = 0.0;
    #jdjpldestartFrac = 0.0;

    # read the whole file at once into lines of an array
    filename = os.path.join(os.path.dirname(__file__), 'data', infilename)

    filedat = np.loadtxt(filename)
    filedim = filedat.shape[0]

    jpldearr = {}
    #load data into x y z arrays
    jpldearr['year'] = filedat[:, 0]
    jpldearr['mon'] = filedat[:, 1]
    jpldearr['day'] = filedat[:, 2]
    jpldearr['rsun1'] = filedat[:, 3]
    jpldearr['rsun2'] = filedat[:, 4]
    jpldearr['rsun3'] = filedat[:, 5]
    jpldearr['rsmag'] = filedat[:, 6]
    jpldearr['rmoon1'] = filedat[:, 8]
    jpldearr['rmoon2'] = filedat[:, 9]
    jpldearr['rmoon3'] = filedat[:, 10]
    for i in range(filedim):
        jdtdb,jdtdbf = jday(jpldearr['year'][i], jpldearr['mon'][i],
                            jpldearr['day'][i], 0, 0, 0.0)
        jpldearr['mjd'][i] = jdtdb + jdtdbf - 2400000.5
    # ---- find epoch date
    jdjpldestart,jdjpldestartFrac = jday(jpldearr['year'][0],
                                         jpldearr['mon'][0],
                                         jpldearr['day'][0], 0, 0, 0.0)
    # doesn't have the -2400000.5 for some reason? -zeg
    # jpldearr['mjd'][0] = jdjpldestart + jdjpldestartFrac
    return jpldearr, jdjpldestart, jdjpldestartFrac

# -----------------------------------------------------------------------------
#
#                           function findjpldeparam
#
#  this routine finds the jplde parameters for a given time. several types of
#  interpolation are available. the cio and iau76 nutation parameters are also
#  read for optimizing the speeed of calculations.
#
#  author        : david vallado                      719-573-2600   12 dec 2005
#
#  inputs          description                               range / units
#    jdtdb         - epoch julian date                     days from 4713 BC
#    jdtdbF        - epoch julian date fraction            day fraction from jdutc
#    interp        - interpolation                        n-none, l-linear, s-spline
#    jpldearr      - array of jplde data records
#    jdjpldestart  - julian date of the start of the jpldearr data (set in initjplde)
#
#
#  outputs       :
#    dut1        - delta ut1 (ut1-utc)                        sec
#    dat         - number of leap seconds                     sec
#    lod         - excess length of day                       sec
#    xp          - x component of polar motion                rad
#    yp          - y component of polar motion                rad
#    ddpsi       - correction to delta psi (iau-76 theory)    rad
#    ddeps       - correction to delta eps (iau-76 theory)    rad
#    dx          - correction to x (cio theory)               rad
#    dy          - correction to y (cio theory)               rad
#    x           - x component of cio                         rad
#    y           - y component of cio                         rad
#    s           -                                            rad
#    deltapsi    - nutation longitude angle                   rad
#    deltaeps    - obliquity of the ecliptic correction       rad
#
#  locals        :
#                -
#
#  coupling      :
#    none        -
#
#  references    :
#    vallado       2013,
#  [rsun, rsmag, rmoon, rmmag] = findjpldeparam( jdtdb, jdtdbF, 'l', jpldearr, jdjpldestart);
# --------------------------------------------------------------------------- */

def findjpldeparam(jdtdb: float, jdtdbF: float, interp: str, jpldearr,
                   jdjpldestart: float):
    """this routine finds the jplde parameters for a given time. several types of
    interpolation are available. the cio and iau76 nutation parameters are also
    read for optimizing the speeed of calculations.

    Parameters
    ----------
    jdtdb : float
        epoch julian date: days from 4713 bc
    jdtdbF : float
        epoch julian date fraction: day fraction from jdutc
    interp : str
        interpolation: 'n' - none, 'l' - linear, 's' spline
    jpldearr : dictionary of arrays
        jplde information from initjplde
    jdjpldestart : float
        julian date start of the jpldearr data set in initjplde

    Returns
    -------
    rsun : ndarray
        position vector of the sun
    rsmag : float
        magnitude of sun position vector
    rmoon : ndarray
        position vector of the moon
    rmmag : float
        magnitude of moon position vector
    """

    #Int32 recnum;
    #Int32 off1, off2;
    #double fixf, jdjpldestarto, mjd, jdb, mfme;
    #jpldedataClass jplderec, nextjplderec;
    #rsun = new double[3];
    #rmoon = new double[3];

    # the ephemerides are centered on jdtdb, but it turns out to be 0.5, or 0000 hrs.
    # check if any whole days in jdF
    jdb = int(math.floor(jdtdb + jdtdbF)) + 0.5

    mfme = (jdtdb + jdtdbF - jdb) * 1440.0
    if (mfme < 0.0):
        mfme = 1440.0 + mfme

    #printf("jdtdb #lf  #lf  #lf  #lf \n ", jdtdb, jdtdbF, jdb, mfme);x[recnum]

    # ---- read data for day of interest
    jdjpldestarto = int(math.floor(jdtdb + jdtdbF - jdjpldestart))
    recnum = int(math.floor(jdjpldestarto)) - 1
    # check for out of bound values
    rsun = np.empty(3)
    rmoon = np.empty(3)
    if ((recnum >= 0) and (recnum <= 51829)):
        # ---- set non-interpolated values
        rsun[0] = jpldearr['rsun1'][recnum]
        rsun[1] = jpldearr['rsun2'][recnum]
        rsun[2] = jpldearr['rsun3'][recnum]
        rsmag = jpldearr['rsmag']
        mjd = jpldearr['mjd'][recnum]
        rmoon[0] = jpldearr['rmoon1'][recnum]
        rmoon[1] = jpldearr['rmoon2'][recnum]
        rmoon[2] = jpldearr['rmoon3'][recnum]
        rmmag = smu.mag(rmoon)
        # ---- find nutation parameters for use in optimizing speed
        # ---- do linear interpolation
        if (interp == 'l'):
            fixf = mfme / 1440.0
            rsun[0] = (jpldearr['rsun1'][recnum]
                       + (jpldearr['rsun1'][recnum + 1]
                          - jpldearr['rsun1'][recnum]) * fixf)
            rsun[1] = (jpldearr['rsun2'][recnum]
                       + (jpldearr['rsun2'][recnum + 1]
                          - jpldearr['rsun2'][recnum]) * fixf)
            rsun[2] = (jpldearr['rsun3'][recnum]
                       + (jpldearr['rsun3'][recnum + 1]
                          - jpldearr['rsun3'][recnum]) * fixf)
            rsmag = smu.mag(rsun)
            rmoon[0] = (jpldearr['rmoon1'][recnum]
                        + (jpldearr['rmoon1'][recnum + 1]
                           - jpldearr['rmoon1'][recnum]) * fixf)
            rmoon[1] = (jpldearr['rmoon2'][recnum]
                        + (jpldearr['rmoon2'][recnum + 1]
                           - jpldearr['rmoon2'][recnum]) * fixf)
            rmoon[2] = (jpldearr['rmoon3'][recnum]
                        + (jpldearr['rmoon3'][recnum + 1]
                           - jpldearr['rmoon3'][recnum]) * fixf)
            rmmag = smu.mag(rmoon)
            #printf("sunm #i rsmag #lf fixf #lf n #lf nxt #lf \n", recnum, rsmag, fixf, jpldearr['rsun'](1), jpldearr['rsun'](1));
            #printf("recnum l #i fixf #lf #lf rsun #lf #lf #lf \n", recnum, fixf, jpldearr['rsun'](1), rsun(1), rsuny, rsunz);
        # ---- do spline interpolations
        if (interp == 's'):
            off1 = 1
            off2 = 2
            fixf = mfme / 1440.0
            # setup so the interval is in between points 2 and 3
            rsun[0] = smu.cubicinterp(jpldearr['rsun1'][recnum - off1],
                                      jpldearr['rsun1'][recnum],
                                      jpldearr['rsun1'][recnum + off1],
                                      jpldearr['rsun1'][recnum + off2],
                                      jpldearr['mjd'][recnum - off1],
                                      jpldearr['mjd'][recnum],
                                      jpldearr['mjd'][recnum + off1],
                                      jpldearr['mjd'][recnum + off2],
                                      jpldearr['mjd'][recnum] + fixf)
            rsun[1] = smu.cubicinterp(jpldearr['rsun2'][recnum - off1],
                                      jpldearr['rsun2'][recnum],
                                      jpldearr['rsun2'][recnum + off1],
                                      jpldearr['rsun2'][recnum + off2],
                                      jpldearr['mjd'][recnum - off1],
                                      jpldearr['mjd'][recnum],
                                      jpldearr['mjd'][recnum + off1],
                                      jpldearr['mjd'][recnum + off2],
                                      jpldearr['mjd'][recnum] + fixf)
            rsun[2] = smu.cubicinterp(jpldearr['rsun3'][recnum - off1],
                                      jpldearr['rsun3'][recnum],
                                      jpldearr['rsun3'][recnum + off1],
                                      jpldearr['rsun3'][recnum + off2],
                                      jpldearr['mjd'][recnum - off1],
                                      jpldearr['mjd'][recnum],
                                      jpldearr['mjd'][recnum + off1],
                                      jpldearr['mjd'][recnum + off2],
                                      jpldearr['mjd'][recnum] + fixf)
            rsmag = smu.mag(rsun)
            rmoon[0] = smu.cubicinterp(jpldearr['rmoon1'][recnum - off1],
                                       jpldearr['rmoon1'][recnum],
                                       jpldearr['rmoon1'][recnum + off1],
                                       jpldearr['rmoon1'][recnum + off2],
                                       jpldearr['mjd'][recnum - off1],
                                       jpldearr['mjd'][recnum],
                                       jpldearr['mjd'][recnum + off1],
                                       jpldearr['mjd'][recnum + off2],
                                       jpldearr['mjd'][recnum] + fixf)
            rmoon[1] = smu.cubicinterp(jpldearr['rmoon2'][recnum - off1],
                                       jpldearr['rmoon2'][recnum],
                                       jpldearr['rmoon2'][recnum + off1],
                                       jpldearr['rmoon2'][recnum + off2],
                                       jpldearr['mjd'][recnum - off1],
                                       jpldearr['mjd'][recnum],
                                       jpldearr['mjd'][recnum + off1],
                                       jpldearr['mjd'][recnum + off2],
                                       jpldearr['mjd'][recnum] + fixf)
            rmoon[2] = smu.cubicinterp(jpldearr['rmoon3'][recnum - off1],
                                       jpldearr['rmoon3'][recnum],
                                       jpldearr['rmoon3'][recnum + off1],
                                       jpldearr['rmoon3'][recnum + off2],
                                       jpldearr['mjd'][recnum - off1],
                                       jpldearr['mjd'][recnum],
                                       jpldearr['mjd'][recnum + off1],
                                       jpldearr['mjd'][recnum + off2],
                                       jpldearr['mjd'][recnum] + fixf)
            rmmag = smu.mag(rmoon)
            #printf("recnum s #i mfme #lf days rsun #lf #lf #lf \n", recnum, mfme, rsunx, rsuny, rsunz);
            #printf(" #lf #lf #lf #lf \n", jpldearr(recnum - off2).mjd, jpldearr(recnum - off1.mjd, jpldearr[recnum].mjd, jpldearr[recnum + off1].mjd);
        # set default values
    else:
        rsun[0] = 0.0
        rsun[1] = 0.0
        rsun[2] = 0.0
        rsmag = 0.0
        rmoon[0] = 0.0
        rmoon[1] = 0.0
        rmoon[2] = 0.0
        rmmag = 0.0

    #  findjpldeparam
    return rsun, rsmag, rmoon, rmmag

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




























