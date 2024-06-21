import os
import numpy as np
import math
import orbit_utils as obu
from space_constants import *
import space_conversions as sc
import spacemath_utils as smu
import spacetime_utils as stu

# script testmat.py
#
# This script tests the SGP4 propagator.

# Author:
#   Jeff Beck
#   beckja@alumni.lehigh.edu

# Version Info:
#   1.0 (051019) - Initial version from Vallado C++ version.
#   1.0 (aug 14, 2006) - update for paper
#   2.0 (apr 2, 2007) - update for manual operations
#   3.0 (3 jul, 2008) - update for opsmode operation afspc or improved
#   3.1 (2 dec, 2008) - fix tsince/1440.0 in jd update

# sgp4fix consolidate call to getgravconst in sgp4init
# these are set in sgp4init
# global tumin mu radiusearthkm xke j2 j3 j4 j3oj2

directory = os.path.join(os.path.dirname(__file__), "data")
outdir = os.path.join(os.path.dirname(__file__), 'testoutput')
print('output directory set to: %s, change in testmat.py if needed \n'
      % (directory))
# global opsmode

# add operation smode for afspc (a) or improved (i)
# opsmode= input('input opsmode afspc (a), improved i ', 's');
# this is the standard method of operation
opsmode = 'a'
#         //typerun = 'c' compare 1 year of full satcat data
#         //typerun = 'v' verification run, requires modified elm file with
#         //typerun = 'm' maunual operation- either mfe, epoch, or dayof yr
#         //              start stop and delta times


# options are c, v, m, u
# m not implemented -zeg
typerun = 'u'
if (typerun == 'm'):
    # input mfe, epoch (YMDHMS), or dayofyr approach, m, e, d: ', 's')
    typeinput = 'e'
else:
    typeinput = 'e'

# whichconst = input('input constants 721, (72), 84 ');
# this is the standard method of operation
whichconst = 72
# placeholder
ateme = np.array([0.0, 0.0, 0.0])
#         // ---------------- setup files for operation ------------------
#         // input 2-line element set file
infilename = os.path.join(directory, "TestTLE.dat")
infile = open(infilename, 'r')
if (infile == - 1):
    print('Failed to open file: %s\n' % (infilename))
    #return reci, veci, aeci

if (typerun == 'c'):
    outfile = open(os.path.join(outdir, 'tmatall.out'), 'wt')
elif (typerun == 'v'):
    outfile = open(os.path.join(outdir, 'tmatver.out'), 'wt')
else:
    outfile = open(os.path.join(outdir, 'tmat.out'), 'wt')

if (outfile == - 1):
    print('Failed to open outfile in %s\n' % (outdir))
    #return reci, veci, aeci

idebug = True
dbgfile = -1
satrec = {}
#        // ----------------- test simple propagation -------------------
longstr2 = 'Not EOF \n'
while longstr2[-1] == '\n':

    longstr1 = infile.readline()
    while not longstr1[0].isnumeric():
        longstr1 = infile.readline()
    longstr2 = infile.readline()

# while (not feof(infile) ):
#     longstr1 = fgets(infile, 130)
#     while ((longstr1(1) == '#') and (feof(infile) == 0)):
#         longstr1 = fgets(infile, 130)
#     if (feof(infile) == 0):
#         longstr2 = fgets(infile, 130)

    # sgp4fix additional parameters to store from the TLE
    satrec['classification'] = 'U'
    satrec['intldesg'] = '        '
    satrec['ephtype'] = 0
    satrec['elnum'] = 0
    satrec['revnum'] = 0
    if idebug:
        catno = longstr1[2:7]
        dbgfile = open(os.path.join(outdir, 'sgp4test.dbg.' + catno + '.out'), 'wt')
        dbgfile.write('this is the debug output\n\n')
    # convert the char string to sgp4 elements
    # includes initialization of sgp4
    startmfe, stopmfe, deltamin, satrec = sc.twoline2rv(longstr1, longstr2, typerun, typeinput, opsmode, whichconst)
    outfile.write('%d xx\n' % (satrec['satnum']))
    print(' %d\n' % (satrec['satnum']))
    # call the propagator to get the initial state vector value
    satrec, rteme, vteme = obu.sgp4(satrec, 0.0)
    outfile.write(' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n'
                  % (satrec['t'], rteme[0], rteme[1], rteme[2], vteme[0],
                     vteme[1], vteme[2]))
    # fprintf(outfile, ' #16.7f #16.7f #16.7f #16.7f #12.8f #12.8f #12.8f\n', ...
    #     satrec['t'], ro(1), ro(2), ro(3), vo(1), vo(2), vo(3));
    # fprintf(1, ' #16.8f #16.8f #16.8f #16.8f #12.9f #12.9f #12.9f\n', ...
    #     satrec['t'], ro(1), ro(2), ro(3), vo(1), vo(2), vo(3));
    tsince = startmfe
    # // check so the first value isn't written twice
    if (abs(tsince) > smalle8):
        tsince = tsince - deltamin
    eqeterms = 2
    # // loop to perform the propagation
    while ((tsince < stopmfe) and (satrec['error'] == 0)):
        tsince = tsince + deltamin
        if (tsince > stopmfe):
            tsince = stopmfe
        satrec, rteme, vteme = obu.sgp4(satrec, tsince)
        if (satrec['error'] > 0):
            print('# *** error: t:= %f *** code = %3i\n'
                  % (tsince, satrec['error']))
        if (satrec['error'] == 0):
            if ((typerun != 'v') and (typerun != 'c')):
                jdutc = satrec['jdsatepoch']
                jdutcfrac = satrec['jdsatepochf'] + tsince / 1440.0
                if jdutcfrac < 0.0:
                    jdutc = jdutc - 1.0
                    jdutcfrac = jdutcfrac + 1.0
                year, mon, day, hr, minute, sec = stu.invjday(jdutc, jdutcfrac)
                outfile.write(' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f '
                              '%12.9f %5i%3i%3i %2i:%2i:%9.6f  \n'
                              % (tsince, rteme[0], rteme[1], rteme[2],
                                 vteme[0], vteme[1], vteme[2],
                                 year, mon, day, hr, minute, sec))
                # #16.8f#16.8f#16.8#12.9f#12.9f#12.9f and \r
                # demo getting lat lon
                # set the EOP parameters. normally this would be at each
                # time instant in the loop, but
                # we'll keep it constant here for illustration
                lod = 0.0015563
                xp = - 0.140682 * arcsec2rad
                yp = 0.333309 * arcsec2rad
                dut1 = - 0.4399619
                jdut1 = jdutc
                jdut1frac = jdutcfrac + dut1 / 86400.0
                ttt = (jdut1 - 2451545.0) / 36525.0
                recef, vecef, aecef = sc.teme2ecef(rteme.T, vteme.T, ateme.T,
                                                   ttt, jdut1, lod, xp, yp,
                                                   eqeterms)
                #[latgc, latgd, lon, hellp] = ecef2ll ( recef );
                #fprintf(1, ' lat lon #11.7f  #11.7f \n', latgc * rad2deg, lon * rad2deg);
            else:
                jdutc = satrec['jdsatepoch']
                jdutcfrac = satrec['jdsatepochf'] + tsince / 1440.0
                if jdutcfrac < 0.0:
                    jdutc = jdutc - 1.0
                    jdutcfrac = jdutcfrac + 1.0
                year, mon, day, hr, minute, sec = stu.invjday(jdutc, jdutcfrac)
                outfile.write(' %16.8f %16.8f %16.8f %16.8f %12.8f %12.8f '
                              '%12.8f'
                              % (tsince, rteme[0], rteme[1], rteme[2],
                                 vteme[0], vteme[1], vteme[2]))
                # fprintf(1, ' #16.8f #16.8f #16.8f #16.8f #12.9f #12.9f #12.9f \n', ...
                # tsince, ro(1), ro(2), ro(3), vo(1), vo(2), vo(3));
                p, a, ecc, incl, node, argp, nu, m, arglat, truelon, lonper = \
                    sc.rv2coe(rteme, vteme)
                outfile.write(' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f '
                              '%10.5f %5i%3i%3i %2i:%2i:%9.6f \n'
                              % (a, ecc, incl * rad2deg, node * rad2deg,
                                 argp * rad2deg, nu * rad2deg, m * rad2deg,
                                 year, mon, day, hr, minute, sec))
                # demo getting lat lon.
                # set some EOP parameters. normally this would be at each
                # time instant in the loop, but
                # we'll keep it constant here for illustration
                lod = 0.0015563
                xp = - 0.140682 * arcsec2rad
                yp = 0.333309 * arcsec2rad
                dut1 = - 0.4399619
                jdut1 = jdutc
                jdut1frac = jdutcfrac + dut1 / 86400.0
                ttt = (jdut1 - 2451545.0) / 36525.0
                recef, vecef, aecef = sc.teme2ecef(rteme.T, vteme.T, ateme.T,
                                                   ttt, jdut1, lod, xp, yp,
                                                   eqeterms)
                #[latgc, latgd, lon, hellp] = ecef2ll ( recef );
#fprintf(1, ' lat lon #11.7f  #11.7f \n', latgc * rad2deg, lon * rad2deg);

    if (idebug and dbgfile != -1):
        dbgfile.close()


# fclose(infile)
# fclose(outfile)


# sgp4fix demonstrate method of running SGP4 directly from orbital element values
#1 08195U 75081A   06176.33215444  .00000099  00000-0  11873-3 0   813
#2 08195  64.1586 279.0717 6877146 264.7651  20.2257  2.00491383225656

xpdotp = 1440.0 / (2.0 * np.pi)

whichconst = 72
opsmode = 'a'
satrec['satnum'] = 8195
satrec['jdsatepoch'] = 2453911.0
satrec['jdsatepochf'] = 0.8321544402
satrec['no_kozai'] = 2.00491383
satrec['ecco'] = 0.6877146
satrec['inclo'] = 64.1586
satrec['nodeo'] = 279.0717
satrec['argpo'] = 264.7651
satrec['mo'] = 20.2257
satrec['nddot'] = 0.0
satrec['bstar'] = 0.00011873
satrec['ndot'] = 9.9e-07
satrec['elnum'] = 813
satrec['revnum'] = 22565
satrec['classification'] = 'U'
satrec['intldesg'] = '75081A'
satrec['ephtype'] = 0

# convert units and initialize
satrec['no_kozai'] = satrec['no_kozai'] / xpdotp
satrec['ndot'] = satrec['ndot'] / (xpdotp * 1440.0)
satrec['nddot'] = satrec['nddot'] / (xpdotp * 1440.0 * 1440)
satrec['inclo'] = satrec['inclo'] * deg2rad
satrec['nodeo'] = satrec['nodeo'] * deg2rad
satrec['argpo'] = satrec['argpo'] * deg2rad
satrec['mo'] = satrec['mo'] * deg2rad

# set start/stop times for propagation
startmfe = 0.0
stopmfe = 2880.0
deltamin = 120.0
satrec = obu.sgp4init(whichconst, opsmode, satrec, satrec['jdsatepoch']
                      + satrec['jdsatepochf']
                      - 2433281.5, satrec['bstar'],
                      satrec['ndot'], satrec['nddot'],
                      satrec['ecco'], satrec['argpo'],
                      satrec['inclo'], satrec['mo'],
                      satrec['no_kozai'], satrec['nodeo'])
tsince = startmfe
while ((tsince < stopmfe) and (satrec['error'] == 0)):
    tsince = tsince + deltamin
    if (tsince > stopmfe):
        tsince = stopmfe
    satrec, rteme, vteme = obu.sgp4(satrec, tsince)
    jdutc = satrec['jdsatepoch']
    jdutcfrac = satrec['jdsatepochf'] + tsince / 1440.0
    if jdutcfrac < 0.0:
        jdutc = jdutc - 1.0
        jdutcfrac = jdutcfrac + 1.0
    year, mon, day, hr, minute, sec = stu.invjday(jdutc, jdutcfrac)
    print(' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f \n'
          % (tsince, rteme[0], rteme[1], rteme[2],
             vteme[0], vteme[1], vteme[2]))

