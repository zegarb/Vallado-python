import numpy as np
import os
from space_constants import *
import spacemath_utils as smu

st = 'y'
outpath = os.path.join(os.path.dirname(__file__), 'testoutput', 'fig712.out')
outfile = open(outpath)
for kt in range(1, 8):
    # ---- use for fixed r and r locations, and fig 7-12, 16, and 20 ---- }
    #      readln( infile, ro(1), ro(2), ro(3), rtgt(1), rtgt(2), rtgt(3), dt, dmy, direc );
    #      fprintf( kt:3,'fig 7-12 ro rtgt ', ro(4):11:7, rtgt(4):11:7 );
    kepmov = 'n'
    if kt == 1:
        rinto = np.array([1.04, 0.0, 0.0]) * re
        rtgto = np.array([-1.04, 0.2, 0.0]) * re
        dt = 20 * tumin
        direc = 's'
    if kt == 2:
        rinto = np.array([2.0, 0.0, 0.0]) * re
        rtgto = np.array([-2.0, 0.2, 0.0]) * re
        dt = 20 * tumin
        direc = 's'
    if kt == 3:
        rinto = np.array([4.0, 0.0, 0.0]) * re
        rtgto = np.array([-4.0, 0.2, 0.0]) * re
        dt = 20 * tumin
        direc = 's'
    if kt == 4:
        rinto = np.array([6.6107, 0.0, 0.0]) * re
        rtgto = np.array([-6.6107, 0.2, 0.0]) * re
        dt = 20 * tumin
        direc = 's'
    # make up a velocity vector (circular orbit) for dv calcs...
    if (kt > 0) and (kt < 5):
        vinto = np.array([0.0, math.sqrt(mu / smu.mag(rinto)), 0.0])
        vtgto = np.array([0.0, math.sqrt(mu / smu.mag(rtgto)), 0.0])
    if kt == 5:
        rinto = np.array([-6518.1083,-2403.8479,-22.1722])
        vinto = np.array([2.604057,-7.105717,-0.263218])
        rtgto = np.array([6697.4756, 1794.5831, 0.0])
        vtgto = np.array([-1.962372, 7.323674, 0.0])
        dt = 20 * tumin
        direc = 's'
        kepmov = 'n'
    if kt == 6:
        rinto = np.array([5328.7862, 4436.1273, 101.472])
        vinto = np.array([-4.864779, 5.816486, 0.240163])
        rtgto = np.array([6697.4756, 1794.5831, 0.0])
        vtgto = np.array([-1.962372, 7.323674, 0.0])
        dt = 20 * tumin
        direc = 's'
        kepmov = 'y'
    if kt == 7:
        trangle = 90.0 * deg2rad
        rinto = np.array([1.0, 0.0, 0.0])
        vinto = np.array([2.604057,-7.105717,- 0.263218])
        rtgto = np.array([1.0 * math.cos(trangle), 1.0 * math.sin(trangle), 0.0])
        vtgto = np.array([-1.962372, 7.323674, 0.0])
        dt = 20 * tumin
        direc = 's'
        kepmov = 'n'
        #    mu = 4.0*pi*pi;
    direc = 's'
    nrev = 0
    dolam
    direc = 'l'
    nrev = 0
    dolam
    direc = 's'
    nrev = 1
    dolam
    direc = 'l'
    nrev = 1
    dolam
    direc = 's'
    nrev = 2
    dolam
    direc = 'l'
    nrev = 2
    dolam


print('done with fig 7-13 data ')