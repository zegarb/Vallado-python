import math
import numpy as np
from space_conversions import coe2rv, rv2coe, eq2rv, rv2eq
from spacemath_utils import mag
from space_constants import *


#     -----------------------------------------------------------------
#
#                              Ex2_5.m
#
#  this file demonstrates example 2-5. it also includes some stressing
#  cases for the coe and rv conversions for all orbit types.
#
#                          companion code for
#             fundamentals of astrodynamics and applications
#                                 2013
#                            by david vallado
#
#     (h)               email davallado@gmail.com
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#            16 feb 19  david vallado
#                         update for new constants
#  changes :
#            13 feb 07  david vallado
#                         original baseline
#
#     *****************************************************************



print('coe quick test ----------------------------' )

#r=np.array([ 6524.834], [6862.875], [6448.296]])
#v=np.array([ 4.901327], [5.533756], [-1.976341]])

r=np.array([6524.834, 6862.875, 6448.296])
v=np.array([4.901327, 5.533756, -1.976341])

print('start %15.9f %15.9f %15.9f' % (r[0], r[1], r[2]))
print(' v %15.10f %15.10f %15.10f\n' % (v[0], v[1], v[2]))

# --------  rv2coe   - position and velocity vectors to classical elements

p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(r, v)

ideg = incl * rad2deg if incl else -0.0
odeg = omega * rad2deg if omega else -0.0
adeg = argp * rad2deg if argp else -0.0
ndeg = nu * rad2deg if nu else -0.0
mdeg = m * rad2deg if m else -0.0
argdeg = arglat * rad2deg if arglat else -0.0
tdeg = truelon * rad2deg if truelon else -0.0
ldeg = lonper * rad2deg if lonper else -0.0

print('rv2coe:  p km        a km         ecc            incl deg     raan deg     argp deg    nu deg       m deg     arglat       truelon     lonper')
print('coes    %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n'
       % (p, a, ecc, ideg, odeg, adeg, ndeg, mdeg, argdeg, tdeg, ldeg))

#pause



# alt test various combinations of coe/eq and rv
for j in range(1):
    if j == 0:
        print('coe tests -----------------------------------------------------------\n' )
    else:
        print('\n\neq tests -----------------------------------------------------------\n' )
    #end
    for i in range(20):
        print(f"iteration: {i}")
        if i == 0:
            print('coe test ----------------------------%d\n'%i )
            r=np.array([6524.834, 6862.875, 6448.296])
            v=np.array([4.901327, 5.533756, -1.976341])
        #end
        if i == 1:
            print('coe test ----------------------------%d\n'%i )
            r=np.array([6524.834, 6862.875, 6448.296])
            v=np.array([4.901327, 5.533756, -1.976341])
        #end

        # ------- elliptical orbit tests -------------------
        if i == 2:
            print('coe test elliptical ----------------------------%d\n'%i )
            r=np.array([1.1372844, -1.0534274, -0.8550194])*6378.137
            v=np.array([0.6510489, 0.4521008, 0.0381088])*7.905366149846
        #end
        if i == 3:
            print('coe test elliptical ----------------------------%d\n'%i )
            r=np.array([1.0561942, -0.8950922, -0.0823703])*6378.137
            v=np.array([-0.5981066, -0.6293575, 0.1468194])*7.905366149846
        #end

        # ------- circular inclined orbit tests -------------------
        if i == 4:
            print('coe test near circular inclined ----------------------------%d\n'%i )
            r=np.array([-0.4222777, 1.0078857, 0.7041832])*6378.137
            v=np.array([-0.5002738, -0.5415267, 0.4750788])*7.905366149846
        #end
        if i == 5:
            print('coe test near circular inclined ----------------------------%d\n'%i )
            r=np.array([-0.7309361, -0.6794646, -0.8331183])*6378.137
            v=np.array([-0.6724131, 0.0341802, 0.5620652])*7.905366149846
        #end

        if i == 6: # -- CI u = 45 deg
            print('coe test circular inclined ----------------------------%d\n'%i )
            r=np.array([-2693.34555010128,  6428.43425355863, 4491.37782050409])
            v=np.array([   -3.95484712246016, -4.28096585381370, 3.75567104538731])
        #end
        if i == 7: # -- CI u = 315 deg
            print('coe test circular inclined ----------------------------%d\n'%i )
            r=np.array([-7079.68834483379, 3167.87718823353, -2931.53867301568])
            v=np.array([    1.77608080328182, 6.23770933190509,  2.45134017949138])
        #end

        # ------- elliptical equatorial orbit tests -------------------
        if i == 8:
            print('coe test elliptical near equatorial ----------------------------%d\n'%i )
            r=np.array([ 21648.6109280739, -14058.7723188698, -0.0003598029])
            v=np.array([ 2.16378060719980, 3.32694348486311, 0.00000004164788 ])
        #end
        if i == 9:
            print('coe test elliptical near equatorial ----------------------------%d\n'%i )
            r=np.array([  7546.9914487222,  24685.1032834356, -0.0003598029])
            v=np.array([ 3.79607016047138, -1.15773520476223, 0.00000004164788 ])
        #end

        if i == 10: # -- EE w = 20 deg
            print('coe test elliptical equatorial ----------------------------%d\n'%i )
            r=np.array([-22739.1086596208 , -22739.1086596208  ,   0.0])
            v=np.array([    2.48514004188565,  -2.02004112073465 , 0.0])
        #end
        if i == 11: # -- EE w = 240 deg
            print('coe test elliptical equatorial ----------------------------%d\n'%i )
            r=np.array([ 28242.3662822040  ,  2470.8868808397   , 0.0])
            v=np.array([    0.66575215057746 , -3.62533022188304,  0.0])
        #end

        # ------- circular equatorial orbit tests -------------------
        if i == 12:
            print('coe test circular near equatorial ----------------------------%d\n'%i )
            r=np.array([ -2547.3697454933, 14446.8517254604, 0.000 ])
            v=np.array([  -5.13345156333487, -0.90516601477599, 0.00000090977789 ])
        #end
        if i == 13:
            print('coe test circular near equatorial ----------------------------%d\n'%i )
            r=np.array([  7334.858850000, -12704.3481945462,   0.000 ])
            v=np.array([  -4.51428154312046, -2.60632166411836, 0.00000090977789 ])
        #end

        if i == 14: # -- CE l = 65 deg
            print('coe test circular equatorial ----------------------------%d\n'%i )
            r=np.array([ 6199.6905946008, 13295.2793851394,      0.0])
            v=np.array([ -4.72425923942564, 2.20295826245369,    0.0])
        #end
        if i == 15: # -- CE l = 65 deg i = 180 deg
            print('coe test circular equatorial ----------------------------%d\n'%i )
            r=np.array([ 6199.6905946008, -13295.2793851394,      0.0])
            v=np.array([ -4.72425923942564, -2.20295826245369,    0.0])
        #end

        # ------- parabolic orbit tests -------------------
        if i == 16:
            print('coe test parabolic ----------------------------%d\n'%i )
            r=np.array([  0.5916109, -1.2889359, -0.3738343])*6378.137
            v=np.array([   1.1486347, -0.0808249, -0.1942733])*7.905366149846
        #end

        if i == 17:
            print('coe test parabolic ----------------------------%d\n'%i )
            r=np.array([-1.0343646, -0.4814891,  0.1735524])*6378.137
            v=np.array([ 0.1322278, 0.7785322, 1.0532856  ])*7.905366149846
        #end

        if i == 18:
            print('coe test hyperbolic ---------------------------%d\n'%i )
            r=np.array([0.9163903, 0.7005747, -1.3909623  ])*6378.137
            v=np.array([0.1712704, 1.1036199, -0.3810377  ])*7.905366149846
        #end

        if i == 19:
            print('coe test hyperbolic ---------------------------%d\n'%i )
            r=np.array([12.3160223, -7.0604653, -3.7883759])*6378.137
            v=np.array([-0.5902725, 0.2165037, 0.1628339  ])*7.905366149846
        #end


        print('start %15.9f %15.9f %15.9f' % (r[0], r[1], r[2]))
        print(' v %15.10f %15.10f %15.10f\n' % (v[0], v[1], v[2] ))

        if j == 0:
            # --------  rv2coe       - position and velocity vectors to classical elements
            p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper = rv2coe(r, v)
            ideg = incl * rad2deg if incl else -0.0
            odeg = omega * rad2deg if omega else -0.0
            adeg = argp * rad2deg if argp else -0.0
            ndeg = nu * rad2deg if nu else -0.0
            mdeg = m * rad2deg if m else -0.0
            argdeg = arglat * rad2deg if arglat else -0.0
            tdeg = truelon * rad2deg if truelon else -0.0
            ldeg = lonper * rad2deg if lonper else -0.0
            print('rv2coe:  p km        a km         ecc            incl deg     raan deg     argp deg    nu deg       m deg     arglat       truelon     lonper')
            print('coes    %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f\n'
                % (p, a, ecc, ideg, odeg, adeg, ndeg, mdeg, argdeg, tdeg, ldeg))


            # --------  coe2rv       - classical elements to position and velocity
            rn, vn = coe2rv(p, ecc, incl, omega, argp, nu, arglat, truelon, lonper)
            print('rn    %15.9f %15.9f %15.9f' % (rn[0], rn[1], rn[2]))
            print(' vn %15.10f %15.10f %15.10f\n' % (vn[0], vn[1], vn[2]))
            dr = np.zeros([3])
            dr[0] = r[0] - rn[0]
            dr[1] = r[1] - rn[1]
            dr[2] = r[2] - rn[2]
            if mag(dr) > 0.01:
                print('ERROR in rv/coe case dr = %11.7f \n' % mag(dr))
            #end
        else:
            # --------  rv2eq       - position and velocity vectors to classical elements
            a, n, af, ag, chi, psi, meanlonM, meanlonNu, fr  = rv2eq (r, v)
            print('rv2eq:        a km       n rad     af        ag       chi        psi      meanlonnu deg   meanlonm deg ')
            print('eqs %11.4f %11.4f %13.9f %13.7f %11.5f %11.5f %11.5f %11.5f \n' % \
                (a, n, af, ag, chi, psi, meanlonNu * rad2deg, meanlonM * rad2deg))

            # --------  eq2rv       - classical elements to position and velocity
            [rn, vn] = eq2rv(a, af, ag, chi, psi, meanlonM, fr)
            print('rn    %15.9f %15.9f %15.9f' % (rn[0], rn[1], rn[2]))
            print(' vn %15.10f %15.10f %15.10f\n' % (vn[0], vn[1], vn[2]))
            dr = np.zeros([3])
            dr[0] = r[0] - rn[0]
            dr[1] = r[1] - rn[1]
            dr[2] = r[2] - rn[2]
            if mag(dr) > 0.01:
                print('ERROR in rv/eq case dr = %11.7f \n' % mag(dr))
            #end
        #end

    #end;  % for
#end; % for through coe/eq tests




