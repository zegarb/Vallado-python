import numpy as np
import math
from elliptic12 import elliptic12
#ARCLENGTH_ELLIPSE Calculates the arclength of ellipse.
#
#   ARCLENGTH_ELLIPSE(A, B, THETA0, THETA1) Calculates the arclength of ellipse
#   using the precise formulas based on the representation of
#   the arclength by the Elliptic integral of the second kind.
#
#   Ellipse parameters:
#       T - measured in radians from 0 in the positive direction,
#           Period: 2*Pi
#       A - major axis
#       B - minor axis
#
#   Parametric equations:
#       x(t) = a.cos(t)
#       y(t) = b.sin(t)
#
#   Cartesian equation:
#   x^2/a^2 + y^2/b^2 = 1
#
#   Eccentricity:
#       e = Sqrt(1 - (a/b)^2)
#
#   Focal parameter:
#       b^2/Sqrt(a^2 - b^2)
#
#   Foci:
#       (-Sqrt(a^2 - b^2), 0)   OR   (Sqrt(a^2 - b^2), 0)
#
#   Arclength:
#       b*EllipticE( t, 1 - (a/b)^2 )
#
#   Mathematica Test 1:
#       In:= b = 10; a = 5;
#            SetPrecision[b*EllipticE[2Pi, 1.0- a^2/b^2],20]
#      Out:= 48.442241102738385905
#
#   Mathematica Test 2:
#       In:= b = 10; a = 5;
#            SetPrecision[b*(EllipticE[Pi/2-Pi/10, 1.0- a^2/b^2]-EllipticE[Pi/10, 1.0- a^2/b^2]),20]
#      Out:= 7.3635807913930495516
#
#   MATLAB Test 1:
#       # full ellipse
#       arclength = arclength_ellipse(5,10)
#       arclength =
#           48.442241102738436
#
#   MATLAB Test 2:
#       # arclength ellipse
#       arclength = arclength_ellipse(5,10,pi/10,pi/2)
#       arclength =
#           7.363580791393055
#
#   References:
#   @see http://mathworld.wolfram.com/Ellipse.html
#   @see http://www.wolframalpha.com/input/?i=ellipse+arc+length&lk=1&a=ClashPrefs_*PlaneCurve.Ellipse.PlaneCurveProperty.ArcLength-
#

# Special thanks to for bug correction
#    drbitboy (Brian Carcich) https://github.com/drbitboy
# 2015-07-14 (New Horizons flyby of Pluto)
#
# 1) Old code returned values that were in error
# 1.1)  arclength_ellipse(1., .5, pi*.001, pi*.002) returned 0
# 1.2)  arclength_ellipse(1., .5, pi*.002, pi*.001) returned -.0003*pi instead of pi correct .0005*pi
# 1.3)  arclength_ellipse(1., .5, theta0, theta1) did not return the negative of the same call with the thetas reversed
# 2) Angles theta0 and theta1 were always interpreted as measured from the semi-minor axis
#
# 3) Corrected code:
# 3.1) Angle theta is measured from the positive a axis
# 3.2) The standard form of the b*E(phi,m) arc length integral has m = 1 - (a/b)^2
# 3.2.1) N.B. That only only works if b>a
# 3.3) If a>b, then an alternate formula is used:  a*E(PI/2 - phi, m') where m' = 1 - (b/a)^2
# 3.4) A few simple cases will show that the new code is correct
#        arclength_ellipse(1, .5, pi*.001, pi*.002) ~  pi*.0005
#        arclength_ellipse(1, .5, pi*.002, pi*.001) = -arclength(1, .5, pi*.001, pi*.002) ~ -pi*.0005
#        arclength_ellipse(1., 2., pi*.001, pi*.002) ~ pi*.002
#        arclength_ellipse(1, .5, pi/2 - pi*.002, pi/2 - pi*.001) ~ -pi*.001
#        arclength_ellipse(1, 2., pi/2 - pi*.002, pi/2 - pi*.001) ~ -pi*.001
#        etc.

# Copyright Elliptic Project 2011
# For support,
#     moiseev.igor[at]gmail.com
#     Moiseev Igor

def arclength_ellipse(a: np.ndarray, b: np.ndarray, theta0: np.ndarray = None,
                      theta1: np.ndarray = None):

    if not isinstance(a, np.ndarray):
        a = np.array([a])

    if not isinstance(b, np.ndarray):
        b = np.array([b])

    if not isinstance(theta0, np.ndarray):
        theta0 = np.array([theta0])

    if not isinstance(theta1, np.ndarray):
        theta1 = np.array([theta1])

    arclength = a * (theta1 - theta0)

    if not theta0.any() or not theta1.any():
        if not theta0.any() and not theta1.any():
            theta0 = np.full(a.shape, 0)
            theta1 = np.full(a.shape, 2*math.pi)
        else:
            print("Error: requires both theta0 and theta1 set or neither!")
            return

    theta0a = np.zeros(len(a))
    theta1a = np.zeros(len(a))
    ab = np.zeros(len(a))
    for i in range(len(a)):
        if a[i] < b[i]:
            ab[i] = 1 - (a[i] / b[i]) ** 2
            theta0a[i] = theta0[i]
            theta1a[i] = theta1[i]
        elif a[i] > b[i]:
            ab[i] =  1 - (b[i] / a[i]) ** 2
            theta0a[i] = math.pi / 2 - theta0[i]
            theta1a[i] = math.pi / 2 - theta1[i]

    F1, E1, Z1 = elliptic12(theta1a, ab)
    F0, E0, Z0 = elliptic12(theta0a, ab)

    arclength = np.zeros(len(a))
    for i in range(len(a)):
        if a[i] < b[i]:
            arclength[i] = b[i] * (E1[i] - E0[i])
        elif a[i] > b[i]:
            arclength[i] = a[i] * (E0[i] - E1[i])

    return arclength

if __name__ == '__main__':
    # test input arrays. 0-4 is tests from lines 78-82, 5 is "MATLAB test 1",
    # and 6 is "MATLAB test 2." All outputs match the answers except MATLAB test 2.
    answers = np.array([.0005, -.0005, .002, .001, .001, 48.4422411 / math.pi,
                        7.36358/math.pi])
    a = np.array([1, 1, 1, 1, 1, 5, 5])
    b = np.array([.5, .5, 2, .5, 2, 10, 10])
    theta0 = np.array([.001, .002, .001, .5 - .002, .5 - .002, 0, 0.1])
    theta1 = np.array([.002, .001, .002, .5 - .001, .5 - .001, 2, 0.5])
    theta0 = theta0 * math.pi
    theta1 = theta1 * math.pi

    arclength = arclength_ellipse(a, b, theta0, theta1)

    print(f'{answers}')
    print(f'{arclength / math.pi}')
    print(f'{(arclength / math.pi) / answers}')
