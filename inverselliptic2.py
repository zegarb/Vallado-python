import sys
import numpy as np
from elliptic12 import elliptic12
# INVERSELLIPTIC2 evaluates the value of the INVERSE Incomplete Elliptic Integrals
# of the Second Kind.

# INVE = INVERSELLIPTIC2(E,M,TOL) where E is a value of the integral to
# be inverted, 0<M<1 is the module and TOL is the tolerance (optional).
# Default value for the tolerance is eps = 2.220e-16.

# INVERSELLIPTIC2 uses the method described by Boyd J. P.
# to determine the value of the inverse Incomplete Elliptic Integrals
# of the Second Kind using the Empirical initialization to
# the Newtons iteration method [1].

# NOTICE. Please pay attention to the definition of the elliptic functions
# which follows the Abramovitz et al [2], for more theory on elliptic
# functions please consult the Lawden book [3].

# Elliptic integral of the second kind:

# E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);

# Empirical initialization [1]:

# T0(z,m) = pi/2 + sqrt(r)/(theta ? pi/2)

# where
# z \in [?E(pi/2,m), E(pi/2,m)]x[0, 1], value of the entire parameter space
# r = sqrt((1-m)^2 + zeta^2)
# zeta = 1 - z/E(pi/2,m)
# theta = atan((1 - m)/zeta)


# Example:
# # modulus and phase in degrees
# [phi,alpha] = meshgrid(0:5:90, 0:2:90);
# # values of integrals
# [F,E] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);
# # values of inverse
# invE = inverselliptic2(E, sin(pi/180*alpha).^2);
# # the difference between phase phi and invE should close to zero
# phi - invE * 180/pi

# See also ELLIPKE, ELLIPTIC12.

# References:
# [1] J. P. Boyd, "Numerical, Perturbative and Chebyshev Inversion
# of the Incomplete Elliptic Integral of the Second Kind", Applied Mathematics and Computation (January 2012)
# [2] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
# Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
# [3] D. F. Lawden, "Elliptic Functions and Applications"
# Springer-Verlag, vol. 80, 1989

# Copyright (C) 2011 by Elliptic Project. All rights reserved.

# GNU GENERAL PUBLIC LICENSE Version 2, June 1991
# http://www.gnu.org/licenses/gpl.html
# Everyone is permitted to copy and distribute verbatim copies of this
# script under terms and conditions of GNU GENERAL PUBLIC LICENSE.

# For support, please reply to
# moiseev.igor[at]gmail.com
# Moiseev Igor

# ELLIPTIC PROJECT: http://elliptic.googlecode.com
# Group:

def inverselliptic2(E: np.ndarray, m: np.ndarray, tol: float = None):
    if not isinstance(E, np.ndarray):
        E = np.array([E])

    if not isinstance(m, np.ndarray):
        m = np.array([m])

    if tol == None:
        tol = sys.float_info.epsilon

    if not np.isreal(E).all() or not np.isreal(m).all() :
        raise Exception('Input arguments must be real.')

    # if len(m) == 1:
    #     m = m(np.ones(E.shape))

    # if len(E) == 1:
    #     E = E(np.ones(m.shape))

    # if not m.shape == E.shape:
    #     raise Exception('E and M must be the same size.')

    invE = np.zeros(E.shape)
    # make a row vector
    m = m
    E = E

    if np.any(m < 0) or np.any(m > 1):
        raise Exception('M must be in the range 0 <= M <= 1.')

    # cdav change for small eccentricities
    # for i in range(len(m)):
    #     if abs(m[i]) < 1e-07:
    #         m[i] = 1e-07


    # inputs
    z = E
    mu = 1 - m
    # complete integral initialization
    # E1 = special.ellipe(m, tol)
    u = np.full(m.shape, np.pi / 2)
    __, E1, _ = elliptic12(u, m)
    zeta = 1 - z / E1
    r = np.sqrt(zeta * zeta + mu * mu)
    theta = np.arctan(mu / (z + sys.float_info.epsilon))
    # Empirical initialization [1]
    invE = np.pi / 2 + np.sqrt(r) * (theta - (np.pi / 2))
    for _ in range(4):
        __, E, _ = elliptic12(invE, m, tol)
        invE = invE -(E - z) / np.sqrt(1 - m * np.sin(invE) ** 2)

    return invE
