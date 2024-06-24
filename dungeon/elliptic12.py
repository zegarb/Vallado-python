import numpy as np
import math
import sys

# ELLIPTIC12 evaluates the value of the Incomplete Elliptic Integrals
# of the First, Second Kind and Jacobi's Zeta Function.

#   [F, E, Z] = ELLIPTIC12(U, M, TOL) where U is a phase in radians, 0<M<1 is
#   the module and TOL is the tolerance (optional). Default value for
#   the tolerance is eps = 2.220e-16.

#   ELLIPTIC12 uses the method of the Arithmetic-Geometric Mean
#   and Descending Landen Transformation described in [1] Ch. 17.6,
#   to determine the value of the Incomplete Elliptic Integrals
#   of the First, Second Kind and Jacobi's Zeta Function [1], [2].

#       F(phi, m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);
#       E(phi, m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
#       Z(phi, m) = E(u, m) - E(m)/K(m)*F(phi, m).

#   Tables generating code ([1], pp. 613-621):
#       [phi, alpha] = meshgrid(0:5:90, 0:2:90);                  # modulus and phase in degrees
#       [F, E, Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  # values of integrals

#   See also ELLIPKE, ELLIPJ, ELLIPTIC12I, ELLIPTIC3, THETA, AGM.

#   References:
#   [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical Functions",
#       Dover Publications", 1965, Ch. 17.1 - 17.6 (by L.M. Milne-Thomson).
#   [2] D. F. Lawden, "Elliptic Functions and Applications"
#       Springer-Verlag, vol. 80, 1989

# Copyright Elliptic Project 2011
# For support, please reply to
#     moiseev.igor[at]gmail.com
#     Moiseev Igor,

# The code is optimized for ordered inputs produced by the functions
# meshgrid, ndgrid. To obtain maximum performace (up to 30#) for singleton,
# 1-dimensional and random arrays remark call of the function unique(.)
# and edit further code.

def elliptic12(u, m, tol = None):
    if tol == None:
        tol = sys.float_info.epsilon

    if not np.isreal(u).all() or not np.isreal(m).all():
        raise Exception('Input arguments must be real. Use ELLIPTIC12i for complex arguments.')

    if len(m) == 1:
        m = m(np.ones(u.shape))

    if len(u) == 1:
        u = u(np.ones(m.shape))

    if not m.shape == u.shape:
        raise Exception('U and M must be the same size.')

    F = np.zeros(u.shape)
    E = F.copy()
    Z = E.copy()

    m = m.T
    u = u.T

    if np.any(m < 0) or np.any(m > 1):
        raise Exception('M must be in the range 0 <= M <= 1.')

    # cdav change for small eccentricities
    # for i in range(len(m)):
    #     if abs(m[i]) < 1e-07:
    #         m[i] = 1e-07

    I = np.nonzero(np.logical_and(m != 1, m != 0))[0]
    if len(I):
        # mu, J, K = unique(m(I))
        mu, K = np.unique(m[I], False, True)
        # K = uint32(K)
        mumax = len(mu)
        signU = np.sign(u[I])
        # pre-allocate space and augment if needed
        chunk = 7
        a = np.zeros((chunk, mumax))
        c = a.copy()
        b = a.copy()
        a[0, :] = np.ones(mumax)
        c[0, :] = np.sqrt(mu)
        b[0, :] = np.sqrt(1 - mu)
        n = np.zeros(mumax)
        i = 0
        while (abs(c[i, :]) > tol).any(): # Arithmetic-Geometric Mean of A, B and C

            i = i + 1
            if i > a.shape[0]:
                np.append(a, np.zeros((1, mumax)), axis=0)
                np.append(b, np.zeros((1, mumax)), axis=0)
                np.append(c, np.zeros((1, mumax)), axis=0)

            a[i, :] = 0.5 * (a[i - 1, :] + b[i - 1, :])
            b[i, :] = np.sqrt(a[i - 1, :] * b[i - 1, :])
            c[i, :] = 0.5 * (a[i - 1, :] - b[i - 1, :])
            in_ = np.nonzero(np.logical_and((abs(c[i, :]) <= tol),
                                            (abs(c[i - 1, :]) > tol)))[0]
            if len(in_):
                n[in_] = i - 1

        mmax = len(I)
        mn = int(np.max(n))
        phin = np.zeros((1, mmax))
        C = np.zeros(mmax)
        Cp = C.copy()
        e = C.copy()
        phin = signU * u[I]
        i = 0
        c2 = c ** 2
        while i < mn: # Descending Landen Transformation

            in_ = np.nonzero(n[K] > i)[0]
            if len(in_):
                phin[in_] = (np.arctan(b[i, K[in_]] / a[i, K[in_]]
                                         * np.tan(phin[in_]))
                               + np.pi * np.ceil(phin[in_] / np.pi - 0.5)
                               + phin[in_])
                e[in_] = 2.0 ** (i - 1)
                C[in_] = C[in_] + e[in_[0]] * c2[i, K[in_]]
                Cp[in_] = Cp[in_] + c[i + 1, K[in_]] * np.sin(phin[in_])
            i = i + 1

        Ff = phin / (a[mn, K] * e * 2)
        F[I] = Ff * signU
        Z[I] = Cp * signU
        E[I] = (Cp + ((1 - 1 / 2 * C) *  Ff)) * signU

    # Special cases: m == {0, 1}
    m0 = np.nonzero(m == 0)[0]
    if len(m0):
        F[m0] = u[m0]
        E[m0] = u[m0]
        Z[m0] = 0

    m1 = np.nonzero(m == 1)[0]
    um1 = abs(u[m1])
    if len(m1):
        N = np.floor((um1 + np.pi / 2) / np.pi)
        M = np.nonzero(um1 < np.pi / 2)[0]
        F[m1[M]] = np.log(np.tan(np.pi / 4 + u[m1[M]] / 2))
        F[m1[um1 >= np.pi / 2]] = math.inf *  np.sign(u[m1[um1 >= np.pi / 2]])
        E[m1] = (((-1) ** N * np.sin(um1)) + 2 * N) * np.sign(u[m1])
        Z[m1] = -1 ** N * np.sin(u[m1])

    return F, E, Z


m = np.array([0, .1, .25, .4, .5, .4, .25, .5, 1, 0, 1])
u = np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.9, 20])
F, E, Z = elliptic12(u, m)
print(f'{F}\n\n{E}\n\n{Z}')