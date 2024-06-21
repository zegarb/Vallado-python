

    # ELLIPTIC12 evaluates the value of the Incomplete Elliptic Integrals
# of the First, Second Kind and Jacobi's Zeta Function.

    #   [F,E,Z] = ELLIPTIC12(U,M,TOL) where U is a phase in radians, 0<M<1 is
#   the module and TOL is the tolerance (optional). Default value for
#   the tolerance is eps = 2.220e-16.

    #   ELLIPTIC12 uses the method of the Arithmetic-Geometric Mean
#   and Descending Landen Transformation described in [1] Ch. 17.6,
#   to determine the value of the Incomplete Elliptic Integrals
#   of the First, Second Kind and Jacobi's Zeta Function [1], [2].

    #       F(phi,m) = int(1/sqrt(1-m*sin(t)^2), t=0..phi);
#       E(phi,m) = int(sqrt(1-m*sin(t)^2), t=0..phi);
#       Z(phi,m) = E(u,m) - E(m)/K(m)*F(phi,m).

    #   Tables generating code ([1], pp. 613-621):
#       [phi,alpha] = meshgrid(0:5:90, 0:2:90);                  # modulus and phase in degrees
#       [F,E,Z] = elliptic12(pi/180*phi, sin(pi/180*alpha).^2);  # values of integrals

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

def elliptic12(u = None,m = None,tol = None):
    if len(varargin) < 3:
        tol = eps

    if len(varargin) < 2:
        raise Exception('Not enough input arguments.')

    if not True  or not True :
        raise Exception('Input arguments must be real. Use ELLIPTIC12i for complex arguments.')

    if len(m) == 1:
        m = m(np.ones((u.shape,u.shape)))

    if len(u) == 1:
        u = u(np.ones((m.shape,m.shape)))

    if not m.shape==u.shape :
        raise Exception('U and M must be the same size.')

    F = np.zeros((u.shape,u.shape))
    E = F
    Z = E
    m = np.transpose(m)

    u = np.transpose(u)
    if np.any(m < 0) or np.any(m > 1):
        raise Exception('M must be in the range 0 <= M <= 1.')

    # cdav change for small eccentricities
    if np.abs(m) < 1e-07:
        m = 1e-07

    I = uint32(find(m != np.logical_and(1,m) != 0))
    if not len(I)==0 :
        mu,J,K = unique(m(I))
        K = uint32(K)
        mumax = len(mu)
        signU = np.sign(u(I))
        # pre-allocate space and augment if needed
        chunk = 7
        a = np.zeros((chunk,mumax))
        c = a
        b = a
        a[1,:] = np.ones((1,mumax))
        c[1,:] = np.sqrt(mu)
        b[1,:] = np.sqrt(1 - mu)
        n = uint32(np.zeros((1,mumax)))
        i = 1
        while np.any(np.abs(c(i,:)) > tol):

            i = i + 1
            if i > a.shape[1-1]:
                a = np.array([[a],[np.zeros((2,mumax))]])
                b = np.array([[b],[np.zeros((2,mumax))]])
                c = np.array([[c],[np.zeros((2,mumax))]])
            a[i,:] = 0.5 * (a(i - 1,:) + b(i - 1,:))
            b[i,:] = np.sqrt(np.multiply(a(i - 1,:),b(i - 1,:)))
            c[i,:] = 0.5 * (a(i - 1,:) - b(i - 1,:))
            in_ = uint32(find(np.logical_and((np.abs(c(i,:)) <= tol),(np.abs(c(i - 1,:)) > tol))))
            if not len(in_)==0 :
                mi,ni = in_.shape
                n[in_] = np.ones((mi,ni)) * (i - 1)

        mmax = len(I)
        mn = double(np.amax(n))
        phin = np.zeros((1,mmax))
        C = np.zeros((1,mmax))
        Cp = C
        e = uint32(C)
        phin = np.multiply(signU,u(I))
        i = 0
        c2 = c ** 2
        while i < mn:

            i = i + 1
            in_ = uint32(find(n(K) > i))
            if not len(in_)==0 :
                phin[in_] = np.arctan(np.multiply(b(i,K(in_)) / a(i,K(in_)),np.tan(phin(in_)))) + np.multiply(np.pi,np.ceil(phin(in_) / np.pi - 0.5)) + phin(in_)
                e[in_] = 2.0 ** (i - 1)
                C[in_] = C(in_) + double(e(in_(1))) * c2(i,K(in_))
                Cp[in_] = Cp(in_) + np.multiply(c(i + 1,K(in_)),np.sin(phin(in_)))

        Ff = phin / (np.multiply(a(mn,K),double(e)) * 2)
        F[I] = np.multiply(Ff,signU)
        Z[I] = np.multiply(Cp,signU)
        E[I] = np.multiply((Cp + np.multiply((1 - 1 / 2 * C),Ff)),signU)

    # Special cases: m == {0, 1}
    m0 = find(m == 0)
    if not len(m0)==0 :
        F[m0] = u(m0)
        E[m0] = u(m0)
        Z[m0] = 0

    m1 = find(m == 1)
    um1 = np.abs(u(m1))
    if not len(m1)==0 :
        N = int(np.floor((um1 + np.pi / 2) / np.pi))
        M = find(um1 < np.pi / 2)
        F[m1[M]] = np.log(np.tan(np.pi / 4 + u(m1(M)) / 2))
        F[m1[um1 >= np.pi / 2]] = np.multiply(Inf,np.sign(u(m1(um1 >= np.pi / 2))))
        E[m1] = np.multiply((np.multiply((- 1) ** N,np.sin(um1)) + 2 * N),np.sign(u(m1)))
        Z[m1] = np.multiply((- 1) ** N,np.sin(u(m1)))

    return F,E,Z
