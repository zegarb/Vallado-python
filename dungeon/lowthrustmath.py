

if __name__ == '__main__':
    import math
    x = 0.5
    a0 = 1.38629436112
    a1 = 0.09666344259
    a2 = 0.03590092383
    a3 = 0.03742563713
    a4 = 0.01451196212

    b0 = 0.5
    b1 = 0.12498593597
    b2 = 0.06880248576
    b3 = 0.03328355346
    b4 = 0.00441787012


    K1 = (a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)
    K2 = (b0 + b1*x + b2*x**2 + b3*x**3 + b4*x**4)*math.log(1/x)
    K = K1 + K2
    P = math.sqrt(1-x)*K
    Pdot1 = (-a0 + 2*a1 - 3*a1*x + 4*a2*x - 5*a2*x**2 + 6*a3*x**2 - 7*a3*x**2 + 8*a4*x**3 - 9*a4*x**4) / (2*math.sqrt(1-x))
    Pdot2 = -(K1*math.sqrt(1-x))/x - math.log(x) * Pdot1
    Pdot = Pdot1 + Pdot2

    c1 = 0.44325141463
    c2 = 0.06260601220
    c3 = 0.04757383546
    c4 = 0.01736506451

    d1 = 0.24998368310
    d2 = 0.09200180037
    d3 = 0.04069697526
    d4 = 0.00526449639

    E1 = (1 + c1*x + c2*x**2 + c3*x**3 + c4*x**4)
    E2 = (d1*x + d2*x**2 + d3*x**3 + d4*x**4)*math.log(1/x)
    E = E1 + E2
    R1 = (1/x)*E
    R2 = ((x-1)/math.sqrt(x))*K
    R = R1 + R2

    Rdot1a = (-1 + c2*x**2 + 2*c3*x**3 + 3*c4*x**4)/(x**2)
    Rdot1b = -((1/x)*(d1 + d2*x + d3*x**2 + d4*x**3) + (d2 + 2*d3*x + 3*d4*x**2)*math.log(x))
    Rdot1 = Rdot1a + Rdot1b

    Rdot2a = (a0 + a0*x - a1*x + 3*a1*x**2 - 3*a2*x**2 + 5*a2*x**3 - 5*a3*x**3 + 7*a3*x**4 - 7*a4*x**4 + 9*a4*x**5) / (2*x*math.sqrt(x))
    #Deriv Calc
    Rdot2b = - ((-2*b0 + (2*b0-2*b1)*x + (2*b1 - 2*b2)*x**2 + (2*b2 - 2*b3)*x**3 + (2*b3 - 2*b4)*x**4 + 2*b4*x**5) + (b0 + (b0 - b1)*x + (3*b1 - 3*b2)*x**2 + (5*b2 - 5*b3)*x**3 + (7*b3 - 7*b4)*x**4 + 9*b4*x**5)*math.log(x)) / (2*x**(3/2))
    #Symbolab
    #Rdot2bsym = -((((x-1)*(b0 + b1*x + b2*x**2 + b3*x**3 + b4*x**4))/(x) + math.log(x) * (5*b4*x**4 + 4*b3*x**3 - 4*b4*x**3 + 3*b2*x**2 - 3*b3*x**2 + 2*b1*x - 2*b2*x + b0 - b1)) * (math.sqrt(x)/x) - ((b0 + b1*x + b2*x**2 + b3*x**3 + b4*x**4)*math.log(x)*(x-1))/(2*x**(3/2)))

    #print(Rdot2b)
    #print(Rdot2bsym)
    Rdot2 = Rdot2a + Rdot2b

