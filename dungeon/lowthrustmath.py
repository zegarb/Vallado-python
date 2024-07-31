import math
import numpy as np
import spacemath_utils as smu

def lowthrust():
    #x = 0.446269201
    x = 0.367174867668404
    x1 = 1 - x
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


    K1 = (a0 + a1*x1 + a2*x1**2 + a3*x1**3 + a4*x1**4)
    K2 = (b0 + b1*x1 + b2*x1**2 + b3*x1**3 + b4*x1**4)*math.log(1/x1)
    K = K1 + K2
    K1dot = (a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4)
    K2dot = (b0 + b1*x + b2*x**2 + b3*x**3 + b4*x**4)*math.log(1/x)
    Kdot = K1dot + K2dot

    # K1dot = a1 + 2*a2*x + 3*a3*x**2 + 4*a4*x**3
    # K2dot = (-(1/x)*(b0 + b1*x + b2*x**2 + b3*x**3 + b4*x**4) + math.log(1/x)*(b1 + 2*b2*x + 3*b3*x**2 + 4*b4*x**3))
    # Kdot = K1dot + K2dot
    # Pdot1 = (-a0 + 2*a1 - 3*a1*x1 + 4*a2*x1 - 5*a2*x1**2 + 6*a3*x1**2 - 7*a3*x1**2 + 8*a4*x1**3 - 9*a4*x1**4) / (2*math.sqrt(1-x1))
    # Pdot2 = -((b0 + b1*x1 + b2*x1**2 + b3*x1**3 + b4*x1**4) *math.sqrt(1-x1))/x1 - (math.log(x1)*(-b0 + 2*b1 - 3*b1*x1 + 4*b2*x1 - 5*b2*x1**2 + 6*b3*x1**2 - 7*b3*x1**2 + 8*b4*x1**3 - 9*b4*x1**4))/(2*math.sqrt(1-x1))
    # Pdot = Pdot1 + Pdot2


    c1 = 0.44325141463
    c2 = 0.06260601220
    c3 = 0.04757383546
    c4 = 0.01736506451

    d1 = 0.24998368310
    d2 = 0.09200180037
    d3 = 0.04069697526
    d4 = 0.00526449639

    E1 = (1 + c1*x1 + c2*x1**2 + c3*x1**3 + c4*x1**4)
    E2 = (d1*x1 + d2*x1**2 + d3*x1**3 + d4*x1**4)*math.log(1/x1)
    E = E1 + E2

    E1dot = (1 + c1*x + c2*x**2 + c3*x**3 + c4*x**4)
    E2dot = (d1*x + d2*x**2 + d3*x**3 + d4*x**4)*math.log(1/x)
    Edot = E1dot + E2dot

    Kf1, Ef1, _ = smu.elliptic12(math.pi/2, x)
    Kf2, Ef2, _ = smu.elliptic12(math.pi/2, x1)

    print(Kf1)
    print(K)
    print(Kf2)
    print(Kdot)
    print(Ef1)
    print(E)
    print(Ef2)
    print(Edot)

    # Rdot1a = (-1 + c2*x1**2 + 2*c3*x1**3 + 3*c4*x1**4)/(x1**2)
    # Rdot1b = -((1/x1)*(d1 + d2*x1 + d3*x1**2 + d4*x1**3) + (d2 + 2*d3*x1 + 3*d4*x1**2)*math.log(x1))
    # Rdot1 = Rdot1a + Rdot1b
    # Rdot2a = (a0 + a0*x1 - a1*x1 + 3*a1*x1**2 - 3*a2*x1**2 + 5*a2*x1**3 - 5*a3*x1**3 + 7*a3*x1**4 - 7*a4*x1**4 + 9*a4*x1**5) / (2*x1*math.sqrt(x1))
    #Deriv Calc
    # Rdot2b = - ((-2*b0 + (2*b0-2*b1)*x1 + (2*b1 - 2*b2)*x1**2 + (2*b2 - 2*b3)*x1**3 + (2*b3 - 2*b4)*x1**4 + 2*b4*x1**5) + (b0 + (b0 - b1)*x1 + (3*b1 - 3*b2)*x1**2 + (5*b2 - 5*b3)*x1**3 + (7*b3 - 7*b4)*x1**4 + 9*b4*x1**5)*math.log(x1)) / (2*x1**(3/2))
    #Symbolab
    #Rdot2bsym = -((((x1-1)*(b0 + b1*x1 + b2*x1**2 + b3*x1**3 + b4*x1**4))/(x1) + math.log(x1) * (5*b4*x1**4 + 4*b3*x1**3 - 4*b4*x1**3 + 3*b2*x1**2 - 3*b3*x1**2 + 2*b1*x1 - 2*b2*x1 + b0 - b1)) * (math.sqrt(x1)/x1) - ((b0 + b1*x1 + b2*x1**2 + b3*x1**3 + b4*x1**4)*math.log(x1)*(x1-1))/(2*x1**(3/2)))
    #print(Rdot2b)
    #print(Rdot2bsym)

    P = math.sqrt(1-x)*K
    Pdot = math.sqrt(1-x) * Kdot + -(1/(2*math.sqrt(1-x))) * K

    R1 = (1/x)*E
    R2 = ((x-1)/(math.sqrt(x)))*K
    R = R1 + R2
    Rdot1 = (1/x)*Edot + (-1/x**2)*E
    Rdot2 = ((x-1)/(math.sqrt(x)))*Kdot + ((x+1)/(2*x*math.sqrt(x)))*K
    Rdot = Rdot1 + Rdot2

    # acurr = 6.6105
    acurr = 1.0314

    lambacofa = (4/math.pi)*math.sqrt(acurr**3/1.0) * P
    lambacofi = (2/math.pi)*math.sqrt(acurr/1.0) * R
    dotlambacofa = (4/math.pi)*math.sqrt(acurr**3/1.0) * Pdot
    dotlambacofi = (2/math.pi)*math.sqrt(acurr/1.0) * Rdot

    solution = np.linalg.solve([[lambacofa, lambacofi], [dotlambacofa,dotlambacofi]], [-1, 0])
    print(solution)

    acurr = 1.0314
    lambacofi = -0.54
    bigX = (math.pi/(2*lambacofi)) * math.sqrt(1/acurr)
    bigZ = bigX**-2

    alpha0 = 0.0
    alpha1 =  2.467410607
    alpha2 = -1.907470562
    alpha3 = 35.892442177
    alpha4 = -214.67979624
    alpha5 = 947.773272608
    alpha6 = -2114.861134906
    alpha7 = 2271.240058672
    alpha8 = -1127.457440108
    alpha9 = 192.953875268
    alpha10 = 8.577733773
    alphaco = np.array([alpha0,alpha1,alpha2,alpha3,alpha4,alpha5,alpha6,alpha7,alpha8,alpha9,alpha10])

    beta0 = 1.0
    beta1 = 0.4609698838
    beta2 = 13.7756315324
    beta3 = -69.1245316678
    beta4 = 279.067183250
    beta5 = -397.6628952136
    beta6 = -70.0139935047
    beta7 = 528.0334266841
    beta8 = -324.9303836520
    beta9 = 20.5838245170
    beta10 = 18.8165370778
    betaco = np.array([beta0,beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8, beta9, beta10])

    cv = 0.0
    sumalpha =  0.0
    sumbeta = 0.0
    for i in range(11):
        sumalpha = sumalpha + (alphaco[i] * bigZ**i)
        sumbeta = sumbeta + (betaco[i] * bigZ**i)

    cv = sumalpha * sumbeta
    print(cv)
    
if __name__ == '__main__':
    import numpy as np
    #(-5.0,9.88),(-5.5,8.17),(-6.0,6.87),(-6.5,5.85),(-7.0,5.06),(-7.5,4.36),(-8.0,3.87),(-8.5,3.42),(-9.0,3.06),(-9.5,2.75),(-10.0,2.46),(-11.0,2.05),(-12.0,1.71),(-13.0,1.46),(-14.0,1.25),(-15.0,1.1),(-16.0,1.0)
    #(9.88,-5.0),(8.17,-5.5),(6.87,-6.0),(5.85,-6.5),(5.06,-7.0),(4.36,-7.5),(3.87,-8.0),(3.42,-8.5),(3.06,-9.0),(2.75,-9.5),(2.46,-10.0),(2.05,-11.0),(1.71,-12.0),(1.46,-13.0)(1.25,-14.0),(1.1,-15.0),(1.0,-16.0)
    xs = [-5.0,-5.5,-6.0,-6.5,-7.0,-7.5,-8.0,-8.5,-9.0,-9.5,-10.0,-11.0,-12.0,-13.0,-14.0,-15.0,-16.0]
    ys = [9.88,8.17,6.87,5.85,5.06,4.36,3.87,3.42,3.06,2.75,2.46,2.05,1.71,1.46,1.25,1.1,1.0]
    p = np.poly1d(np.polyfit(xs, ys, deg=4))
    print(p)




