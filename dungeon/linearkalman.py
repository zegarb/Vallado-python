import numpy as np
import space_conversions as sc
import spacemath_utils as smu
import orbit_utils as obu
from space_constants import *

# Algorithm 69 implementation. shelving it until there's some way to test it
# a few errors with matrix size mismatches - zeg
def linearkalman(xnom: np.ndarray, pnom: np.ndarray,
                 Q: np.ndarray, R: np.ndarray, obsarr: np.ndarray,
                 sitelatgd: float, sitelon: float, sitealt: float, ttt: float,
                 xnomjd: float, lod: float, xp: float, yp: float, ddpsi: float,
                 ddeps: float):

    sinlat, coslat, sinlon, coslon = smu.getsincos(sitelatgd, sitelon)
    ecef2sez = np.array([[sinlat*coslon, sinlat*sinlon, -coslat],
                         [-sinlon, coslon, 0],
                         [coslat*coslon, coslat*sinlon, sinlat]])
    recef, _, _ = sc.eci2ecef(xnom[[0, 1, 2]], xnom[[3, 4, 5]],
                              np.zeros(3), ttt, xnomjd, lod, xp, yp, 'c',
                              ddpsi, ddeps)
    rmag = smu.mag(recef)
    F = np.zeros([6, 6])
    F[0, 3] = 1
    F[1, 4] = 1
    F[2, 5] = 1
    F[3, 0] = (-mu/rmag**3) + (3 * mu * recef[0]**2) / rmag**5
    F[3, 1] = 3 * mu * recef[0] * recef[1] / rmag**5
    F[3, 2] = 3 * mu * recef[0] * recef[2] / rmag**5
    F[4, 0] = F[3, 1]
    F[4, 1] = (-mu/rmag**3) + (3 * mu * recef[1]**2) / rmag**5
    F[4, 2] = 3 * mu * recef[1] * recef[2] / rmag**5
    F[5, 0] = F[3, 2]
    F[5, 1] = F[4, 2]
    F[5, 2] = (-mu/rmag**3) + (3 * mu * recef[2]**2) / rmag**5

    H = []
    dxbar = []
    dxhat = []
    xbar = []
    xhat = []
    pbar = []
    phat = []
    z = []
    b = []
    K = []

    i = 0
    for obs in obsarr:
        rhosez, _ = sc.razel2sez(obs['rng'], obs['az'], obs['el'], 0, 0, 0)
        rho = smu.mag(rhosez)
        temp = np.array([[rhosez[0]/rho, rhosez[1]/rho, rhosez[2]/rho],
                        [0, 0, 0],
                        [0, 0, 0]])

        temp[1, 0] = rhosez[1] / ((rhosez[1]**2 / rhosez[0]**2 + 1) * rhosez[0]**2)
        temp[1, 1] = -1 / ((rhosez[1]**2 / rhosez[0]**2 + 1) * rhosez[0])
        # temp[1, 2] = 0
        temp[2, 0] = -(rhosez[0] * rhosez[2]) / (rho**3
                                        * math.sqrt(-rhosez[2]**2 / rho**2 + 1))
        temp[2, 1] = -(rhosez[1] * rhosez[2]) / (rho**3
                                        * math.sqrt(-rhosez[2]**2 / rho**2 + 1))
        temp[2, 2] = (rho**2 - rhosez[2]**2) / (rho**3
                                            * math.sqrt(-rhosez[2]**2 / rho**2 + 1))

        H.append(np.cross(temp, ecef2sez))
        time = obs['time'] + obs['timef']
        dt = (time - xnomjd) * 86400
        stm = np.eye([6, 6]) + F
        temp1 = 10
        j = 1
        while temp1 > 1e-12:
            temp1 = dt**j / math.factorial(j)
            stm = stm + F**j * temp1

        rest = xnom[[0, 1 ,2]]
        vest = xnom[[3, 4, 5]]
        rbar, vbar = obu.pkepler(rest, vest, dt, 0, 0)
        xbar.append(np.append(rbar, vbar))

        z.append(np.array([obs['rng'], obs['az'], obs['el']]))
        #initial xhat = 0
        if i == 0:
            xbar.append(stm @ np.zeros((6, 1)))
            pbar.append(stm @ pnom @ stm.T + Q)
        else:
            xbar.append(stm @ xhat[i - 1])
            pbar.append(stm @ phat[i - 1] @ stm.T + Q)


        xrho, xaz, xel, _, _, _ = sc.rv2razel(xbar[[0, 1, 2]], xbar[[3, 4, 5]],
                                sitelatgd, sitelon, sitealt, ttt, time, lod,
                                xp, yp, 0, ddpsi, ddeps)
        b.append(z - H[i] @ np.array([xrho, xaz, xel]))
        # Pretty sure H is 3x3 and pbar is 6x6. -zeg
        K.append(pbar[i] @ H[i].T @ (H[i] @ pbar[i] @ H[i].T + R)**-1)
        dxhat.append(dxbar[i] + K[i] @ (b[i] - H[i] @ dxbar[i]))
        phat.append(pbar[i] - K[i] @ H[i] @ pbar[i])
        xhat.append(xbar[i] + dxhat[i])
        i = i + 1

    return xhat, phat

