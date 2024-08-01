import numpy as np
import math

# from the book page 798: 4th edition, or 814: 5th edition

# Note: numpy does NOT like 'vertical slice' arrays (shape(x, 1)).
# it will automatically convert a vertical slice to a 1D array
# (s[0:3, i], pts[0:6, i]), and will expect a 1D array when
# assigning to a vertical slice (sigmapts[0:3, jj], y[0:6, i]).
# Trying to assign a shape(3, 1) array to sigmapts[0:3, jj] will crash! -zeg

def unscentedkalman(reci: np.ndarray, veci: np.ndarray, covariance: np.ndarray):
    s = math.sqrt(6) * np.linalg.cholesky(covariance)
    sigmapts = np.zeros([6, 12])
    #reci, veci need to be 1D
    reci = reci.flatten()
    veci = veci.flatten()
    for i in range(6):
        # jj modified for starting at index 0 instead of 1
        jj = i * 2
        sigmapts[0:3, jj] = reci + s[0:3, i]
        sigmapts[3:6, jj] = veci + s[3:6, i]
        sigmapts[0:3, jj+1] = reci - s[0:3, i]
        sigmapts[3:6, jj+1] = veci - s[3:6, i]

    return sigmapts

def unscentedreassemble(pts: np.ndarray):
    yu = np.zeros(6)
    for i in range(12):
        yu = yu + pts[0:6, i]
    yu = yu / 12

    y = np.zeros([6, 12])
    for i in range(12):
        y[0:6, i] = pts[0:6, i] - yu
    tmp = (y @ y.T) / 12
    cov = (tmp + tmp.T) / 2
    return cov