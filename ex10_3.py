#     -----------------------------------------------------------------
#
#                              Ex10_3
#
#  this file demonstrates example 10-3.
#
#                          companion code for
#             fundamentals of astrodyanmics and applications
#                                 2007
#                            by david vallado
#
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#             7 jun 07  david vallado
#                         original
#  changes :
#            13 feb 07  david vallado
#                         original baseline
#
#     *****************************************************************

# problem 3
import numpy as np
print('Example 10-3 ---------')
xo = np.array([1.0, 2.0, 3.0, 4.0])
yo = np.array([2.5, 8.0, 19.0, 50.0])
print('xo:\n', (xo))
print('yo:\n', (yo))
beta = np.log(yo[3] / yo[2]) / np.log(4 / 3)
alpha = yo[2] / 3 ** beta
print(f'alpha = {alpha}, beta = {beta}')

# do first time manually to get intermediate values
yn = np.zeros(4)
parynalp = np.zeros(4)
parynbet = np.zeros(4)
for i in range(4):
    yn[i] = alpha * xo[i] ** beta
    parynalp[i] = xo[i] ** beta
    parynbet[i] = alpha * np.log(xo[i]) * xo[i] ** beta

ata = np.array([[np.sum(parynalp**2), np.sum(parynalp*parynbet)],
                [np.sum(parynalp*parynbet), np.sum(parynbet**2)]])
atb = np.array([[np.sum(parynalp * (yo - yn))], [np.sum(parynbet * (yo - yn))]])
atai = np.linalg.inv(ata)
ans = atai @ atb
alpha = alpha + ans[0, 0]
beta = beta + ans[1, 0]
b = np.empty(4)
for j in range(4):
    b[j] = yo[j] - alpha * (xo[j] ** beta)

print('yn = \n', yn)
print('parynalp = \n', parynalp)
print('parynbet = \n', parynbet)
print('ata = \n', ata)
print('atb = \n', atb)
print('atai =\n', atai)
print('b = \n', b)

rms = np.sqrt(np.sum(b ** 2) / 4)
rmsold = rms
print('0 dx %11.7f  %11.7f ans %11.7f  %11.7f rms %11.7f \n'
      % (ans[0, 0], ans[1, 0], alpha, beta, rms))
for loop in range(3):
    for i in range(4):
        yn[i] = alpha * xo[i] ** beta
        parynalp[i] = xo[i] ** beta
        parynbet[i] = alpha * np.log(xo[i]) * xo[i] ** beta
    ata = np.array([[np.sum(parynalp**2), np.sum(parynalp*parynbet)],
                [np.sum(parynalp*parynbet), np.sum(parynbet**2)]])
    atb = np.array([[np.sum(parynalp * (yo - yn))], [np.sum(parynbet * (yo - yn))]])
    atai = np.linalg.inv(ata)
    ans = atai @ atb
    alpha = alpha + ans[0, 0]
    beta = beta + ans[1, 0]
    for j in range(4):
        b[j] = yo[j] - alpha * (xo[j] ** beta)
    print(f' b = \n{b}')
    rms = np.sqrt(np.sum(b ** 2) / 4)
    rmsdel = (rmsold - rms) / rmsold
    print('%2i dx %11.7f  %11.7f ans %11.7f  %11.7f rms %11.7f %11.7f \n'
          % (loop + 1, ans[0, 0], ans[1, 0], alpha, beta, rms, rmsdel))
    rmsold = rms
