#     -----------------------------------------------------------------
#
#                              Ex10_2
#
#  this file demonstrates example 10-2.
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

# problem 2
import numpy as np
print('example 10-2 -------------------------\n')
xo = np.array([1, 2, 3, 4, 5, 6, 7, 8])
yo = np.array([1, 1, 2, 3, 3, 4, 7, 6])
print(f'xo: {xo}\n')
print(f'yo: {yo}:\n')
ata = np.array([[8, sum(xo)], [sum(xo), sum(xo * xo)]])
atb = np.array([[sum(yo)], [sum(xo * yo)]])
atai = np.linalg.inv(ata)
ans = atai @ atb
print(f'ata = \n{ata}\n'
      f'atb = \n{atb}\n'
      f'atai = \n{atai}\n'
      f'ans=\n{ans}')
yc = ans[0] + ans[1] * xo
print('yc:\n', (yc))
res = yo - yc
print('res:\n', (res))
rms = np.sqrt(1 / 8 * sum(res * res))
conv = np.sqrt(sum(res * res) / 7)
print('sum res 2 %11.7f  rms %11.7f  sum res /7 %11.7f \n'
      % (sum(res * res), rms, conv))
print('adjusted cov %11.7f  %11.7f \n'
      % (conv * np.sqrt(atai[0, 0]), conv * np.sqrt(atai[1, 1])))

print('Example 10-2 with observation 7 removed -------------------------\n' % ())
#clear('all')
xo = np.array([1, 2, 3, 4, 5, 6, 8])
yo = np.array([1, 1, 2, 3, 3, 4, 6])
print('xo:\n', (xo))
print('yo:\n', (yo))
ata = np.array([[7, sum(xo)], [sum(xo), sum(xo * xo)]])
atb = np.array([[sum(yo)], [sum(xo * yo)]])
atai = np.linalg.inv(ata)
ans = atai @ atb
print(f'ata = \n{ata}\n'
      f'atb = \n{atb}\n'
      f'atai = \n{atai}\n'
      f'ans=\n{ans}')
yc = ans[0] + ans[1] * xo
res = yo - yc
rms = np.sqrt(1 / 7 * sum(res * res))
conv = np.sqrt(sum(res * res) / 6)
print('sum res 2 %11.7f  rms %11.7f  sum res /7 %11.7f \n'
      % (sum(res * res), rms, conv))
print('adjusted cov %11.7f  %11.7f \n'
      % (conv * np.sqrt(atai[0, 0]), conv * np.sqrt(atai[1, 1])))

# problem 3 example with weighting
print('Section 10.2.2 Linear Weighted Least Squares ---------------------\n')
xo = np.array([1, 2, 3, 4, 5, 6, 7, 8])
yo = np.array([1, 1, 2, 3, 3, 4, 4, 6])
print('xo:\n', (xo))
print('yo:\n', (yo))
w1 = 1.0 / 0.1

w2 = 1.0 / 0.02
w = np.zeros((8, 8))
w[0, 0] = w1
w[1, 1] = w2
w[2, 2] = w1
w[3, 3] = w2
w[4, 4] = w1
w[5, 5] = w2
w[6, 6] = w1
w[7, 7] = w2
i = np.array([1, 1, 1, 1, 1, 1, 1, 1])
a = np.array([i, xo]).T
b = np.array([yo]).T
atw = a.T @ w
atwa = atw @ a
atwb = atw @ b
atwai = np.linalg.inv(atwa)
ans = atwai @ atwb
print(f'atwa = \n{atwa}\n'
      f'atwb = \n{atwb}\n'
      f'atwai = \n{atwai}\n'
      f'ans=\n{ans}')
#    ata=[8 sum(xo); sum(xo) sum(xo*xo')]
#    atb=[sum(yo); sum(xo*yo')]
#    atai = np.linalg.inv(ata)
#    ans= atai*atb

yc = ans[0] + ans[1] * xo
print('yc:\n', (yc))
res = yo - yc
print('res:\n', (res))
rms = np.sqrt(1 / 8 * sum(res * res))
conv = np.sqrt(sum(res * res) / 7)
print('sum res 2 %11.7f  rms %11.7f  sum res /7 %11.7f \n'
      % (sum(res * res), rms, conv))
print('adjusted cov %11.7f  %11.7f \n'
      % (np.sqrt(atwai[0, 0]), np.sqrt(atwai[1, 1])))

# # example
# xo = np.array([1, 2, 3, 4, 5, 6, 7, 8])
# yo = np.array([[1, 1, 2, 3, 3, 4, 4, 6]]).T
# i = np.array([1, 1, 1, 1, 1, 1, 1, 1])
# a = np.array([i, xo]).T
# b = np.array(yo)
# w1 = 1.0 / 0.1

# w2 = 1.0 / 0.02

# #  w1 = 1.0;
# #  w2 = 1.0;

# w[0, 0] = w1
# w[1, 1] = w2
# w[2, 2] = w1
# w[3, 3] = w2
# w[4, 4] = w1
# w[5, 5] = w2
# w[6, 6] = w1
# w[7, 7] = w2
# a.T @ w @ a
# covinv = np.linalg.inv(a.T @ w @ a)
# a.T @ b
# np.linalg.inv(a.T @ w @ a)
# np.linalg.inv(a.T @ w @ a) @ a.T @ b
# np.sqrt(covinv[0, 0])
# np.sqrt(covinv[1, 1])
