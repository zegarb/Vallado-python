#     -----------------------------------------------------------------
#
#                              Ex10_1
#
#  this file demonstrates example 10-1.
#
#                          companion code for
#             fundamentals of astrodyanmics and applications
#                                 2013
#                            by david vallado
#
#     (h)               email davallado@gmail.com
#     (w) 719-573-2600, email dvallado@agi.com
#
#     *****************************************************************
#
#  current :
#            19 feb 19  david vallado
#                         update for new constants
#  changes :
#            13 feb 07  david vallado
#                         original baseline
#
#     *****************************************************************

# problem 1
import numpy as np
print('problem 1 --------------------------\n')
xo = np.array([1, 2, 3, 4, 5, 6, 7, 8])
yo = np.array([1, 1, 2, 3, 3, 4, 4, 6])
print(f'xo: {xo}\n')
print(f'yo: {yo}\n')
ata = np.array([[8, sum(xo)], [sum(xo), sum(xo * xo)]])
atb = np.array([[sum(yo)], [sum(xo * yo)]])
atai = np.linalg.inv(ata)
ans = atai @ atb
print(f'ata = \n{ata}\n atb = \n{atb}\n atai = \n{atai}\n ans = \n{ans}')