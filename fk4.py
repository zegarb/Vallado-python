import numpy as np

# ----------------------------------------------------------------------------
#
#                           function fk4
#
#  this function converts vectors between the b1950 and j2000 epochs.  be
#    aware that this process is not exact. there are different secular rates
#    for each system, and there are differences in the central location. the
#    matrices are multiplied directly for speed.
#
#  author        : david vallado                  719-573-2600   21 jun 2002
#
#  revisions
#                -
#
#  inputs          description                    range / units
#    rj2000      - initial j2000 eci position vec er, km, etc
#    vj2000      - initial j2000 eci velocity vec
#    direction   - which set of vars to output    from  too
#
#  outputs       :
#    fk4m        - conversion b1950 to j2000
#    fk4mi       - conversion j2000 to b1950
#
#  locals        :
#
#  coupling      :
#
#  references    :
#    vallado       2001, 227-228
#
# [fk4m] = fk4 ( );
# rj2000 = fk4m * rb1950;
# ----------------------------------------------------------------------------
def fk4(option: str):
    match option:
        case 'a':
            # in book
            # 1950 - 2000
            return np.array([[0.9999256794956877, -0.0111814832204662, -0.0048590038153592],
                             [0.0111814832391717, 0.9999374848933135, -0.0000271625947142],
                             [0.0048590037723143, -0.0000271702937440, 0.9999881946043742]])

        case 'b':
            # stk approach
            # New way is formed by multiplying the matrices on pages
            # 173 and 174 and adding in the correction to equinox given
            # on page 168 of the supplement to the astronomical almanac
            # 1950 - 2000
            return np.array([[0.999925678612394, -0.011181874556714, -0.004858284812600],
                             [0.011181874524964, 0.999937480517880, -0.000027169816135],
                             [0.004858284884778, -0.000027156932874, 0.999988198095508]])

        case default:
            # from Exp supp to Ast Almanac pg 185 6x6
            # 1950 - 2000
            return np.array([[0.9999256782, -0.0111820611, -0.0048579477,
                              0.00000242395018, -0.00000002710663, -0.00000001177656],
                             [0.0111820610, 0.9999374784, -0.0000271765,
                              0.00000002710663, 0.00000242397878, -0.00000000006587],
                             [0.0048579479, -0.0000271474, 0.9999881997,
                              0.00000001177656, -0.00000000006582, 0.00000242410173],
                             [-0.000551, -0.238565, 0.435739,
                              0.99994704, -0.01118251, -0.00485767],
                             [0.238514, -0.002667, -0.008541,
                              0.01118251, 0.99995883, -0.00002718],
                             [-0.435623, 0.012254, 0.002117,
                              0.00485767, -0.00002714, 1.00000956]])

def fk4i(option: str):
    match option:
        case 'a':
            # in book
            # 1950 - 2000
            return np.array([[9.99925679e-01, 1.11814832e-02, 4.85900377e-03],
                             [-1.11814832e-02, 9.99937485e-01, -2.71702937e-05],
                             [-4.85900382e-03, -2.71625947e-05, 9.99988195e-01]])
        case 'b':
            # stk approach
            # New way is formed by multiplying the matrices on pages
            # 173 and 174 and adding in the correction to equinox given
            # on page 168 of the supplement to the astronomical almanac
            # 1950 - 2000
            return np.array([[9.99925679e-01, 1.11818745e-02, 4.85828488e-03],
                             [-1.11818746e-02,  9.99937481e-01, -2.71569329e-05],
                             [-4.85828481e-03, -2.71698161e-05,  9.99988198e-01]])
        case default:
            # from Exp supp to Ast Almanac pg 185 6x6
            # 1950 - 2000
            return np.array([[9.99925678e-01, 1.11820610e-02, 4.85794790e-03,
                              -5.51000000e-04, 2.38514000e-01, -4.35623000e-01],
                             [-1.11820611e-02, 9.99937478e-01, -2.71474000e-05,
                              -2.38565000e-01, -2.66700000e-03, 1.22540000e-02],
                             [-4.85794770e-03, -2.71765000e-05, 9.99988200e-01,
                              4.35739000e-01, -8.54100000e-03, 2.11700000e-03],
                             [2.42395018e-06, 2.71066300e-08, 1.17765600e-08,
                              9.99947040e-01, 1.11825100e-02, 4.85767000e-03],
                             [-2.71066300e-08, 2.42397878e-06, -6.58200000e-11,
                              -1.11825100e-02, 9.99958830e-01, -2.71400000e-05],
                             [-1.17765600e-08, -6.58700000e-11,  2.42410173e-06,
                              -4.85767000e-03, -2.71800000e-05, 1.00000956e+00]])
