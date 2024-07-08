import os
import numpy as np
import spacetime_utils as stu
# -----------------------------------------------------------------------------
#
#                           function initjplde
#
#  this function initializes the jpl planetary ephemeris data. the input
#  data files are from processing the ascii files into a text file of sun
#  and moon positions.
#
#  author        : david vallado                  719-573-2600   22 jan 2018
#
#  revisions
#
#  inputs          description                    range / units
#    infilename  - input data file
#
#  outputs       :
#    jpldearr    - array of jplde data records
#    jdjpldestart- julian date of the start of the jpldearr data
#
#  locals        :
#                -
#
#  coupling      :
#
#  references    :
#  [jpldearr, jdjpldestart, jdjpldestartFrac] = initjplde('D:\Codes\LIBRARY\DataLib\sunmooneph_430t.txt');
#  -------------------------------------------------------------------------- */


def initjplde(infilename: str):
    """this function initializes the jpl planetary ephemeris data. the input
       data files are from processing the ascii files into a text file of sun
       and moon positions.

    Parameters
    ----------
    infilename : str
        name of the data file to use

    Returns
    -------
    jpldearr : dictionary of arrays
        arrays of jplde data records
    jdjpldestart : float
        julian date of the start of the jpldearr data
    """
    #double jdtdb, jdtdbf;
    #string pattern;
    #Int32 i;
    #jdjpldestart = 0.0;
    #jdjpldestartFrac = 0.0;

    # read the whole file at once into lines of an array
    filename = os.path.join(os.path.dirname(__file__), 'data', infilename)

    filedat = np.loadtxt(filename)
    filedim = filedat.shape[0]

    jpldearr = {}
    #load data into x y z arrays
    jpldearr['year'] = filedat[:, 0]
    jpldearr['mon'] = filedat[:, 1]
    jpldearr['day'] = filedat[:, 2]
    jpldearr['rsun1'] = filedat[:, 3]
    jpldearr['rsun2'] = filedat[:, 4]
    jpldearr['rsun3'] = filedat[:, 5]
    jpldearr['rsmag'] = filedat[:, 6]
    jpldearr['rmoon1'] = filedat[:, 8]
    jpldearr['rmoon2'] = filedat[:, 9]
    jpldearr['rmoon3'] = filedat[:, 10]
    for i in range(filedim):
        jdtdb,jdtdbf = stu.jday(jpldearr['year'][i], jpldearr['mon'][i],
                                jpldearr['day'][i], 0, 0, 0.0)
        jpldearr['mjd'][i] = jdtdb + jdtdbf - 2400000.5
    # ---- find epoch date
    jdjpldestart,jdjpldestartFrac = stu.jday(jpldearr['year'][0],
                                             jpldearr['mon'][0],
                                             jpldearr['day'][0], 0, 0, 0.0)
    # doesn't have the -2400000.5 for some reason? -zeg
    # jpldearr['mjd'][0] = jdjpldestart + jdjpldestartFrac
    return jpldearr, jdjpldestart, jdjpldestartFrac