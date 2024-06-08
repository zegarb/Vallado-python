#
# calculate algorithm 34 quantities
#

import numpy as np
import spacemath_utils as smu
import spacetime_utils as stu
import orbit_utils as obu
from space_constants import sunradius, re, au

#au is mean earth distance from the sun
angumb = np.arctan((sunradius - re) / au)
angpen = np.arctan((sunradius + re) / au)


reci1 = np.array([- 41221.79149309,8864.59854079,0.0])
veci1 = np.array([- 0.646416796,- 3.005940793,- 0.0])
# +50 and -80 seem to work here to get the proper angles
dtsec = 21000
reci,veci,error = obu.kepler(reci1,veci1,dtsec)
year = 2008
mon = 3
day = 16
hr = 6
min = 13
sec = 0.0
jd,jdfrac = stu.jday(year,mon,day,hr,min,sec)
rsun,rtasc,decl = obu.sun(jd + jdfrac)


# now for algorithm 34
pen,umb = obu.shadow(reci, rsun, angumb, angpen)

print('U %r, P %r' % (pen,umb))
#s = 2.0 * smu.mag(reci) * ang1
#print(' %11.7f  %11.7f \n' % (s,(s / smu.mag(veci)) / 60.0))

